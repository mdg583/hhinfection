#include <Rcpp.h>
#include <math.h>
#include "lik_obs.h"
#include "utilities.h"
#include "sec_genfun.h"

//' Log-probabilities for observed infections given number of observed and secondary cases
//'
//' This is the same as function h(m=m,i=obs,j=act), with two differences:
//' 1. We don't actually sum the values, but instead fill values in an array, to group summation together later
//' 2. The arguement lfact contains (log) factors to multiply against each element, from further up summation.
//' This is to allow a single summation, which uses the log-sum-exp method
//'
//' @param pbuf memory space into which min(act,obs)-max(0,obs+act-m) + 1 values will be stored
//' @param m number of potential secondary cases
//' @param act number of secondary cases among potential secondary cases
//' @param obs number of observed secondary cases among potential secondary cases
//' @param lfact log-factor which will be added to each probability
//' @l_Sn log of test sensitivity
//' @l_Sp log of test specificity
//' @l_1mSn log of one minus test sensitivity
//' @l_1mSp log of one minus test specificity
//' @return
//' Values will be written to pbuf[0]...pbuf[min(act,obs)-max(0,obs+act-m)]. The value returned
//' is the maximum value among those written to pbuf
double lprob_obs_sec(double *pbuf, int m, int act, int obs, double lfact, double l_Sn, double l_Sp, double l_1mSn, double l_1mSp){
  int astart = std::max(0,obs+act-m);
  int aend = std::min(act,obs);
  
  // Performance analysis using RcppClock
  //   ticker  time   prop
  //   coef    74.98  52.6%
  //   loop    67.35  47.3%
  // Performance comments: binomial coefs are O(n), loop is O(n)
  long long coef1 = choose(obs,astart);
  long long coef2 = choose(m-obs,act-astart);
  double maxpbuf = 0.0;
  for(int i = 0; i < aend-astart+1; i++){
    int a = astart + i;
    int b = act-a;
    int c = obs-a;
    int d = m-a-b-c;
    if(i>0){
      coef1 = choose_rinc(coef1, obs, a);
      coef2 = choose_rdec(coef2, m-obs, act-a);
    }
    pbuf[i] = lfact + log(coef1) + log(coef2) + mul(a,l_Sn) + mul(b,l_1mSn) + mul(c,l_1mSp) + mul(d,l_Sp);
    if(i==0 || pbuf[i] > maxpbuf) maxpbuf = pbuf[i];
  }
  return maxpbuf;
}

//' Log-probability for observed secondary infections given number of primary infections
//'
//' @param obs Number of observed infections among potential secondary cases
//' @param mis Number of missed observations among potential secondary cases
//' @param pri Number of primary infections
//' @param m Number of potential secondary cases
//' @param gdat pre-computed log-probabilities for number of secondary infections
//' @l_Sn log of test sensitivity
//' @l_Sp log of test specificity
//' @l_1mSn log of one minus test sensitivity
//' @l_1mSp log of one minus test specificity
//' @return Log-probability for observed secondary infections
double lprob_sec_obs_pri(int obs, int mis, int pri, int m, Gdat gdat,
                         double l_Sn, double l_Sp, double l_1mSn, double l_1mSp){
  double *r = g(gdat,pri,m); // r[0...m] are log-probabilities for j=0...m secondary infections

  // Performance analysis using RcppClock
  //   ticker  time    prop
  //   all     2439.2  100.0%
  //   loop    2160.8   88.6%
  //   mse      125.3    5.1%
  // Peformance comments: Cost is in the loop. Cost of mse is not too much, nor malloc (~6.3%).

  // How much memory is needed? Loop over summations to count memory needed.
  int buf_size = 0;
  for(int j = 0; j <= m; j++){
    for(int e = std::max(0,j-(m-mis)); e <= std::min(mis,j); e++){
      buf_size += std::min(j-e,obs) - std::max(0,obs+j-e-(m-mis)) + 1;
    }
  }
  double *pbuf = (double *)malloc(buf_size*sizeof(double));
  int n = 0; // index marker for pbuf
  double maxpbuf = 0.0;
  
  // number of secondary infections
  for(int j = 0; j <= m; j++){
    int estart = std::max(j-(m-mis),0);
    int eend = std::min(mis,j);
    long long coef = choose(mis,estart);
    double lmj = log(choose(m,j));
    // positives among missed observations
    for(int i = 0; i < eend-estart+1; i++){
      int e = estart + i;
      if(i > 0) coef = choose_rinc(coef, mis, e);

      // probability of observed secondary cases
      double maxpbufi = lprob_obs_sec(pbuf+n, m-mis, j-e, obs, r[j] - lmj + log(coef), l_Sn, l_Sp, l_1mSn, l_1mSp);

      if(n==0 || maxpbufi > maxpbuf) maxpbuf = maxpbufi;
      n += std::min(j-e,obs)-std::max(0,obs+j-e-(m-mis))+1;
    }
  }
  // shortcut: probability 0
  if(maxpbuf==-1.0/0.0){
    free(pbuf);
    return maxpbuf;
  }
  // log-sum-exp: subtract maximum value, sum exponents, take log, add maximum value
  double prob = 0.0;
  for(int i = 0; i < n; i++){
    prob += exp(pbuf[i] - maxpbuf);
  }
  free(pbuf);
  return log(prob) + maxpbuf;
}

//' Log-probability for observed infections given primary infections
//'
//' @param obs Array of observations. -1=missed, 0=observed negative, 1=observed positive
//' @param pri Array of primary infections (0 or 1)
//' @param hh_size Size of household
//' @param gdat pre-computed log-probabilities for secondary infection numbers
//' @param l_Sn log of test sensitivity
//' @param l_Sp log of test specificity
//' @param l_1mSn log of one minus test sensitivity
//' @param l_1mSp log of one minus test specificity
//' @return log-probability
double lprob_obs_pri(int* obs, int* pri, const int hh_size, Gdat gdat,
                      const double l_Sn, const double l_Sp, const double l_1mSn, const double l_1mSp){
  double pprim; // probability of primary infection observation
  // If secondary infection rate (q) is 0, then only primary infections are infected. So gat
  // entire household observation likelihood
  if(gdat.q == 0.0){
    int a=0,b=0,c=0,d=0;
    for(int i = 0; i < hh_size; i++){
      if(obs[i]!=-1){ // if individual is not a missed observation
        a +=  obs[i] &&  pri[i];
        b += !obs[i] &&  pri[i];
        c +=  obs[i] && !pri[i];
        d += !obs[i] && !pri[i];
      }
    }
    pprim = mul(a,l_Sn) + mul(b,l_1mSn) + mul(c,l_1mSp) + mul(d,l_Sp);
    return pprim;
  }
  // observations among primary infections
  int a=0,b=0;
  for(int i = 0; i < hh_size; i++){
    if(obs[i]!=-1){ // if individual is not a missed observation
      a +=  obs[i] && pri[i];
      b += !obs[i] && pri[i];
    }
  }
  pprim = mul(a,l_Sn) + mul(b,l_1mSn);
  if(pprim == -1.0 / 0.0) return pprim; // shortcut, prob of primary observations is 0

  // Counts needed for secondary observation probability
  int k=0,mis=0,obs2=0,m=0;
  for(int i = 0; i < hh_size; i++){
    if(pri[i]){
      k++; // primary cases
    }else{
      m++; // potential secondary case
      mis  += (obs[i]==-1); // missed observations among potential secondary cases
      obs2 += (obs[i]==1);  // observed cases among potential secondary cases
    }
  }
  double psec = lprob_sec_obs_pri(obs2, mis, k, m, gdat, l_Sn, l_Sp, l_1mSn, l_1mSp);
  return pprim + psec;
}

//' Probability of primary infection given independent individual probabilities
//'
//' @param pri Array of primary infections (0 or 1)
//' @param lprim log-probabilities of primary infection
//' @param l_1mprim log of one minus probabilities of primary infection
//' @param n household size
//' @return log-probability
double lprob_pri(int *pri, double *lprim, double *l_1mprim, int n){
  double lprob = 0;
  for(int i = 0; i < n; i++){
    lprob += mul(pri[i],lprim[i]) + mul(1-pri[i],l_1mprim[i]);
  }
  return lprob;
}

//' Draw a sample for primary infections
//'
//' @param pri memory into which boolean primary infection will be written
//' @param p probabilities of primary infection
//' @param n household size
void pri_sample(int *pri, double *p, int n){
  double *draws = (Rcpp::runif(n)).begin();
  for(int i = 0; i < n; i++){
    pri[i] = (draws[i] < p[i]);
  }
}

//' Log-probability for observations in a household
//'
//' T < 0: never simulate
//' T = 0: always simulate
//' T > 0: simulate for households over size T
//'
//' @param obs Array of observations. -1=missed, 0=observed negative, 1=observed positive
//' @param prim Array of primary infection probabilities
//' @param hh_size Size of household
//' @param q Household secondary infection rate parameter
//' @param Sn test sensitivity
//' @param Sp test specificity
//' @param T Household size threshold for primary infection simulations.
//' @param M Number of primary infection simulations
//' @return log-probability
double lprob_obs(int* obs, double *prim, const int hh_size, double q,
                const double Sn, const double Sp, const int T, int M){
  // Performance analysis using RcppClock
  //   ticker    time    prop
  //   all      145.8  100.0%
  //   lgdat      7.5    5.1%
  //   sim_lse    0.4    0.3%
  // Peformance comments: lgdat is relatively cheap. The cost of having independent q per household
  // is not much. Cost of log-sum-exponent is very cheap. Most of the cost is in the loop.

  // generate log-probabilities for secondary infection
  Gdat gdat = lgdat(hh_size,q);

  // allocate for primary infection array, and log of primary inf prob
  int *pri = (int *)malloc((hh_size)*sizeof(int));
  memset(pri,0,hh_size*sizeof(int));

  // Either iterate or simulate
  if(T < 0 || hh_size < T){
  	// Iterate
    // Prepare log-probabilities of primary infection
    double *lprim = (double *)malloc((hh_size)*sizeof(double));
    double *l_1mprim = (double *)malloc((hh_size)*sizeof(double));
    for(int i = 0; i < hh_size; i++){
      lprim[i] = log(prim[i]);
      l_1mprim[i] = log(1-prim[i]);
    }

    // prepare memory for probabilities
    double *pbuf = (double *)malloc(((int64_t)1<<hh_size)*sizeof(double));
    double max_prob = 0;
    // iterate over all possible primary infections, using binary integer, bits correspond to members
    for(int64_t priv = 0; priv < ((int64_t)1<<hh_size); priv++){
    	// Convert binary integer to boolean array
      for(int i = 0; i < hh_size; i++){
        pri[i] = (priv >> i) & 1;
      }
      // probability of primary infection
      double pprim = lprob_pri(pri, lprim, l_1mprim, hh_size);
      // probability of observation given primary infections
      pbuf[priv] = pprim + lprob_obs_pri(obs, pri, hh_size, gdat, log(Sn), log(Sp), log(1-Sn), log(1-Sp));
      if(pri==0 || pbuf[priv] > max_prob) max_prob = pbuf[priv];
    }
    free(pri);
    gfree(gdat);
    free(lprim);
    free(l_1mprim);

    // Shortcut: prob is 0
    if(max_prob == -1.0/0.0){
      free(pbuf);
      return max_prob;
    }
    // log-sum-exp
    double prob = 0;
    for(int64_t priv = 0; priv < ((int64_t)1<<hh_size); priv++){
      prob += exp(pbuf[priv] - max_prob);
    }
    free(pbuf);
    return log(prob) + max_prob;
  }else{
  	// Simulate
  	// prepare memory for probabilities
    double *pbuf = (double *)malloc(M*sizeof(double));
    double max_prob = 0;
    for(int i = 0; i < M; i++){
      // generate a primary infection sample into pri
      pri_sample(pri, prim, hh_size);
      // probability of observation given primary infection
      pbuf[i] = lprob_obs_pri(obs, pri, hh_size, gdat, log(Sn), log(Sp), log(1-Sn), log(1-Sp));
      if(pri==0 || pbuf[i] > max_prob) max_prob = pbuf[i];
    }
    free(pri);
    gfree(gdat);

    // shortcut: prob 0
    if(max_prob == -1.0/0.0){
      free(pbuf);
      return max_prob;
    }
    // log-sum-exp
    double prob = 0;
    for(int i = 0; i < M; i++){
      prob += exp(pbuf[i] - max_prob);
    }
    free(pbuf);
    return log(prob) + max_prob - log(M);
  }
}

//' Log-probability for observations for multiple households
//'
//' T < 0: never simulate
//' T = 0: always simulate
//' T > 0: simulate for households over size T
//'
//' @param obs Array of size sum(hh_sizes) of observations. -1=missed, 0=observed negative, 1=observed positive
//' @param prim Array of size sum(hh_sizes) of primary infection probabilities
//' @param hh_size Array of size n of household sizes
//' @param q Array of size n of secondary infection rates per household
//' @param n Number of households
//' @param Sn test sensitivity
//' @param Sp test specificity
//' @param T Household size threshold for primary infection simulations.
//' @param M Number of primary infection simulations
//' @return log-probability
double lprob_obs_hhs(int *obs, double *prim, int *hh_sizes, double *q, const int n,
                     const double Sn, const double Sp, const int T, const int M){
  double lprob = 0.0;
  int i = 0;
  for(int hh = 0; hh < n; hh++){
    lprob += lprob_obs(obs+i, prim+i, hh_sizes[hh], q[hh], Sn, Sp, T, M);
    i += hh_sizes[hh];
  }
  return lprob;
}
