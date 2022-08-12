// ## [[Rcpp::depends(RcppClock)]]
// #include <RcppClock.h>
#include <Rcpp.h>
#include <math.h>
#include "lik_obs.h"
#include "utilities.h"
#include "sec_genfun.h"

// Performance Analysis: uncomment to get performance (see RcppClock documentation)
// Rcpp::Clock *clock1;
// Rcpp::Clock *clock2;
// Rcpp::Clock *clock3;

//' Fill an array of log-probabilities for observed infections given the number of observed and secondary cases
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
//' @return maximum value among the elements filled-in to pbuf
double lprob_obs_sec(double *pbuf, int m, int act, int obs, double lfact, double l_Sn, double l_Sp, double l_1mSn, double l_1mSp){
  int astart = std::max(0,obs+act-m);
  int aend = std::min(act,obs);

  //Rprintf("m=%d,j=%d,i=%d  lfact=%.4f\n",m,act,obs,lfact);
  //Rprintf("a=%d ... %d\n", astart,aend);

  // Performance analysis using RcppClock
  //   ticker  time   prop
  //   coef    74.98  52.6%
  //   loop    67.35  47.3%
  // Performance comments: unsuprising as coef is O(n) as is the loop, since all ops in the loop are O(1)

  //clock3->tick("coef");
	long long coef1 = choose(obs,astart);
	long long coef2 = choose(m-obs,act-astart);
	double maxpbuf = 0.0;
	//clock3->tock("coef");
	//clock3->tick("loop");
	//Rprintf("[ a, b, c, d] c1 c2 = h\n");
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
		//Rprintf("[%2d,%2d,%2d,%2d] %2d %2d = %0.3f\n", a, b, c, d, coef1, coef2, exp(pbuf[i] - lfact));
		if(i==0 || pbuf[i] > maxpbuf) maxpbuf = pbuf[i];
	}
	//Rprintf("\n");
	//clock3->tock("loop");
	return maxpbuf;
}

//' Get the log-probability for observed secondary infections given the number of primary infections
//'
//' The probability depends only on the number of observed secondary cases, not which individuals have secondary infection.
//'
//' @param obs Number of observed infections among potential secondary cases
//' @param mis Number of missed observations among potential secondary cases
//' @param pri Number of primary infections
//' @param m Number of potential secondary cases
//' @param gdat pre-computed log-probabilities for secondary infection numbers
//' @l_Sn log of test sensitivity
//' @l_Sp log of test specificity
//' @l_1mSn log of one minus test sensitivity
//' @l_1mSp log of one minus test specificity
//' @return maximum value among the elements filled-in to pbuf
double lprob_sec_obs_pri(int obs, int mis, int pri, int m, Gdat gdat,
												 double l_Sn, double l_Sp, double l_1mSn, double l_1mSp){
	// r[0...m] are log-probabilities for j=0...m secondary infections
	double *r = g(gdat,pri,m);
  // Performance analysis using RcppClock
  //   ticker  time    prop
  //   all     2439.2  100.0%
  //   loop    2160.8   88.6%
  //   mse      125.3    5.1%
  // Peformance comments: cost of malloc seems to be 6.3%, which is not too much. mse is also cheap. Cost is in the loop.

  // clock2->tick("all");
	// How much memory is needed?
	int buf_size = 0;
	for(int j = 0; j <= m; j++){
		for(int e = std::max(0,j-(m-mis)); e <= std::min(mis,j); e++){
			buf_size += std::min(j-e,obs) - std::max(0,obs+j-e-(m-mis)) + 1;
		}
	}
	double *pbuf = (double *)malloc(buf_size*sizeof(double));
	int n = 0; // index marker for pbuf
	double maxpbuf = 0.0;
	// clock2->tick("loop");
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
			//Rprintf("k=%d,m=%d,j=%d,e=%d\n",pri,m,j,i);
			//Rprintf("  g[%d](%d,%d)=%.4f\n",pri,m,j,exp(r[j]));
			// probability of observed secondary cases
			double maxpbufi = lprob_obs_sec(pbuf+n, m-mis, j-e, obs, r[j] - lmj + log(coef), l_Sn, l_Sp, l_1mSn, l_1mSp);

			// double midprob = 0.0;
			// for(int ii = 0; ii < std::min(j-e,obs)-std::max(0,obs+j-e-(m-mis))+1; ii++){
			//   Rprintf("   [%.4f]\n",exp(pbuf[n+ii]-(r[j]-lmj+log(coef))));
			//   midprob += exp(pbuf[n+ii]-(r[j]-lmj+log(coef)));
			// }
			// Rprintf("  h(%d,%d,%d) = %.4f\n",m,obs,j,midprob);

			if(n==0 || maxpbufi > maxpbuf) maxpbuf = maxpbufi;
			n += std::min(j-e,obs)-std::max(0,obs+j-e-(m-mis))+1;
		}
	}
	// clock2->tock("loop");
	// clock2->tick("mse");
	if(maxpbuf==-1.0/0.0){
	  free(pbuf);
	  return maxpbuf;
	}
  double prob = 0.0;
  for(int i = 0; i < n; i++){
    prob += exp(pbuf[i] - maxpbuf);
  }
  // clock2->tock("mse");
  // clock2->tock("all");
  free(pbuf);
  return log(prob) + maxpbuf;
}

//' Get the log-probability for observed infections given primary infections
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
	// probability of primary infection observation
	// if SAR is 0, then only primary infections are infected, so get entire household observation likelihood
	double pprim;
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
	}else{
	  int a=0,b=0;
	  for(int i = 0; i < hh_size; i++){
	    if(obs[i]!=-1){ // if individual is not a missed observation
	      a +=  obs[i] && pri[i];
	      b += !obs[i] && pri[i];
	    }
	  }
	  pprim = mul(a,l_Sn) + mul(b,l_1mSn);
	  //Rprintf("P(Y_(o) cap Y_(p) | Y_(p)): %.4f\n",exp(pprim));
	  if(pprim == -1.0 / 0.0) return pprim; // shortcut
	}

	// Counts needed for secondary observation probability
	int k=0,mis=0,obs2=0,m=0;
	for(int i = 0; i < hh_size; i++){
		if(pri[i]){
			k++; // primary cases
		}else{
			m++; // potential secondary case
			mis  += (obs[i]==-1);
			obs2 += (obs[i]==1);
		}
	}
	double psec = lprob_sec_obs_pri(obs2, mis, k, m, gdat, l_Sn, l_Sp, l_1mSn, l_1mSp);
	return pprim + psec;
}

//' Probability of primary infection
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

//' Draw a sample for primary infection
//'
//' @param pri memory into which primary infection will be written as binary array
//' @param p probabilities of primary infection
//' @param n household size
void pri_sample(int *pri, double *p, int n){
	double *draws = (Rcpp::runif(n)).begin();
	for(int i = 0; i < n; i++){
		pri[i] = (draws[i] < p[i]);
	}
}

//' Get the log-probability for observations in a household
//'
//' @param obs Array of observations. -1=missed, 0=observed negative, 1=observed positive
//' @param prim Array of primary infection probabilities
//' @param hh_size Size of household
//' @param q Household secondary infection rate parameter
//' @param Sn test sensitivity
//' @param Sp test specificity
//' @param T Threshold for simulations. T<0: never simulate. T=0: always simulate. Otherwise, simulate for households over size T
//' @param M Number of simulations to compute over
//' @return log-probability
double lprob_obs(int* obs, double *prim, const int hh_size, double q,
								const double Sn, const double Sp, const int T, int M){
	// T < 0: never simulate
	// T = 0: always simulate
	// T > 0: simulate for households over size T
	// M: number of simulation samples

	// Performance analysis using RcppClock
	//   ticker    time    prop
	//   hh_prob  145.8  100.0%
	//   lgdat      7.5    5.1%
	//   sim_lse    0.4    0.3%
	//   simulate 137.3   94.2%
	// Peformance comments: Cost of lgdat per household is relatively low - worthwhile to keep. Cost is in simulate, which is best place

	// generate log-probabilities for secondary infection
	// clock1->tick("lgdat");
	Gdat gdat = lgdat(hh_size,q);
	// clock1->tock("lgdat");

	// allocate for primary infection array, and log of primary inf prob
	int *pri = (int *)malloc((hh_size)*sizeof(int));
	memset(pri,0,hh_size*sizeof(int));

	// Either simulate or iterate
	if(T < 0 || hh_size < T){
		// clock1->tick("iterate");
		// iterate over possible primary infections, using binary vector

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
		for(int64_t priv = 0; priv < ((int64_t)1<<hh_size); priv++){
			for(int i = 0; i < hh_size; i++){
				pri[i] = (priv >> i) & 1;
			}

			//Rprintf("Y_(p): ");
		  //print_vec(pri,hh_size);

			// probability of primary infection
			double pprim = lprob_pri(pri, lprim, l_1mprim, hh_size);

			//Rprintf("P(Y_(p)): %.4f\n", exp(pprim));

			// probability of observation given primary infections
			pbuf[priv] = pprim + lprob_obs_pri(obs, pri, hh_size, gdat, log(Sn), log(Sp), log(1-Sn), log(1-Sp));
			if(pri==0 || pbuf[priv] > max_prob) max_prob = pbuf[priv];
			//Rprintf("P(Y_(o) | Y_(p)) = %.4f\n", exp(pbuf[priv]));
			//Rprintf("\n");
		}
		// clock1->tock("iterate");
		free(pri);
		gfree(gdat);
		free(lprim);
		free(l_1mprim);

		// log-sum-exp
		// clock1->tick("it_lse");
		if(max_prob == -1.0/0.0){
		  free(pbuf);
		  return max_prob;
		}
		double prob = 0;
		for(int64_t priv = 0; priv < ((int64_t)1<<hh_size); priv++){
			prob += exp(pbuf[priv] - max_prob);
		}
		// clock1->tock("it_lse");
		free(pbuf);
		return log(prob) + max_prob;
	}else{
		// clock1->tick("simulate");
		// simulate infections
		double *pbuf = (double *)malloc(M*sizeof(double));
		double max_prob = 0;
		for(int i = 0; i < M; i++){
			// generate a sample primary infection into pri
			pri_sample(pri, prim, hh_size);
			// probability of observation given primary infection
			pbuf[i] = lprob_obs_pri(obs, pri, hh_size, gdat, log(Sn), log(Sp), log(1-Sn), log(1-Sp));
			if(pri==0 || pbuf[i] > max_prob) max_prob = pbuf[i];
		}
		// clock1->tock("simulate");
		free(pri);
		gfree(gdat);

		// log-sum-exp
		// clock1->tick("sim_lse");
		if(max_prob == -1.0/0.0){
		  free(pbuf);
		  return max_prob;
		}
		double prob = 0;
		for(int i = 0; i < M; i++){
			prob += exp(pbuf[i] - max_prob);
		}
		// clock1->tock("sim_lse");
		free(pbuf);
		return log(prob) + max_prob - log(M);
	}
}

//' Get the log-probability for observations across multiple households
//'
//' @param obs Array of size sum(hh_sizes) of observations. -1=missed, 0=observed negative, 1=observed positive
//' @param prim Array of size sum(hh_sizes) of primary infection probabilities
//' @param hh_size Array of size n of household sizes
//' @param q Array of size n of secondary infection rates per household
//' @param n Number of households
//' @param Sn test sensitivity
//' @param Sp test specificity
//' @param T Threshold for simulations. T<0: never simulate. T=0: always simulate. Otherwise, simulate for households over size T
//' @param M Number of simulations to compute over
//' @return log-probability
double lprob_obs_hhs(int *obs, double *prim, int *hh_sizes, double *q, const int n,
										 const double Sn, const double Sp, const int T, const int M){
	double lprob = 0.0;
	int i = 0;
	// clock1 = new Rcpp::Clock();
	// clock2 = new Rcpp::Clock();
	// clock3 = new Rcpp::Clock();

	for(int hh = 0; hh < n; hh++){
		// clock1->tick("hh_prob");
		lprob += lprob_obs(obs+i, prim+i, hh_sizes[hh], q[hh], Sn, Sp, T, M);
		// clock1->tock("hh_prob");
		i += hh_sizes[hh];
	}
	// clock1->stop("t1");
	// clock2->stop("t2");
	// clock3->stop("t3");
	// delete clock1;
	// delete clock2;
	// delete clock3;
	return lprob;
}
