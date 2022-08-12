#include <Rcpp.h>
#include <math.h>
#include "lik_obs.h"
#include "pfi.h"
#include "pfi_rg.h"

void print_vec(int *v, int n){
	for(int i = 0; i < n; i++){
		Rprintf("%d ",v[i]);
	}
	Rprintf("\n");
}

void print_bivec(int v, int n){
	for(int i = 0; i < n; i++){
		Rprintf("%d ",(v >> i) & 1);
	}
	Rprintf("\n");
}

void print_flvec(double *v, int n){
	for(int i = 0; i < n; i++){
		Rprintf("%.2f ",v[i]);
	}
	Rprintf("\n");
}

// n choose r hopefully without overflow
long choose(long n, long r){
  if(2*r > n) r = n-r; // swap to easier case
  // n-r choose 0 = n-r
  long res = n-r;
  // (n-r+1 choose 0+1) = (n-r choose 0) * ((n-r+1)/(0+1))
  for(int i = 1; i <= r; i++){
    res = (res * (n-r+i))/i;
  }
  return res;
}

// compute a trinomial coefficient, (b+c+d) choose b,c,d
long trinom(long b, long c, long d){
  return choose(b+c,b)*choose(b+c+d,b+c);
}

void bcoefs(long *coefs, long n){
	coefs[0] = 1;
	coefs[1] = n;
	for(int i = 2; i <= n/2; i++){
		coefs[i] = coefs[i-1] * (n+1-i) / i;
	}
	for(int i = (n+1)/2; i <= n; i++){
		coefs[i] = coefs[n-i];
	}
}

// count bits in long integer
// Note that most compilers (GCC, MS, vlang) on 64-bit platforms have a function to do this using a single assembly instruction.
// However R could support some other platform like ARM.
// Adapted from https://stackoverflow.com/questions/51387998/count-bits-1-on-an-integer-as-fast-as-gcc-builtin-popcountint#answer-51388543
int count_bits(int64_t x){
	int count = 0;
	unsigned char *ptr = (unsigned char *)&x;
		for (int i=0;i<(int)sizeof(int64_t);i++) {
				count += bitcount[ptr[i]];
		}
		return count;
}

// Multiply an int by a double, with the product set to 0 for n=0 for all x
// This is to ensure that 0*Inf = 0. The goal here is to minize operations, especially avoid jump
// This helps somewhat with performance but may not be much. But it is critical code, running
// millions or billions of times.
inline double mul1(int n, double x){
	double nx = n*x;
	if(n==0) nx = 0;
	return nx;
}
inline double mul2(int n, double x){
	// bitwise all 1's when n != 0
	int64_t n_zero = 0 - (int64_t)(n!=0);
	double nx = n*x;
	int64_t nx2 = (*(int64_t*)&nx) & n_zero;
	return (*(double*)&nx2);
}
inline double mul3(int n, double x){
	int64_t n_zero,nx2;
	double nx,nx3;
	// bitwise all 1's when n != 0
	n_zero = 0 - (int64_t)(n!=0);
	nx = n*x;
	memcpy(&nx2, &nx, sizeof nx); // (*(long long*)&nx)
	nx2 = nx2 & n_zero; // set bits to 0 if n==0
	memcpy(&nx3, &nx2, sizeof nx2); // (*(long long*)&nx)
	return nx3; // interpret bits as double
}
#define MUL mul3

// Log probability of observed outcome given actual outcome, using sensitivity and specificity
// obs and act are 64-bit integers, each bit corresponds to an individual. 1=infected, 0=not infected
// mis: unobserved individuals
inline double lprob_obs_act(int64_t obs, int64_t act, int64_t mis, int n,
                            double l_Sn, double l_Sp, double l_1mSn, double l_1mSp){
	// Prob obs|act = Sn, Prob obs|~act	= 1-Sp, Prob ~obs|act	= 1-Sn, Prob ~obs|~act = Sp
	double prob = 0;
	int mask = ((1<<n)-1) & ~mis;
	int x1 = count_bits((obs & act) & mask);
	int x2 = count_bits((~obs & act) & mask);
	int x3 = count_bits((obs & ~act) & mask);
	int x4 = count_bits((~obs & ~act) & mask);
	prob = MUL(x1,l_Sn) + MUL(x2,l_1mSn) + MUL(x3,l_1mSp) + MUL(x4,l_Sp);
	return prob;
}

double lprob_sec_obs_pri(int obs, int pri, int mis, int m, double *lg_data, int max_hh,
                         double l_Sn, double l_Sp, double l_1mSn, double l_1mSp){
  double *r = g(pri,m,lg_data,max_hh);
  // an upper bound on computations is m*mis*obs <= min(mis*m^2,obs*m^2)
  // count how much memory is needed
  int xsize=0;
  for(int j = 0; j <= m; j++){
    for(int i = 0; i <= mis; i++){
      if(j+i < obs){
        xsize += ((j+i)*(j+i+1))/2;
      }else{
        xsize += (obs*(obs+1))/2;
      }
    }
  }
  double *lprobs = (double *)malloc(xsize * sizeof(double));
  double maxlprobs = 0;
  int ii = 0;
  // secondary infections
  for(int j = 0; j <= m; j++){
    // true positives among missed observations
    for(int i = 0; i <= mis; i++){
      long lmisi = log(choose(mis,i));// update rolling coefficient?
      // iterator for groups
      long mc = trinom(j+i,obs,m-(j+i)-obs);
      int alim = obs;
      if(j+i < obs) alim = j+i;
      for(int a = 0; a <= alim; a++){
        int b = j+i-a;
        int c = obs-a;
        int d = m-a-b-c;
        lprobs[ii] = r[j] + lmisi + log(mc) + l_Sn*a + l_1mSp*b + l_1mSn*c + l_Sp*d;
        if(ii==0 || lprobs[ii] > maxlprobs) maxlprobs = lprobs[ii];
        ii++;
        // next multinomial coefficient
        mc = (mc * (b-1)*(c-1))/((a+1)*(d+1));
      }
    }
  }
  // max-sum-exp
  double lprob = 0.0;
  for(int i = 0; i < ii; i++){
    lprob += exp(lprobs[i] - maxlprobs);
  }
  free(lprobs);
  return log(lprob) + maxlprobs;
}

// Log probabiliy of observation given primary infection and parameters
double lprob_obs_pri_sum(int64_t obs, int64_t pri, int64_t mis, int hh_size, double *lg_data, int max_hh,
                        double l_Sn, double l_Sp, double l_1mSn, double l_1mSp){



  int k = count_bits(pri); // num primary infection
  int m = hh_size - k;		 // sum susceptible individuals
  // separate obs according to primary infection
  int obs_p = 0;
  int mis_p = 0;
  int obs_np = 0;
  int mis_np = 0;
  for(int i = 0, ip = 0, inp = 0; i < hh_size; i++){
    // if primary, add bit to obs_p
    int prii = ((pri >> i) & 1);
    int obsi = ((obs >> i) & 1);
    int misi = ((mis >> i) & 1);
    obs_p += (1 << ip) * obsi * prii;
    mis_p += (1 << ip) * misi * prii;
    ip += prii;
    // if not primary, add bit to obs_np
    obs_np += (1 << inp) * obsi * (1-prii);
    mis_np += (1 << inp) * misi * (1-prii);
    inp += (1-prii);
  }
  // create a vector m - n 1's for primary infection
  pri = (1 << (hh_size-m)) - 1;

  // probability of observing primary infections
  double pprim = lprob_obs_act(obs_p, pri, mis_p, k, l_Sn, l_Sp, l_1mSn, l_1mSp);
  if(pprim == -1.0 / 0.0) return pprim; // shortcut

  // log prob of secondary infection numbers
  double *r = g(k,m,lg_data,max_hh);

  // calculate log binomial coefficients
  long		binco[64];
  double lbinco[64];
  bcoefs(binco,m);
  for(int i = 0; i <= m; i++) { lbinco[i] = log(binco[i]); }

  // get max possible value for Pr(vo,v|vp) : g(n-m, m, lg_data, max_hh) + lprob_obs_act(obs_np, sec)
  // This is to reduce floating point error for log-sum-exp calculation
  double max_lprob = r[0];
  for(int i = 0; i <= m; i++){
    if(r[i] > max_lprob) max_lprob = r[i] - lbinco[i];
  }
  max_lprob = max_lprob + pprim + (m-count_bits(mis_np))*fmax(fmax(l_Sn,l_Sp),fmax(l_1mSn,l_1mSp));

  // iterate over all possible secondary infections
  // for m=20, this is 1 million iterations. For m=30, its 1 billion. Hopefully we don't get that many.
  double prob = 0;
  if(pri == 0){
    // no primary infections, so r[j]=-Inf for j>0 (ie sec>0)
    int j = count_bits(0);
    prob += exp(r[j] - lbinco[j] + pprim +
      lprob_obs_act(obs_np, 0, mis_np, m, l_Sn, l_Sp, l_1mSn, l_1mSp) - max_lprob);
  }else{
    for(int64_t sec = 0; sec < (1 << m); sec++){
      // probability of secondary infection and probability of observation
      int j = count_bits(sec);
      prob += exp(r[j] - lbinco[j] + pprim +
        lprob_obs_act(obs_np, sec, mis_np, m, l_Sn, l_Sp, l_1mSn, l_1mSp) - max_lprob);
    }
  }
  return log(prob) + max_lprob;
}

// Log probabiliy of observation given primary infection and parameters
double lprob_obs_pri_it(int64_t obs, int64_t pri, int64_t mis, int hh_size, double *lg_data, int max_hh,
                        double l_Sn, double l_Sp, double l_1mSn, double l_1mSp){
	int k = count_bits(pri); // num primary infection
	int m = hh_size - k;		 // sum susceptible individuals
	// separate obs according to primary infection
	int obs_p = 0;
	int mis_p = 0;
	int obs_np = 0;
	int mis_np = 0;
	for(int i = 0, ip = 0, inp = 0; i < hh_size; i++){
		// if primary, add bit to obs_p
		int prii = ((pri >> i) & 1);
		int obsi = ((obs >> i) & 1);
		int misi = ((mis >> i) & 1);
		obs_p += (1 << ip) * obsi * prii;
		mis_p += (1 << ip) * misi * prii;
		ip += prii;
		// if not primary, add bit to obs_np
		obs_np += (1 << inp) * obsi * (1-prii);
		mis_np += (1 << inp) * misi * (1-prii);
		inp += (1-prii);
	}
	// create a vector m - n 1's for primary infection
	pri = (1 << (hh_size-m)) - 1;

	// probability of observing primary infections
	double pprim = lprob_obs_act(obs_p, pri, mis_p, k, l_Sn, l_Sp, l_1mSn, l_1mSp);
	if(pprim == -1.0 / 0.0) return pprim; // shortcut

	// log prob of secondary infection numbers
	double *r = g(k,m,lg_data,max_hh);

	// calculate log binomial coefficients
	long		binco[64];
	double lbinco[64];
	bcoefs(binco,m);
	for(int i = 0; i <= m; i++) { lbinco[i] = log(binco[i]); }

	// get max possible value for Pr(vo,v|vp) : g(n-m, m, lg_data, max_hh) + lprob_obs_act(obs_np, sec)
	// This is to reduce floating point error for log-sum-exp calculation
	double max_lprob = r[0];
	for(int i = 0; i <= m; i++){
		if(r[i] > max_lprob) max_lprob = r[i] - lbinco[i];
	}
	max_lprob = max_lprob + pprim + (m-count_bits(mis_np))*fmax(fmax(l_Sn,l_Sp),fmax(l_1mSn,l_1mSp));

	// iterate over all possible secondary infections
	// for m=20, this is 1 million iterations. For m=30, its 1 billion. Hopefully we don't get that many.
	double prob = 0;
	if(pri == 0){
		// no primary infections, so r[j]=-Inf for j>0 (ie sec>0)
		int j = count_bits(0);
		prob += exp(r[j] - lbinco[j] + pprim +
		  lprob_obs_act(obs_np, 0, mis_np, m, l_Sn, l_Sp, l_1mSn, l_1mSp) - max_lprob);
	}else{
		for(int64_t sec = 0; sec < (1 << m); sec++){
			// probability of secondary infection and probability of observation
			int j = count_bits(sec);
			prob += exp(r[j] - lbinco[j] + pprim +
			  lprob_obs_act(obs_np, sec, mis_np, m, l_Sn, l_Sp, l_1mSn, l_1mSp) - max_lprob);
		}
	}
	return log(prob) + max_lprob;
}

double lprob_obs_pri_sim(int64_t obs, int64_t pri, int64_t mis, const int hh_size, const double q,
                         const double l_Sn, const double l_Sp, const double l_1mSn, const double l_1mSp,
                         double *pbuf2, const int M2){
  int priv[hh_size];
  int secv[hh_size];
  int M[hh_size*hh_size]; // infection matrix

  // convert binary vectors into arrays
  for(int i = 0; i < hh_size; i++){
    priv[i] = (pri >> i) & 1;
  }

  double max_lprob = 0;
  for(int k = 0; k < M2; k++){
    memset(secv,0,hh_size*sizeof(int));
    // generate a random graph
    randomgraph(M,hh_size,q);

    // infect secondary cases based on graph
    for(int i=0; i<hh_size; i++){
      if(priv[i] == 1) cluster(M,hh_size,i,secv);
    }

    // convert secondary cases back to binary vector
    int64_t sec = 0;
    for(int i = 0; i < hh_size; i++){
      sec += (1<<i) * secv[i];
    }

    // get probability of observation
    pbuf2[k] = lprob_obs_act(obs, sec, mis, hh_size, l_Sn, l_Sp, l_1mSn, l_1mSp);
    if(k==0 || pbuf2[k] > max_lprob) max_lprob = pbuf2[k];
  }
  // log-sum-exp
  double prob = 0;
  for(int k = 0; k < M2; k++){
    prob += exp(pbuf2[k] - max_lprob);
  }

  // get probability of observation based on simulated secondary infection
  return log(prob) - log(M2) + max_lprob;
}

// get probability of primary infeciton
double lprob_pri(int64_t pri, double *lprim, double *l_1mprim, int n){
	double lprob = 0;
	for(int i = 0; i < n; i++){
		int prii = (pri >> i) & 1;
		lprob += prii*lprim[i] + (1-prii)*l_1mprim[i];
	}
	return lprob;
}

// Generate a binary integer representing with first n bits being 1 with ind. probs in p
int64_t rbool(double *p, int n){
	double *draws = (Rcpp::runif(n)).begin();
	int64_t x = 0;
	for(int i = 0; i < n; i++){
		x += (draws[i] < p[i])*(1 << i);
	}
	return x;
}

double lprob_obs_it(int64_t obsv, int64_t misv, double *lprim, double *l_1mprim, const int hh_size,
                    double *lg_data, const int max_hh, double q,
                    const double l_Sn, const double l_Sp, const double l_1mSn, const double l_1mSp,
                    double *pbuf1, double *pbuf2, const int M2, int T2){
  // iterate over all possible primary infections
  memset(pbuf1,0,(1 << hh_size)*sizeof(double));
  double max_prob = 0;
  for(int64_t pri = 0; pri < (1<<hh_size); pri++){
    double prob = 0;
    // get probability of primary infection
    prob += lprob_pri(pri, lprim, l_1mprim, hh_size);
    // get probability of observation
    if(T2 != 0 && (1<<hh_size) > T2){
      prob += lprob_obs_pri_sim(obsv, pri, misv, hh_size, q, l_Sn, l_Sp, l_1mSn, l_1mSp, pbuf2, M2);
    }else{
      prob += lprob_obs_pri_it(obsv, pri, misv, hh_size, lg_data, max_hh, l_Sn, l_Sp, l_1mSn, l_1mSp);
    }
    pbuf1[pri] = prob;
    if(pri==0 || prob > max_prob) max_prob = prob;
  }
  // log-sum-exp
  double prob = 0;
  for(int64_t i = 0; i < (1<<hh_size); i++){
    prob += exp(pbuf1[i] - max_prob);
  }
  return log(prob) + max_prob;
}

double lprob_obs_sim(int64_t obsv, int64_t misv, double *prim, const int hh_size,
                    double *lg_data, const int max_hh, double q,
                    const double l_Sn, const double l_Sp, const double l_1mSn, const double l_1mSp,
                    double *pbuf1, double *pbuf2, const int M1, const int M2, int T2){
  // Generate M simulations of primary infection
  memset(pbuf1,0,M1*sizeof(double));
  double max_prob = 0;
  for(int i = 0; i < M1; i++){
    int64_t pri = rbool(prim,hh_size);
    double prob = 0;
    if(T2 != 0 && hh_size > T2){
      prob += lprob_obs_pri_sim(obsv, pri, misv, hh_size, q, l_Sn, l_Sp, l_1mSn, l_1mSp, pbuf2, M2);
    }else{
      prob += lprob_obs_pri_it(obsv, pri, misv, hh_size, lg_data, max_hh, l_Sn, l_Sp, l_1mSn, l_1mSp);
    }
    pbuf1[i] = prob;
    if(pri==0 || prob > max_prob){
      max_prob = prob;
    }
  }
  // log-sum-exp
  double prob = 0;
  for(int i = 0; i < M1; i++){
    prob += exp(pbuf1[i] - max_prob);
  }
  return log(prob) + max_prob - log(M1);
}

// Log probabiliy of household observation, iterative (brute force)
double lprob_obs(int *obs, double *prim, const int hh_size, double *lg_data, const int max_hh, const double q,
                 const double l_Sn, const double l_Sp, const double l_1mSn, const double l_1mSp,
                 double *pbuf1, double *pbuf2, const int M1, const int M2, const int T1, const int T2){
	// Convert boolean vectors into binary integers
	// interpret obs<0 as missing observations
	int64_t obsv = 0;
	int64_t misv = 0;
	for(int i = 0; i < hh_size; i++){
		obsv += ((int64_t)1<<i)*(obs[i] > 0);
		misv += ((int64_t)1<<i)*(obs[i] < 0);
	}

	// T1 is Household size threshold for simulating primary infection
	if(T1 != 0 && (1<<hh_size) > T1){
	  return lprob_obs_sim(obsv, misv, prim, hh_size, lg_data, max_hh, q,
                        l_Sn, l_Sp, l_1mSn, l_1mSp, pbuf1, pbuf2, M1, M2, T2);
	}else{
	  double lprim[64];
	  double l_1mprim[64];
	  for(int i = 0; i < hh_size; i++){
	    lprim[i] = log(prim[i]);
	    l_1mprim[i] = log(1-prim[i]);
	  }
	  return lprob_obs_it(obsv, misv, lprim, l_1mprim, hh_size, lg_data, max_hh, q,
                       l_Sn, l_Sp, l_1mSn, l_1mSp, pbuf1, pbuf2, M2, T2);
	}
}

// Log probabiliy of household observation, iterative (brute force)
double lprob_obs_hhs(int *obs, double *prim, int *hh, const int n, const double q,
                     const double Sn, const double Sp,
                     const int M1, const int M2, const int T1, const int T2){
  // get max hh size, as well as space needed for pbuf1, pbuf2
  int pbuf1_size = 0;
  int pbuf2_size = 0;
	int max_hh = 0;
	int last_hh = 0;
	for(int i = 1; i <= n; i++){
		if(i == n || hh[last_hh] != hh[i]){
		  int hh_size = i - last_hh;
			if(max_hh < hh_size) max_hh = hh_size;
			if(T1 != 0 && (1<<hh_size) > T1)  pbuf1_size = (std::max)(pbuf1_size, M1);
			if(T1 == 0 || (1<<hh_size) <= T1) pbuf1_size = (std::max)(pbuf1_size, (1<<hh_size));
			if(T2 != 0 && (1<<hh_size) > T2)  pbuf2_size = (std::max)(pbuf2_size, M2);
			if(T2 == 0 || (1<<hh_size) <= T2) pbuf2_size = (std::max)(pbuf2_size, (1<<hh_size));
			last_hh = i;
		}
	}
	if(max_hh > 62) Rcpp::stop("Household too large! Max 62 (though 30 more reasonable)");

	// Allocate enough space for pbuf, up to 2^max_hh or M. M==0 means always iterate
	double *pbuf1 = (double *)malloc(pbuf1_size*sizeof(double));
	double *pbuf2 = (double *)malloc(pbuf2_size*sizeof(double));
	// pre-compute generating functions using max hh size, convert to log probabilities
	int xsize = g_ind(max_hh,max_hh,0)+1;
	double *lg_data = (double *)malloc(xsize*sizeof(double));
	g_dat(lg_data, max_hh, q);
	for(int i = 0; i < xsize; i++) lg_data[i] = log(lg_data[i]);
  // sum log probability for each household
	double lprob = 0;
	int hh_i = 0;
	for(int i = 0; i <= n; i++){ // array ends at n-1, but we need to process last hh
		// check if new household to process
		if(i == n || hh[hh_i] != hh[i]){
			int ni = i - hh_i; // household index from hh_i to hh_i+ni-1 (size ni)
			// Compute probability for household
			lprob += lprob_obs(obs + hh_i, prim + hh_i, ni, lg_data, max_hh, q,
                      log(Sn), log(Sp), log(1-Sn), log(1-Sp),
                      pbuf1, pbuf2, M1, M2, T1, T2);
			// set start index for next household
			hh_i = i;
		}
	}
	free(pbuf1);
	free(pbuf2);
	free(lg_data);
	return lprob;
}
