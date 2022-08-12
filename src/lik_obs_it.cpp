#include <Rcpp.h>
#include <math.h>
#include "sec_genfun.h"
#include "utilities.h"
#include "lik_obs_it.h"
#include "lik_obs.h"

//' Get the log-probability of observed infections given actual infections
//'
//' @param obs 64-bit vector of observed infection
//' @param act 64-bit vector of actual infection
//' @param mis 64-bit vector of individuals with no observation
//' @param n household size
//' @l_Sn log of test sensitivity
//' @l_Sp log of test specificity
//' @l_1mSn log of one minus test sensitivity
//' @l_1mSp log of one minus test specificity
//' @return log-probability of observation
inline double lprob_obs_act(int64_t obs, int64_t act, int64_t mis, int n,
														double l_Sn, double l_Sp, double l_1mSn, double l_1mSp){
	// Prob obs|act = Sn, Prob obs|~act	= 1-Sp, Prob ~obs|act	= 1-Sn, Prob ~obs|~act = Sp
	int64_t mask = ((1ULL << n) - 1) & ~mis;
	int x1 = count_bits((obs & act) & mask);
	int x2 = count_bits((~obs & act) & mask);
	int x3 = count_bits((obs & ~act) & mask);
	int x4 = count_bits((~obs & ~act) & mask);
	double prob = mul(x1,l_Sn) + mul(x2,l_1mSn) + mul(x3,l_1mSp) + mul(x4,l_Sp);
	return prob;
}

//' Get the log-probability for observed secondary infections given the number of primary infections, using iterative method
//'
//' This is for performance comparison only
//'
//' @param obs Binary vector of observed infections
//' @param mis Binary vector of missed observations
//' @param k Number of primary infections
//' @param m Number of potential secondary cases
//' @param gdat pre-computed log-probabilities for secondary infection numbers
//' @l_Sn log of test sensitivity
//' @l_Sp log of test specificity
//' @l_1mSn log of one minus test sensitivity
//' @l_1mSp log of one minus test specificity
//' @return log-probability
double lprob_sec_obs_pri_it(int64_t obs, int64_t mis, int k, int m, Gdat gdat,
												double l_Sn, double l_Sp, double l_1mSn, double l_1mSp){
	// r[0...m] are log-probabilities for j=0...m secondary infections
	double *r = g(gdat,k,m);
	// calculate log binomial coefficients
	long lbinco[64];
	bcoefs(lbinco,m);
	for(int j = 0; j <= m; j++) lbinco[j] = log(lbinco[j]);

	// get a reasonable upper bound for log prob, to reduce overflow error
	double max_lprob = r[0] - lbinco[0];
	for(int j = 1; j <= m; j++){
		if(r[j]-lbinco[j] > max_lprob) max_lprob = r[j] - lbinco[j];
	}
	max_lprob = max_lprob + mul(m-count_bits(mis),fmax(fmax(l_Sn,l_Sp),fmax(l_1mSn,l_1mSp)));
	if(max_lprob == -1.0/0.0) return max_lprob;

	// iterate over all possible secondary infections
	double prob = 0.0;
	for(int64_t sec = 0; sec < (1 << (int64_t)m); sec++){
		int j = count_bits(sec);
		prob += exp(r[j] - lbinco[j] + lprob_obs_act(obs, sec, mis, m, l_Sn, l_Sp, l_1mSn, l_1mSp) - max_lprob);
	}
	return log(prob) + max_lprob;
}

//' Get the log-probability for observed infections given primary infections (Iterative)
//'
//' This is for performance comparison only.
//'
//' @param obs Array of observations. -1=missed, 0=observed negative, 1=observed positive
//' @param pri Array of primary infections (0 or 1)
//' @param hh_size Size of household
//' @param gdat pre-computed log-probabilities for secondary infection numbers
//' @l_Sn log of test sensitivity
//' @l_Sp log of test specificity
//' @l_1mSn log of one minus test sensitivity
//' @l_1mSp log of one minus test specificity
//' @return log-probability
double lprob_obs_pri_it(int* obs, int* pri, const int hh_size, Gdat gdat,
												 const double l_Sn, const double l_Sp, const double l_1mSn, const double l_1mSp){
  // probability of primary infection observation
	int a=0,b=0,c=0,d=0;
	for(int i = 0; i < hh_size; i++){
		if(obs[i]!=-1){ // if individual is not a missed observation
			a += obs[i] && pri[i];
			b += obs[i] && pri[i];
			c += obs[i] && pri[i];
			d += obs[i] && pri[i];
		}
	}
	double pprim = mul(a,l_Sn) + mul(b,l_1mSn) + mul(c,l_1mSp) + mul(d,l_Sp);
	if(pprim == -1.0 / 0.0) return pprim; // shortcut

	// Create vectors for secondary observation probability
	int64_t obsv=0, misv=0; int k=0, m=0;
	for(int i = 0; i < hh_size; i++){
		if(pri[i]){
			k++; // primary cases
		}else{
			m++; // potential secondary case
			misv += (1<<i) * (obs[i]==-1);
			obsv += (1<<i) * (obs[i]==1);
		}
	}
	// probability of secondary observations
	double psec = lprob_sec_obs_pri_it(obsv, misv, k, m, gdat, l_Sn, l_Sp, l_1mSn, l_1mSp);
	return pprim + psec;
}

//' Get the log-probability for observations in a household, iterating over secondary cases
//'
//' For performance testing only
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
double lprob_obs_it(int* obs, double *prim, const int hh_size, double q,
								const double Sn, const double Sp, const int T, int M){
	// T < 0: never simulate
	// T = 0: always simulate
	// T > 0: simulate for households over size T
	// M: number of simulation samples

	// generate log-probabilities for secondary infection
	Gdat gdat = lgdat(hh_size,q);

	// allocate for primary infection array, and log of primary inf prob
	int *pri = (int *)malloc((hh_size)*sizeof(int));
	memset(pri,0,hh_size*sizeof(int));

	// Either simulate or iterate
	if(T < 0 || hh_size < T){
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
			// probability of primary infection
			double pprim = lprob_pri(pri, lprim, l_1mprim, hh_size);
			// probability of observation given primary infections
			pbuf[priv] = pprim + lprob_obs_pri_it(obs, pri, hh_size, gdat, log(Sn), log(Sp), log(1-Sn), log(1-Sp));
			if(pri==0 || pbuf[priv] > max_prob) max_prob = pbuf[priv];
		}
		free(pri);
		gfree(gdat);
		free(lprim);
		free(l_1mprim);

		// log-sum-exp
		if(max_prob == -1.0/0.0) return max_prob;
		double prob = 0;
		for(int64_t priv = 0; priv < ((int64_t)1<<hh_size); priv++){
			prob += exp(pbuf[priv] - max_prob);
		}
		free(pbuf);
		return log(prob) + max_prob;
	}else{
		// simulate infections
		double *pbuf = (double *)malloc(M*sizeof(double));
		double max_prob = 0;
		for(int i = 0; i < M; i++){
			// generate a sample primary infection into pri
			pri_sample(pri, prim, hh_size);
			// probability of observation given primary infection
			pbuf[i] = lprob_obs_pri_it(obs, pri, hh_size, gdat, log(Sn), log(Sp), log(1-Sn), log(1-Sp));
			if(pri==0 || pbuf[i] > max_prob) max_prob = pbuf[i];
		}
		free(pri);
		gfree(gdat);

		// log-sum-exp
		if(max_prob == -1.0/0.0) return max_prob;
		double prob = 0;
		for(int i = 0; i < M; i++){
			prob += exp(pbuf[i] - max_prob);
		}
		free(pbuf);
		return log(prob) + max_prob - log(M);
	}
}

//' Get the log-probability for observations across multiple households, iterating over secondary cases
//'
//' For performance test only
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
double lprob_obs_hhs_it(int *obs, double *prim, int *hh_sizes, double *q, const int n,
										 const double Sn, const double Sp, const int T, const int M){
	double lprob = 0.0;
	int i = 0;
	for(int hh = 0; hh < n; hh++){
		i += hh_sizes[hh];
		lprob += lprob_obs_it(obs+i, prim+i, hh_sizes[hh], q[hh], Sn, Sp, T, M);
	}
	return lprob;
}
