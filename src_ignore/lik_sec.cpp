#include <Rcpp.h>
#include <math.h>
#include "lik_sec.h"
#include "pfi.h"

// [[Rcpp::export(.lik_sec)]]
Rcpp::NumericVector lik_sec_R(Rcpp::IntegerVector hhi,
                               Rcpp::NumericVector pi,
                               double sar,
                               int M){
  int n = hhi.size();
  Rcpp::NumericVector fi(n);

  prob_infected(fi.begin(),hhi.begin(),pi.begin(),n,sar,M);
  return(fi);
}

// p: primary cases
// n: household size
// x, maxx: pre-computed generating functions for 2ndary infections
double prob_any_sec(int p, int n, double *x, int max_hh){
	double prob = 0;
	Rprintf("g(%d,%d,%d)\n",p,n,max_hh);
	double *r = g(p,n-p,x,max_hh);
	Rprintf("r=");
	for(int i = 0; i <= n-p; i++) Rprintf("%.4f ",r[i]);Rprintf("\n");
	for(int m = 1; m <= n-p; m++){
		prob = prob + r[m]*m/(n-p);
	}
	return(prob);
}

// probabilities for actual infection based on primary infection probabilities
// fi: probability of final infection (to be filled)
// hh: household identifiers
// pi: probabilities for primary infection
// n: length of hh and pi arrays
// q: secondary attack rate
// M: iterations for simulation
void prob_infected(double *fi, int *hh, double *pi, int n, double q, int M){
	// get max hh size
	int max_hh = 0;
	int last_hh = 0;
	for(int i = 1; i < n; i++){
		if(hh[last_hh] != hh[i]){
			if(max_hh < i - last_hh) max_hh = i - last_hh;
			last_hh = i;
		}
	}
	if(max_hh < n - last_hh) max_hh = n - last_hh;
	Rprintf("max_hh=%d\n",max_hh);
	// initialize fi to 0
	memset(fi,0,n*sizeof(double));
	// get needed storage size for x
	int xsize = g_ind(max_hh,max_hh,0)+1;
	Rprintf("xsize=%d\n",xsize);
	double *x = (double *)malloc(xsize*sizeof(double));
	// pre-compute generating functions using max hh size
	g_dat(x, max_hh, q);
	// allocate storage for simulated primary cases
	int *sim_prim = (int *)malloc(n*sizeof(int));
	// vector for generated random numbers
	Rcpp::NumericVector rn_pi(n);
	for(int i = 0; i < n; i++) Rprintf("%.3f ",pi[i]);Rprintf("\n");
	// Run simulations
	for(int j = 0; j < M; j++){
		// simulate primary infections
		// generate n random numbers in (0,1)
		rn_pi = Rcpp::runif(n);
		// when rn_pi < pi, simulate individual as infected.
		for(int i = 0; i < n; i++){
			sim_prim[i] = rn_pi.begin()[i] < pi[i]; // ilogit
		}
		Rprintf("\nsim_prim=");
		for(int i = 0; i < n; i++) Rprintf("%d ",sim_prim[i]);Rprintf("\n");
		// point to first household
		int hh_i = 0;
		for(int i = 0; i <= n; i++){ // array ends at n-1, but we need to process last hh
			// check if new household to process
			if(i == n || hh[hh_i] != hh[i]){
				Rprintf("hh=%d\n",hh[hh_i]);
				// household size
				int ni = i - hh_i;
				Rprintf("hh_size=%d\n",ni);
				// find out how many primary infections
				int ki = 0;
				for(int ii = 0; ii < ni; ii++){
					if(sim_prim[hh_i + ii]) ki++;
				}
				Rprintf("ki=%d\n",ki);
				if(ki > 0){
					// get the probaiblity of infection for all secondary cases under simulation scenario
					double psec;
					if(ki == ni){
						psec = 0;
					}else{
						psec = prob_any_sec(ki, ni, x, max_hh);
					}
					Rprintf("psec=%.4f\n",psec);
					// update probability using either 1 (for simulated positives) or psec
					// (for potential secondary infections)
					for(int ii = 0; ii < ni; ii++){
						if(sim_prim[hh_i + ii]){
							fi[hh_i + ii] += 1;
						}else{
							fi[hh_i + ii] += psec;
						}
					}
				}
				// point to new household
				hh_i = i;
			}
		}
	}
	free(x);
	free(sim_prim);
	// divide probabilities by M
	for(int i = 0; i < n; i++){
		fi[i] = fi[i] / M;
	}
}
