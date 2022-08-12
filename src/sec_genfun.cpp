#include <Rcpp.h>
#include <stdint.h>
#include "sec_genfun.h"

/**
 * Generating function for secondary infection numbers
 *
 * Strategy: the computation of secondary infection numbers is multiply recursive.
 * We ultimately need to compute probabilities for multiple primary and secondary
 * infection numbers, for a fixed household size and probability of secondary
 * infection. The best strategy seems to be to compute the result for all household
 * sizes (k+m), primary infection numbers (k) and secondary infection numbers (j)
 * up to the desired household size.
 *
 * Data are in 3-dimensions:
 * 1. m (number of susceptibles)
 * 2. k (primary infections)
 * 3. j (number of secondary infections)
 *
 * Let N be the maximum household size
 * m = 0 ...
 * j = 0 ... m
 * k =
 */

/**
 * n: max hh size
 * k: primary infections
 * m: secondary infections
 * Get the memory-pointer offset of the first element of the generating function
 *
 * Memory used for a generating function with mi secondary infections:
 *   mi+1
 * Memory used for generating functions corresponding to ki primary infections:
 *   sum_(mi=0)^(n - ki)(mi + 1)
 * Memory used for all generating functions from ki=0 ... k-1 primary infections:
 *   sum_(ki=0)^(k - 1)[ sum_(mi=0)^(n - ki)(mi + 1) ]
 *
 * By wolfram alpha:
 * sum_(ki=0)^(k - 1)( sum_(mi=0)^(n - ki)(mi + 1)) + sum_(mi=0)^(m - 1)(mi + 1)
 *   = 1/6 k (k^2 - 3 k n - 6 k + 3 n^2 + 12 n + 11) + 1/2 m (m + 1)
 *   = k*(k^2 - 3*k*n - 6*k + 3*n^2 + 12*n + 11)/6 + m*(m + 1)/2
 */
long g_ind(long n, long k, long m){
  if(k+m > n) Rcpp::stop("require k+m <= n");
  if(n > 1200) Rcpp::stop("n too large, require n <= 1200"); // overflow happens somehwere around 1285
  return (k*(k*k - 3*k*n - 6*k + 3*n*n + 12*n + 11))/6 + (m*(m + 1))/2;
}

// Fill in data of pre-computed generating function probabilities
// n: max hh size
// q: parameter for secondary attack rate
void g_dat(double *x, int n, double q){
  double *r;      // will point to polynomial for m
  double *r_prev; // will point to polynomial for m-1
  double *g1;     // recursive reference to gen fun polynomial for 1,m-j
  int k,m;
  // Allocate memory for data. Size is last index + 1
  int size = g_ind(n,n,0)+1;
  memset(x,0,size*sizeof(double));
  // build up from k=0 ... n primary infections
  for(k = 0; k <= n; k++){
    // build up from m=0 ... n-k secondary infections (max hh size m+k=n)
    for(m = 0; m <= (n-k); m++){
      // get pointer to polynomial we are building, of degree m
      r = x + g_ind(n,k,m);
      if(m==0){
        r[0] = 1;
      }else{
        // In this case r_prev is the previous polynomial, of degree m-1
        // Iterate over r_prev
        for(int i = 0; i < m; i++){
          if(k+i==0){
            r[i] = r[i] + r_prev[i];
          }else{
            // Get a pointer to the polynomial for k=1,m=m-i-1, of degree m-i-1
            g1 = x + g_ind(n,1,m-i-1);
            double qki = pow(1-q,k+i);
            // qki * x^j
            r[i] = r[i] + r_prev[i] * qki;
            // (1-qki)*g1*x^[j+1]: iterate over terms in g1
            for(int j = 0; j <= m-i-1; j++){
              r[i+j+1] = r[i+j+1] + r_prev[i]*(1-qki)*g1[j];
            }
          }
        }
      }
      r_prev = r;
    }
  }
}

// Create pre-computed data for secondary infection probabilities for all household sizes up to
// max_hh
Gdat gdat(int max_hh, double q){
	Gdat ret;
	int size = g_ind(max_hh,max_hh,0)+1;
	ret.q = q;
	ret.max_hh = max_hh;
	ret.data = (double *)malloc(size*sizeof(double));
	g_dat(ret.data, max_hh, q);
	return ret;
}

// Create pre-computed data with a log transform
Gdat lgdat(int max_hh, double q){
	Gdat ret = gdat(max_hh, q);
	for(int i = 0; i < gsize(ret); i++){
		ret.data[i] = log(ret.data[i]);
	}
	return ret;
}

// Get the size needed to store gdat data
long gsize(Gdat dat){
  return g_ind(dat.max_hh,dat.max_hh,0)+1;
}

// Get a pointer to the generating function of secondary infection probaiblities
// corresponding to k,m from precomputed data gdat
double *g(Gdat gdat, int k, int m){
  return gdat.data + g_ind(gdat.max_hh, k, m);
}

// Free pre-generated data for secondary infections
void gfree(Gdat gdat){
	free(gdat.data);
}
