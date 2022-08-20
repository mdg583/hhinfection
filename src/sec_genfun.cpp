#include <Rcpp.h>
#include <stdint.h>
#include "sec_genfun.h"

/**
 * Generating function for secondary infection numbers
 *
 * Strategy: the computation of secondary infection numbers is doubly recursive.
 * We ultimately need to compute probabilities for multiple primary and secondary
 * infection numbers, for a fixed household size and probability of secondary
 * infection. The best strategy seems to be to compute the result for all household
 * sizes, primary infection numbers and secondary infection numbers up to a maximum
 * household size.
 */

/**
 * Resulting probabilities will be placed in a block of memory, indexed using 3 variables:
 *  k: number of primary infections
 *  m: number of susceptibles to secondary infection
 *  j: number of secondary infections
 * 
 * Each k and m corresponds to a block of memory of size m+1 holding a degree m polynomial which
 * is the generating function for the given k and m. The coefficients correspond to probabilities
 * for j=0 ... j=m.
 * 
 * A special function g_ind(max_hh,k,m) will give the relative address of a desired generating
 * function in the pre-computed block of memory. Let n=max_hh. Note that:
 * - The memory needed for one g.f. for m=mi and any k is:
 *    mi+1
 * - The memory needed for all g.f.'s mi<m for any k is:
 *    sum_{mi=0}^{m-1} mi+1
 * - The memory needed for all g.f.'s ki<k is:
 *    sum_{ki=0}^{k-1} sum_{mi=0}^{n-ki} mi+1
 * - The memory needed for all g.f.'s ki=k,mi<m and for all ki<k is:
 *    ( sum_{mi=0}^{m-1} mi+1 ) + ( sum_{ki=0}^{k-1} sum_{mi=0}^{n-ki} mi+1 )
 * This is the index of the start of the g.f. for ki=k,mi=m.
 * 
 * Using Wolfram Alpha, this sum simplifies to:
 *   k*(k^2 - 3*k*n - 6*k + 3*n^2 + 12*n + 11)/6 + m*(m + 1)/2
 */ 

// Get the index for the generating function for k,m given maximum household size n
long g_ind(long n, long k, long m){
  if(k+m > n) Rcpp::stop("require k+m <= n");
  if(n > 1000) Rcpp::stop("n too large, require n <= 1000"); // overflow happens somehwere around 1285
  return (k*(k*k - 3*k*n - 6*k + 3*n*n + 12*n + 11))/6 + (m*(m + 1))/2;
}

// Workhorse function: precompute probabilities for all possible k,m up to maximum household size n
// x: pre-allocated space for probabilities
// n: max hh size
// q: parameter for secondary attack rate
void g_dat(double *x, int n, double q){
  double *r;      // will point to polynomial for m
  double *r_prev; // will point to polynomial for m-1
  double *g1;     // recursive reference to gen fun polynomial for 1,m-j
  int k,m;
  // Clear memory. Size is last index + 1
  int size = g_ind(n,n,0)+1;
  memset(x,0,size*sizeof(double));

  // build up from k=0 ... n primary infections
  for(k = 0; k <= n; k++){
    // build up from m=0 ... n-k secondary infections
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
            r[i] = r[i] + r_prev[i] * qki; // qki * x^j
            // To compute (1-qki)*g1*x^(j+1), iterate over terms in g1
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

// Create pre-computed data for secondary infection probabilities for all household sizes up to max_hh
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
