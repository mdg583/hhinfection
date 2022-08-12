#include <Rcpp.h>
#include <math.h>
#include "pfi.h"

// pfi.cpp is a simple implementation of a few household infection models for Glasgow ODL project 2022
// by Matthew Gibson

// [[Rcpp::export(.pfi_lik)]]
Rcpp::NumericVector pfi_lik_R(Rcpp::IntegerVector pi,
                               Rcpp::IntegerVector fi,
                               Rcpp::IntegerVector ni,
                               double sar){
  int m = pi.size();
  Rcpp::NumericVector lik(m);
  // pre-compute all generating functions for secondary infections
  int maxhh = -1;
  for(int i = 0; i < m; i++){
    if(ni[i] > maxhh) maxhh = ni[i];
  }
  Rprintf("Max hh size: %d\n", maxhh);
  int x_size = g_ind(maxhh,maxhh,0)+1;
  //Rprintf("g_dat alloc size: %d\n", x_size);
  double *x = (double *)malloc(x_size*sizeof(double));
  g_dat(x,maxhh,sar);
  // Each line of the vectors corresponds to a single household secondary infection probability calculation.
  for(int i = 0; i < m; i++){
    //Rprintf("g(%d,%d)[%d]",pi[i],ni[i]-pi[i],fi[i]-pi[i]);
    // get the generating function for pi positives and ni-pi susceptibles
    double *r = g(pi[i],ni[i]-pi[i],x,maxhh);
    // lookup likelihood for fi-pi secondaries
    lik[i] = r[fi[i]-pi[i]];
    //Rprintf("= %.4f\n", lik[i]);
  }
  free(x);
  return lik;
}

// [[Rcpp::export(.g_dat)]]
Rcpp::NumericVector g_dat_R(int n, double sar){
  int size = g_ind(n, n, 0)+1;
  Rcpp::NumericVector r(size);
  g_dat(r.begin(),n,sar);
  return r;
}

// [[Rcpp::export(.g_ind)]]
int g_ind_R(int n, int k, int m){
  return g_ind(n,k,m);
}

// Compute probabilities for final infection numbers from 0 .. m for a household
// with pi primary infections, size pi+m
double *g(const int pi, const int m, double *g_data, const int maxn){
  // get generating function for secondary infections
  return g_data + g_ind(maxn, pi, m);
}

// max hh size n, primary infections k, potential secondary infections m
// TODO: figure out the polynomial equivalent
/*int g_ind(int n, int k, int m){
  int ki,mi;
  int r = 0;
  if(k+m > n) Rcpp::stop("require k+m <= n");
  for(ki = 0; ki <= k-1; ki++){
    for(mi = 0; mi <= (n-ki); mi++){
      r = r + mi + 1;
    }
  }
  for(mi = 0;mi <= m-1; mi++){
    r = r + mi + 1;
  }
  return r;
}*/

long g_mem(int max_n){
  return g_ind(max_n,max_n,0)+1;
}

// By wolfram alpha:
// sum_(i=0)^(k - 1)( sum_(j=0)^(n - i)(j + 1)) + sum_(j=0)^(m - 1)(j + 1)
// 1/6 k (k^2 - 3 k n - 6 k + 3 n^2 + 12 n + 11) + 1/2 m (m + 1)
// k*(k^2 - 3*k*n - 6*k + 3*n^2 + 12*n + 11)/6 + m*(m + 1)/2
long g_ind(int n, int k, int m){
  if(k+m > n) Rcpp::stop("require k+m <= n");
  if(n > 1200) Rcpp::stop("require n <= 1200"); // overflow happens somehwere around 1285
  long n2 = (long)n; // just in case
  long k2 = (long)k;
  long m2 = (long)m;
  long r = (k2*(k2*k2 - 3*k2*n2 - 6*k2 + 3*n2*n2 + 12*n2 + 11))/6 + (m2*(m2 + 1))/2;
  return r;
}

// n: max hh size
// q: parameter for secondary attack rate
void g_dat(double *x, int n, double q){
  double *r;
  double *r_prev;
  double *g1;
  int k,m;
  // size is last index + 1
  int size = g_ind(n,n,0)+1;
  // clear x
  memset(x,0,size*sizeof(double));
  for(k = 0; k <= n; k++){
    for(m = 0; m <= (n-k); m++){
      //Rprintf("\nk=%d,m=%d\n",k,m);
      // get pointer to current polynomial, of degree m
      r = x + g_ind(n,k,m);
      if(m==0){
        r[0] = 1;
      }else{
        //Rprintf("r_prev="); for(int jj=0;jj<=m-1;jj++) Rprintf("%f,",r_prev[jj]); Rprintf("\n");
        for(int i = 0; i < m; i++){
          if(k+i==0){
            r[i] = r[i] + r_prev[i];
          }else{
            g1 = x + g_ind(n,1,m-i-1);
            //Rprintf("g1="); for(int jj=0;jj<=m-i-1;jj++) Rprintf("%f,",g1[jj]); Rprintf("\n");
            double qki = pow(1-q,k+i);
            r[i] = r[i] + r_prev[i] * qki;
            // iterate over terms in g1
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