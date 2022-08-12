#include <Rcpp.h>
#include <math.h>
#include <algorithm>
#include <stdint.h>
#include "sec_genfun.h"

// [[Rcpp::export(.Rgdat)]]
Rcpp::NumericVector Rgdat(int max_hh, double q){
  // Generate secondary infection probabilities
  Gdat gd = gdat(max_hh,q);
  // copy into numeric vector
  int n = gsize(gd);
  Rcpp::NumericVector dvec(n);
  double *darr = dvec.begin();
  for(int i = 0; i < n; i++) {
    darr[i] = gd.data[i];
  }
  gfree(gd);
  return dvec;
}

// [[Rcpp::export(.Rg)]]
Rcpp::NumericVector Rg(Rcpp::NumericVector gdat, int max_hh, int k, int m){
  Gdat gdat2;
  // load arguements into Gdat object. Be careful not to free this object
  gdat2.max_hh = max_hh;
  gdat2.data = gdat.begin();
  // get polynomial and store in numeric vector
  Rcpp::NumericVector ret(m+1);
  double *r = g(gdat2,k,m);
  for(int i = 0; i <= m; i++){
    ret[i] = r[i];
  }
  return ret;
}
