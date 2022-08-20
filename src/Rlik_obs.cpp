#include <Rcpp.h>
#include <math.h>
#include <algorithm>
#include <stdint.h>
#include "lik_obs.h"
#include "sec_genfun.h"
#include "utilities.h"

// [[Rcpp::export(.Rlprob_obs_hhs)]]
double Rlprob_obs_hhs(Rcpp::IntegerVector obs, Rcpp::NumericVector prim, Rcpp::IntegerVector hh_sizes, Rcpp::NumericVector q, int n, double Sn, double Sp, int T, int M){
  return lprob_obs_hhs(obs.begin(), prim.begin(), hh_sizes.begin(), q.begin(), n, Sn, Sp, T, M);
}

// [[Rcpp::export(.Rlprob_obs)]]
double Rlprob_obs(Rcpp::IntegerVector obs, Rcpp::NumericVector prim, int hh_size, double q, double Sn, double Sp, int T, int M){
  return lprob_obs(obs.begin(), prim.begin(), hh_size, q, Sn, Sp, T, M);
}

// [[Rcpp::export(.Rlprob_obs_pri)]]
double Rlprob_obs_pri(Rcpp::IntegerVector obs, Rcpp::IntegerVector pri, int hh_size, double q, double Sn, double Sp){
  Gdat gdat = lgdat(hh_size,q);
  double lprob = lprob_obs_pri(obs.begin(), pri.begin(), hh_size, gdat, log(Sn), log(Sp), log(1-Sn), log(1-Sp));
  gfree(gdat);
  return lprob;
}

// [[Rcpp::export(.Rlprob_sec_obs_pri)]]
double Rlprob_sec_obs_pri(int obs, int mis, int pri, int m, double q, double Sn, double Sp){
  Gdat gdat = lgdat(m+pri,q);
  double lprob = lprob_sec_obs_pri(obs, mis, pri, m, gdat, log(Sn), log(Sp), log(1-Sn), log(1-Sp));
  gfree(gdat);
  return lprob;
}
