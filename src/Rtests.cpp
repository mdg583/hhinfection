#include <Rcpp.h>
#include <math.h>
#include <limits.h>
#include <chrono>
#include "lik_obs.h"
#include "utilities.h"

// // [[Rcpp::export(.test_bcoef)]]
// int test_bcoef(){
//   long coefs[11];
//   bcoefs(coefs, 10);
//   long expected[] = {1,10,45,120,210,252,210,120,45,10,1};
//   int equal = 1;
//   for(int i = 0; i <= 10; i++){
//     Rcpp::Rcout << coefs[i] << " ";
//     if(coefs[i] != expected[i]) equal = 0;
//   }
//   Rcpp::Rcout << "\n";
//   return(equal);
// }

// // [[Rcpp::export(.test_count_bits)]]
// int test_count_bits(){
//   int pass = 1;

//   int64_t test[] = {4,7,12,16,19,0,(1<<13)-1,(1LL<<35)-1,(1LL<<62)-1};
//   int64_t expect[] = {1,3,2,1,3,0,13,35,62};

//   for(int i = 0; i < 9; i++){
//     Rcpp::Rcout << "count_bits(" << test[i] << ") = " << count_bits(test[i]) << "\n";
//     if(count_bits(test[i]) != expect[i]) pass = 0;
//   }
//   return(pass);
// }