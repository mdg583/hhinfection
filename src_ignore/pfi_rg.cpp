#include <Rcpp.h>
#include <math.h>
#include "pfi_rg.h"

// pfi.cpp is a simple implementation of a few household infection models for Glasgow ODL project 2022
// by Matthew Gibson

// Two functions are exported: clusters, and pfi_lik. These functions are called by pfi.R,
// and documentation of how to use is in that file.

// Note that by default R's rng is prepared using RNGScope, which takes a few milliseconds of overhead
// [[Rcpp::export(.clusters,rng=false)]]
Rcpp::LogicalVector clusters_R(Rcpp::IntegerVector Ms, int n, Rcpp::LogicalVector xs) {
  int m = Ms.size() / (n*n);
  Rcpp::LogicalVector result(n*m);

  // for each Matrix & initial infected set
  for(int i = 0; i < m; i++){
    // get Matrix
    int *M = Ms.begin() + n*n*i;
    // get initial infected
    int *x = xs.begin() + n*i;
    // get visited nodes
    int *v = result.begin() + n*i;
    // For each infected, update nodes
    for(int j=0; j<n; j++){
      if(x[j]) cluster(M,n,j,v);
    }
  }
  return result;
}

// [[Rcpp::export(.pfi_lik_rg)]]
Rcpp::NumericVector pfi_lik_rg_R(Rcpp::IntegerVector pi,
                               Rcpp::IntegerVector fi,
                               Rcpp::IntegerVector ni,
                               double sar,
                               int iterations){
  int m = pi.size();
  Rcpp::NumericVector lik(m);
  // Each line of the vectors corresponds to a single household secondary infection probability calculation.
  for(int i = 0; i < m; i++){
    lik[i] = pfi_rg(fi[i], pi[i], ni[i], sar, iterations);
  }
  return lik;
}

// Given an adjacency matrix M (n x n), and node i, update v a
// logical array of length n to indicate nodes connected to i
void cluster(int *M, int n, int i, int *v){
  // Mark i as visited
  v[i] = 1;
  // Get adjacent nodes
  int *a = M + n*i;
  // For each node, visit if not visited
  for(int j = 0; j < n; j++){
    if(a[j]==1 && v[j]==0){
      cluster(M,n,j,v); // Updates the visited array
    }
  }
}

int rbern(double p){ return R::runif(0.0,1.0) < p; }

// generate a random household infection graph of size n x n with infection probability p
void randomgraph(int *M, int n, double p){
  // Clear M
  memset(M,0,n*n*sizeof(int));
  // iterate over upper tri
  for(int i=0;i<(n-1);i++){
    for(int j=i+1;j<n;j++){
      // create pair of edges with prob p
      if(rbern(p)){
        M[i*n+j] = 1;
        M[j*n+i] = 1;
      }
    }
  }
}

// Compute probability for a specific number of final infections (fi) based on a given number
// of primary infections (pi) and secondary attack rate (sar) in a household of size n
double pfi_rg(const int fi, const int pi, const int n, const double sar, const int iterations){
  int M[n*n]; // infection matrix
  int v[n];   // final infections

  // known probabilities
  if(fi > n) return(0.0);
  if(pi > fi) return(0.0);
  if(pi == 0 && fi > 0) return(0.0);
  if(pi == 0 && fi == 0) return(1.0);
  if(pi == n && fi == n) return(1.0);

  int countfi = 0;
  for(int i = 0; i < iterations; i++){
    // Create infection matrix
    randomgraph(M,n,sar);

    // clear v
    memset(v,0,n*sizeof(int));

    // Infect first pi individuals. For each of these, call cluster to update final infections  v
    for(int j=0; j<pi; j++){
      cluster(M,n,j,v);
    }

    // count final infections
    int thisfi = 0;
    for(int j = 0; j < n; j++){
      thisfi += v[j];
    }
    // if final infections == fi, increment countfi
    if(thisfi == fi) countfi++;
  }
  return ((double)countfi)/((double)iterations);
}
