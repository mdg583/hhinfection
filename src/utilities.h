#ifndef UTILITIES_H
#define UTILITIES_H

#include <Rcpp.h>
#include <math.h>

//' Multiply an integer by a double
//'
//' This function sets the result to 0 for all n=0, even for non-finite values.
//' The goal is to minimize operations, especially avoiding a jump operation.
//' This helps with performance but may not be much. Motivation is because originally
//' the code would iterate over all possible secondary infections.
//'
//' @param n integer to multiply
//' @param x floating-point value to multiply
//' @return result of n * x
static inline double mul(int n, double x){
	int64_t n_zero,nx2;
	double nx,nx3;
	nx = n*x;
	// memcpy used for syntactically-correct cast between float and int.
	// This is for the compiler and shouldn't result in any extra operation.
	memcpy(&nx2, &nx, sizeof nx);
	n_zero = 0 - (int64_t)(n!=0); // bitwise all 1's when n != 0
	nx2 = nx2 & n_zero; // If n==0, clear the bits of nx
	memcpy(&nx3, &nx2, sizeof nx2);
	return nx3;
}

// Binomial Coefficients

//' Get a binomial coefficient
//'
//' @param n Number of items from which to choose
//' @param r Number of item to choose
//' @return the number of ways of choosing r items from n items
static long long choose(int n, int r){
	if(2*r > n) r = n-r; // swap to easier case
	if(r < 0) Rcpp::stop("n choose r: r must be in 0 ... n\n");
	if(r==0) return 1; // special case: 0 choose 0 = 1
	if(r==1) return n;
	long long res = n-r+1; // n-r+1 choose 1 = n-r+1
	// (n choose r) = (n-1 choose r-1) * (n)/(r)
	// (n-r+i choose i) = (n-r+i-1 choose i-1) * (n-r+i)/(i)
	for(int i = 2; i <= r; i++){
		res = (res * (n-r+i))/i;
	}
	return res;
}

//' Get the binomial coefficient for an incremented value for r
//'
//' @param prev Previous binomial coefficient, for n choose r-1
//' @param n Number of items from which to choose
//' @param new_r Number of item to choose
//' @return the number of ways of choosing new_r items from n items
static inline long long choose_rinc(long long prev, int n, int new_r){
	// (n choose r) = (n choose r-1)(n-r+1)/(r)
	return (prev * (n - new_r + 1)) / new_r;
}

//' Get the binomial coefficient for an incremented value for r
//'
//' @param prev Previous binomial coefficient, for n choose r+1
//' @param n Number of items from which to choose
//' @param new_r Number of item to choose
//' @return The number of ways of choosing new_r items from n items
static inline long long choose_rdec(long long prev, int n, int new_r){
	// (n choose r) = (n choose r+1)(r+1)/(n-r)
	return (prev * (new_r + 1)) / (n - new_r);
}

#endif /* UTILITIES_H */
