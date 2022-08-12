#ifndef UTILITIES_H
#define UTILITIES_H

#include <Rcpp.h>
#include <math.h>

// Bit-Counting using a lookup table

// This is taken from a stackoverflow answer:
// https://stackoverflow.com/questions/51387998/count-bits-1-on-an-integer-as-fast-as-gcc-builtin-popcountint#answer-51388543
// Some compilers (GCC, MS, vlang) on 64-bit platforms have a function to do this using a single
// assembly instruction. But R may support some other platform.
const int bitcount[] = {
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8,
};

//' Count bits in 64-bit integer
//'
//' @param x A 64-bit integer
//' @return the number of 1-bits in x
static inline int count_bits(int64_t x){
	int count = 0;
	unsigned char *ptr = (unsigned char *)&x;
		for (int i=0;i<(int)sizeof(int64_t);i++) {
				count += bitcount[ptr[i]];
		}
		return count;
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
	// n-r+1 choose 1 = n-r+1
	long long res = n-r+1;
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

//' Get all binomial coefficients for a given value of n
//'
//' @param coefs Pointer to space for n+1 binomial coefficients.
//' @param n Number of items from which to choose
static void bcoefs(long *coefs, long n){
	coefs[0] = 1;
	coefs[1] = n;
	for(int i = 2; i <= n/2; i++){
		coefs[i] = choose_rinc(coefs[i-1], n, i);
	}
	for(int i = (n+1)/2; i <= n; i++){
		coefs[i] = coefs[n-i];
	}
}


// Debugging Functions

//' Print the elements of the vector v of length n, including final newline
static void print_vec(int *v, int n){
	for(int i = 0; i < n; i++){
		Rprintf("%d ",v[i]);
	}
	Rprintf("\n");
}

//' Print the bit-values of the 64-bit integer v, up to n bits, without final newline
static void print_bvec(int64_t v, int n){
	for(int i = 0; i < n; i++){
		Rprintf("%d ",(v >> i) & 1);
	}
}

//' Print the floating-point elements of the vector v of length n, including final newline. (2 decimal places)
static void print_fvec(double *v, int n){
	for(int i = 0; i < n; i++){
		Rprintf("%.2f ",v[i]);
	}
	Rprintf("\n");
}


//' Multiply an integer by a double
//'
//' This function sets the result to 0 for all n=0, even for non-finite values.
//' The goal is to minimize operations, especially avoiding a jump operation. Testing shows that
//' this helps with performance but may not be much.
//' @param n integer to multiply by x
//' @param x floating-point value to multiply by n
//' @return result of n * x
static inline double mul(int n, double x){
	int64_t n_zero,nx2;
	double nx,nx3;
	// bitwise all 1's when n != 0
	n_zero = 0 - (int64_t)(n!=0);
	nx = n*x;
	// If n==0, clear the bits of nx
	// memcpy used for syntactically-correct casts for float to/from int
	// The compiler will just change how it handles the data from float to int
	memcpy(&nx2, &nx, sizeof nx); // (*(long long*)&nx)
	nx2 = nx2 & n_zero; // set bits to 0 if n==0
	memcpy(&nx3, &nx2, sizeof nx2); // (*(long long*)&nx)
	return nx3;
}
#endif /* UTILITIES_H */
