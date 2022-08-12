#include <stdint.h>

// bitcount lookup table for 8 bit integers
// from https://stackoverflow.com/questions/51387998/count-bits-1-on-an-integer-as-fast-as-gcc-builtin-popcountint#answer-51388543
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

/**
 * Get n+1 binomial coefficients from (n choose 0) to (n choose n)
 *
 * @param coefs pointer to array of size n+1 for storing binomial coefficients
 * @param n
 */
void bcoefs(long *coefs, long n);

/**
 * Get the number of 1-bits in a 64-bit integer
 *
 * @param int64_t x integer for which to count bits
 *
 * @return int number of 1-bits in x
 */
int count_bits(int64_t x);

/**
 * Probability of observed household outcome given actual outcome (log)
 *
 * @param int64_t obs observed household cases
 * @param int64_t act actual (assumed actual) household cases
 * @param int64_t mis unobserved individuals
 * @param int n number of individuals in household
 * @param double l_Sn log of sensitivity
 * @param double l_Sp log of specificity
 * @param double l_1mSn log of 1 minus sensitivity
 * @param double l_1mSp log of 1 minus specificity
 *
 * @return double probability of observation
 */
inline double lprob_obs_act(int64_t obs, int64_t act, int64_t mis, int n,
                            double l_Sn, double l_Sp, double l_1mSn, double l_1mSp);

/**
 * Probability of observed household outcome given primary infections (log)
 *
 * Probability is computed by iterating over all possible secondary infections
 *
 * @param int64_t obs observed household cases
 * @param int64_t pri primary (assumed primary) household cases
 * @param int64_t mis unobserved individuals
 * @param int hh_size number of individuals in household
 * @param double *lg_data precomputed secondary infection generating functions
 * @param int max_hh maximum household size, used to describe structure of lg_data
 * @param double l_Sn log of sensitivity
 * @param double l_Sp log of specificity
 * @param double l_1mSn log of 1 minus sensitivity
 * @param double l_1mSp log of 1 minus specificity
 *
 * @return double probability of observation
 */
double lprob_obs_pri_it(int64_t obs, int64_t pri, int64_t mis, int hh_size, double *lg_data, int max_hh,
                        double l_Sn, double l_Sp, double l_1mSn, double l_1mSp);

/**
 * Probability of observed household outcome given primary infections (log)
 *
 * Probability is estimated by computing probability of observation based on simulated secondary infections
 *
 * @param int64_t obs observed household cases
 * @param int64_t pri primary (assumed primary) household cases
 * @param int64_t mis unobserved individuals
 * @param int hh_size number of individuals in household
 * @param double q secondary attack rate
 * @param double l_Sn log of sensitivity
 * @param double l_Sp log of specificity
 * @param double l_1mSn log of 1 minus sensitivity
 * @param double l_1mSp log of 1 minus specificity
 * @param double *pbuf2 space for M2 computed probabilities (double)
 * @param int M2 number of secondary simulations to average over
 *
 * @return double probability of observation
 */
double lprob_obs_pri_sim(int64_t obs, int64_t pri, int64_t mis, const int hh_size, const double q,
                         const double l_Sn, const double l_Sp, const double l_1mSn, const double l_1mSp,
                         double *pbuf2, const int M2);
/**
 * Probability of household primary infection
 *
 * @param int64_t pri primary (assumed primary) household cases
 * @param double *lprim log probabilities for primary infection
 * @param double *l_1mprim log of 1 minus probabilities for primary infection
 * @param int n number of individuals in household
 *
 * @return double probability of primary infection
 */
double lprob_pri(int64_t pri, double *lprim, double *l_1mprim, int n);

/**
 * Generate a random primary infection
 *
 * @param double *p probabilities of primary infection
 * @param int n number of individuals in household
 *
 * @return int64_t integer whose binary digits correspond to primary infections
 */
int64_t rbool(double *p, int n);

/**
 * Probability of observed household outcome
 *
 * Probability is computed by iterating over all possible primary infections
 *
 * @param int64_t obsv observed household cases
 * @param int64_t misv missed household cases
 * @param double *lprim log probabilities of primary infection
 * @param double *l_1mprim log of 1 minus probabilities of primary infection
 * @param int hh_size number of individuals in household
 * @param double *lg_data precomputed secondary infection generating functions
 * @param int max_hh maximum household size, used to describe structure of lg_data
 * @param double q secondary attack rate
 * @param double l_Sn log of sensitivity
 * @param double l_Sp log of specificity
 * @param double l_1mSn log of 1 minus sensitivity
 * @param double l_1mSp log of 1 minus specificity
 * @param double *pbuf1 space for M1 computed probabilities (double)
 * @param double *pbuf2 space for M2 computed probabilities (double)
 * @param int M2 number of secondary simulations to average over
 * @param int T2 household size theshold for use of simulation for secondary infection probabilities
 *
 * @return double probability of observation
 */
double lprob_obs_it(int64_t obsv, int64_t misv, double *lprim, double *l_1mprim, const int hh_size,
                    double *lg_data, const int max_hh, double q,
                    const double l_Sn, const double l_Sp, const double l_1mSn, const double l_1mSp,
                    double *pbuf1, double *pbuf2, const int M2, int T2);

/**
 * Probability of observed household outcome
 *
 * Probability is computed by simulating primary infections
 *
 * @param int64_t obsv observed household cases
 * @param int64_t misv missed household cases
 * @param double *prim probabilities of primary infection
 * @param int hh_size number of individuals in household
 * @param double *lg_data precomputed secondary infection generating functions
 * @param int max_hh maximum household size, used to describe structure of lg_data
 * @param double q secondary attack rate
 * @param double l_Sn log of sensitivity
 * @param double l_Sp log of specificity
 * @param double l_1mSn log of 1 minus sensitivity
 * @param double l_1mSp log of 1 minus specificity
 * @param double *pbuf1 space for M1 computed probabilities (double)
 * @param double *pbuf2 space for M2 computed probabilities (double)
 * @param int M1 number of primary simulations to average over
 * @param int M2 number of secondary simulations to average over
 * @param int T2 household size threshold for use of simulation for secondary infection probabilities
 *
 * @return double probability of observation
 */
double lprob_obs_sim(int64_t obsv, int64_t misv, double *prim, const int hh_size,
                     double *lg_data, const int max_hh, double q,
                     const double l_Sn, const double l_Sp, const double l_1mSn, const double l_1mSp,
                     double *pbuf1, double *pbuf2, const int M1, const int M2, int T2);

/**
 * Probability of observed household outcome
 *
 * Probability is computed by simulating primary infections
 *
 * @param int *obs observed household cases. 1=infected,-1=not observed (missed)
 * @param double *prim probabilities of primary infection
 * @param int hh_size number of individuals in household
 * @param double *lg_data pre-computed secondary infection generating functions
 * @param int max_hh maximum household size, used to describe structure of lg_data
 * @param double q secondary attack rate
 * @param double l_Sn log of sensitivity
 * @param double l_Sp log of specificity
 * @param double l_1mSn log of 1 minus sensitivity
 * @param double l_1mSp log of 1 minus specificity
 * @param double *pbuf1 space for M1 computed probabilities (double)
 * @param double *pbuf2 space for M2 computed probabilities (double)
 * @param int M1 number of primary simulations to average over
 * @param int M2 number of secondary simulations to average over
 * @param int T1 household size threshold for use of simulation for primary infection probabilities
 * @param int T2 household size threshold for use of simulation for secondary infection probabilities
 *
 * @return double probability of observation
 */
double lprob_obs(int *obs, double *prim, const int hh_size, double *lg_data, const int max_hh, const double q,
                 const double l_Sn, const double l_Sp, const double l_1mSn, const double l_1mSp,
                 double *pbuf1, double *pbuf2, const int M1, const int M2, const int T1, const int T2);

/**
 * Probability of observed outcomes across multiple households
 *
 * Probability is computed by simulating primary infections
 *
 * @param int *obs observed cases. 1=infected,-1=not observed (missed)
 * @param double *prim probabilities of primary infection
 * @param int *hh household identifier, non-decreasing sequenc
 * @param int n number of observed individuals
 * @param double q probability of secondary infection
 * @param double Sn test sensitivity
 * @param double Sp test specificity
 * @param int M1 number of iterations for primary infection simulation
 * @param int M2 number of iterations for secondary infection simulation
 * @param int T1 household size threshold for use of simulation for primary infection
 * @param int T2 household size threshold for use of simulation for secondary infection
 *
 * @return double probability of observation
 */
double lprob_obs_hhs(int *obs, double *prim, int *hh, const int n, const double q,
                     const double Sn, const double Sp,
                     const int M1, const int M2, const int T1, const int T2);
