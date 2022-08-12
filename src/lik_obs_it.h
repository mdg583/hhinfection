#include "sec_genfun.h"

inline double lprob_obs_act(int64_t obs, int64_t act, int64_t mis, int n, double l_Sn, double l_Sp, double l_1mSn, double l_1mSp);
double lprob_sec_obs_pri_it(int64_t obs, int64_t mis, int k, int m, Gdat gdat, double l_Sn, double l_Sp, double l_1mSn, double l_1mSp);
double lprob_obs_pri_it(int* obs, int* pri, const int hh_size, Gdat gdat, const double l_Sn, const double l_Sp, const double l_1mSn, const double l_1mSp);
double lprob_obs_it(int* obs, double *prim, const int hh_size, double q, const double Sn, const double Sp, const int T, int M);
double lprob_obs_hhs_it(int *obs, double *prim, int *hh_sizes, double *q, const int n, const double Sn, const double Sp, const int T, const int M);