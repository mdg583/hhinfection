#include "sec_genfun.h"

// Function comments are in source code

double lprob_obs_sec(double *pbuf, int m, int act, int obs, double lfact, double l_Sn, double l_Sp, double l_1mSn, double l_1mSp);
double lprob_sec_obs_pri(int obs, int mis, int pri, int m, Gdat gdat, double l_Sn, double l_Sp, double l_1mSn, double l_1mSp);
double lprob_obs_pri(int* obs, int* pri, const int hh_size, Gdat gdat, const double l_Sn, const double l_Sp, const double l_1mSn, const double l_1mSp);
double lprob_pri(int *pri, double *lprim, double *l_1mprim, int n);
void pri_sample(int *pri, double *p, int n);
double lprob_obs(int* obs, double *prim, const int hh_size, double q, const double Sn, const double Sp, const int T, int M);
double lprob_obs_hhs(int *obs, double *prim, int *hh_sizes, double *q, const int n, const double Sn, const double Sp, const int T, const int M);
