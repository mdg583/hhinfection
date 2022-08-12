// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Rlprob_obs_hhs
double Rlprob_obs_hhs(Rcpp::IntegerVector obs, Rcpp::NumericVector prim, Rcpp::IntegerVector hh_sizes, Rcpp::NumericVector q, int n, double Sn, double Sp, int T, int M);
RcppExport SEXP _hhinfection_Rlprob_obs_hhs(SEXP obsSEXP, SEXP primSEXP, SEXP hh_sizesSEXP, SEXP qSEXP, SEXP nSEXP, SEXP SnSEXP, SEXP SpSEXP, SEXP TSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type prim(primSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type hh_sizes(hh_sizesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type Sn(SnSEXP);
    Rcpp::traits::input_parameter< double >::type Sp(SpSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(Rlprob_obs_hhs(obs, prim, hh_sizes, q, n, Sn, Sp, T, M));
    return rcpp_result_gen;
END_RCPP
}
// Rlprob_obs
double Rlprob_obs(Rcpp::IntegerVector obs, Rcpp::NumericVector prim, int hh_size, double q, double Sn, double Sp, int T, int M);
RcppExport SEXP _hhinfection_Rlprob_obs(SEXP obsSEXP, SEXP primSEXP, SEXP hh_sizeSEXP, SEXP qSEXP, SEXP SnSEXP, SEXP SpSEXP, SEXP TSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type prim(primSEXP);
    Rcpp::traits::input_parameter< int >::type hh_size(hh_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type Sn(SnSEXP);
    Rcpp::traits::input_parameter< double >::type Sp(SpSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(Rlprob_obs(obs, prim, hh_size, q, Sn, Sp, T, M));
    return rcpp_result_gen;
END_RCPP
}
// Rlprob_obs_pri
double Rlprob_obs_pri(Rcpp::IntegerVector obs, Rcpp::IntegerVector pri, int hh_size, double q, double Sn, double Sp);
RcppExport SEXP _hhinfection_Rlprob_obs_pri(SEXP obsSEXP, SEXP priSEXP, SEXP hh_sizeSEXP, SEXP qSEXP, SEXP SnSEXP, SEXP SpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type pri(priSEXP);
    Rcpp::traits::input_parameter< int >::type hh_size(hh_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type Sn(SnSEXP);
    Rcpp::traits::input_parameter< double >::type Sp(SpSEXP);
    rcpp_result_gen = Rcpp::wrap(Rlprob_obs_pri(obs, pri, hh_size, q, Sn, Sp));
    return rcpp_result_gen;
END_RCPP
}
// Rlprob_sec_obs_pri
double Rlprob_sec_obs_pri(int obs, int mis, int pri, int m, double q, double Sn, double Sp);
RcppExport SEXP _hhinfection_Rlprob_sec_obs_pri(SEXP obsSEXP, SEXP misSEXP, SEXP priSEXP, SEXP mSEXP, SEXP qSEXP, SEXP SnSEXP, SEXP SpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< int >::type mis(misSEXP);
    Rcpp::traits::input_parameter< int >::type pri(priSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type Sn(SnSEXP);
    Rcpp::traits::input_parameter< double >::type Sp(SpSEXP);
    rcpp_result_gen = Rcpp::wrap(Rlprob_sec_obs_pri(obs, mis, pri, m, q, Sn, Sp));
    return rcpp_result_gen;
END_RCPP
}
// Rlprob_obs_hhs_it
double Rlprob_obs_hhs_it(Rcpp::IntegerVector obs, Rcpp::NumericVector prim, Rcpp::IntegerVector hh_sizes, Rcpp::NumericVector q, int n, double Sn, double Sp, int T, int M);
RcppExport SEXP _hhinfection_Rlprob_obs_hhs_it(SEXP obsSEXP, SEXP primSEXP, SEXP hh_sizesSEXP, SEXP qSEXP, SEXP nSEXP, SEXP SnSEXP, SEXP SpSEXP, SEXP TSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type prim(primSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type hh_sizes(hh_sizesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type Sn(SnSEXP);
    Rcpp::traits::input_parameter< double >::type Sp(SpSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(Rlprob_obs_hhs_it(obs, prim, hh_sizes, q, n, Sn, Sp, T, M));
    return rcpp_result_gen;
END_RCPP
}
// Rlprob_obs_it
double Rlprob_obs_it(Rcpp::IntegerVector obs, Rcpp::NumericVector prim, int hh_size, double q, double Sn, double Sp, int T, int M);
RcppExport SEXP _hhinfection_Rlprob_obs_it(SEXP obsSEXP, SEXP primSEXP, SEXP hh_sizeSEXP, SEXP qSEXP, SEXP SnSEXP, SEXP SpSEXP, SEXP TSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type prim(primSEXP);
    Rcpp::traits::input_parameter< int >::type hh_size(hh_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type Sn(SnSEXP);
    Rcpp::traits::input_parameter< double >::type Sp(SpSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(Rlprob_obs_it(obs, prim, hh_size, q, Sn, Sp, T, M));
    return rcpp_result_gen;
END_RCPP
}
// Rlprob_obs_pri_it
double Rlprob_obs_pri_it(Rcpp::IntegerVector obs, Rcpp::IntegerVector pri, int m, double q, double Sn, double Sp);
RcppExport SEXP _hhinfection_Rlprob_obs_pri_it(SEXP obsSEXP, SEXP priSEXP, SEXP mSEXP, SEXP qSEXP, SEXP SnSEXP, SEXP SpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type pri(priSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type Sn(SnSEXP);
    Rcpp::traits::input_parameter< double >::type Sp(SpSEXP);
    rcpp_result_gen = Rcpp::wrap(Rlprob_obs_pri_it(obs, pri, m, q, Sn, Sp));
    return rcpp_result_gen;
END_RCPP
}
// Rlprob_sec_obs_pri_it
double Rlprob_sec_obs_pri_it(Rcpp::IntegerVector obs, int k, int m, double q, double Sn, double Sp);
RcppExport SEXP _hhinfection_Rlprob_sec_obs_pri_it(SEXP obsSEXP, SEXP kSEXP, SEXP mSEXP, SEXP qSEXP, SEXP SnSEXP, SEXP SpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type Sn(SnSEXP);
    Rcpp::traits::input_parameter< double >::type Sp(SpSEXP);
    rcpp_result_gen = Rcpp::wrap(Rlprob_sec_obs_pri_it(obs, k, m, q, Sn, Sp));
    return rcpp_result_gen;
END_RCPP
}
// Rgdat
Rcpp::NumericVector Rgdat(int max_hh, double q);
RcppExport SEXP _hhinfection_Rgdat(SEXP max_hhSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type max_hh(max_hhSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(Rgdat(max_hh, q));
    return rcpp_result_gen;
END_RCPP
}
// Rg
Rcpp::NumericVector Rg(Rcpp::NumericVector gdat, int max_hh, int k, int m);
RcppExport SEXP _hhinfection_Rg(SEXP gdatSEXP, SEXP max_hhSEXP, SEXP kSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type gdat(gdatSEXP);
    Rcpp::traits::input_parameter< int >::type max_hh(max_hhSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(Rg(gdat, max_hh, k, m));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_hhinfection_Rlprob_obs_hhs", (DL_FUNC) &_hhinfection_Rlprob_obs_hhs, 9},
    {"_hhinfection_Rlprob_obs", (DL_FUNC) &_hhinfection_Rlprob_obs, 8},
    {"_hhinfection_Rlprob_obs_pri", (DL_FUNC) &_hhinfection_Rlprob_obs_pri, 6},
    {"_hhinfection_Rlprob_sec_obs_pri", (DL_FUNC) &_hhinfection_Rlprob_sec_obs_pri, 7},
    {"_hhinfection_Rlprob_obs_hhs_it", (DL_FUNC) &_hhinfection_Rlprob_obs_hhs_it, 9},
    {"_hhinfection_Rlprob_obs_it", (DL_FUNC) &_hhinfection_Rlprob_obs_it, 8},
    {"_hhinfection_Rlprob_obs_pri_it", (DL_FUNC) &_hhinfection_Rlprob_obs_pri_it, 6},
    {"_hhinfection_Rlprob_sec_obs_pri_it", (DL_FUNC) &_hhinfection_Rlprob_sec_obs_pri_it, 6},
    {"_hhinfection_Rgdat", (DL_FUNC) &_hhinfection_Rgdat, 2},
    {"_hhinfection_Rg", (DL_FUNC) &_hhinfection_Rg, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_hhinfection(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
