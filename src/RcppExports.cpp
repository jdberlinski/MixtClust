// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// to_array
arma::cube to_array(NumericVector sigmas, int K, int p);
RcppExport SEXP _MixtClust_to_array(SEXP sigmasSEXP, SEXP KSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type sigmas(sigmasSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(to_array(sigmas, K, p));
    return rcpp_result_gen;
END_RCPP
}
// fix_var
arma::mat fix_var(arma::mat sigma, double tol);
RcppExport SEXP _MixtClust_fix_var(SEXP sigmaSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(fix_var(sigma, tol));
    return rcpp_result_gen;
END_RCPP
}
// mahalanobis
arma::vec mahalanobis(arma::mat x, arma::vec mu, arma::mat sigma, bool ischol);
RcppExport SEXP _MixtClust_mahalanobis(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP ischolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type ischol(ischolSEXP);
    rcpp_result_gen = Rcpp::wrap(mahalanobis(x, mu, sigma, ischol));
    return rcpp_result_gen;
END_RCPP
}
// dMVT
arma::vec dMVT(arma::mat x, arma::vec mu, arma::mat sigma, double nu, bool logans, bool ischol);
RcppExport SEXP _MixtClust_dMVT(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP nuSEXP, SEXP logansSEXP, SEXP ischolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< bool >::type logans(logansSEXP);
    Rcpp::traits::input_parameter< bool >::type ischol(ischolSEXP);
    rcpp_result_gen = Rcpp::wrap(dMVT(x, mu, sigma, nu, logans, ischol));
    return rcpp_result_gen;
END_RCPP
}
// h
arma::vec h(arma::mat x, arma::vec mu, arma::mat sigma, double nu, arma::vec grp, arma::umat Ru);
RcppExport SEXP _MixtClust_h(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP nuSEXP, SEXP grpSEXP, SEXP RuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type grp(grpSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type Ru(RuSEXP);
    rcpp_result_gen = Rcpp::wrap(h(x, mu, sigma, nu, grp, Ru));
    return rcpp_result_gen;
END_RCPP
}
// up_Z
arma::mat up_Z(arma::mat x, arma::mat mus, NumericVector sigmas, arma::vec nus, arma::vec pis, arma::vec grp, arma::umat Ru);
RcppExport SEXP _MixtClust_up_Z(SEXP xSEXP, SEXP musSEXP, SEXP sigmasSEXP, SEXP nusSEXP, SEXP pisSEXP, SEXP grpSEXP, SEXP RuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mus(musSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigmas(sigmasSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nus(nusSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pis(pisSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type grp(grpSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type Ru(RuSEXP);
    rcpp_result_gen = Rcpp::wrap(up_Z(x, mus, sigmas, nus, pis, grp, Ru));
    return rcpp_result_gen;
END_RCPP
}
// up_W
arma::mat up_W(arma::mat x, arma::mat mus, NumericVector sigmas, arma::vec nus, arma::vec grp, arma::umat Ru);
RcppExport SEXP _MixtClust_up_W(SEXP xSEXP, SEXP musSEXP, SEXP sigmasSEXP, SEXP nusSEXP, SEXP grpSEXP, SEXP RuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mus(musSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigmas(sigmasSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nus(nusSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type grp(grpSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type Ru(RuSEXP);
    rcpp_result_gen = Rcpp::wrap(up_W(x, mus, sigmas, nus, grp, Ru));
    return rcpp_result_gen;
END_RCPP
}
// up_pi
arma::vec up_pi(arma::mat z);
RcppExport SEXP _MixtClust_up_pi(SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(up_pi(z));
    return rcpp_result_gen;
END_RCPP
}
// up_mu
arma::mat up_mu(arma::mat x, arma::mat z, arma::mat w, arma::mat A);
RcppExport SEXP _MixtClust_up_mu(SEXP xSEXP, SEXP zSEXP, SEXP wSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(up_mu(x, z, w, A));
    return rcpp_result_gen;
END_RCPP
}
// up_Sigma
arma::cube up_Sigma(arma::mat x, arma::mat z, arma::mat w, arma::mat mus, arma::mat A, String constr);
RcppExport SEXP _MixtClust_up_Sigma(SEXP xSEXP, SEXP zSEXP, SEXP wSEXP, SEXP musSEXP, SEXP ASEXP, SEXP constrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mus(musSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< String >::type constr(constrSEXP);
    rcpp_result_gen = Rcpp::wrap(up_Sigma(x, z, w, mus, A, constr));
    return rcpp_result_gen;
END_RCPP
}
// up_nu
NumericVector up_nu(NumericMatrix z, NumericMatrix w, NumericVector nus, NumericVector ps, bool constr, bool approx);
RcppExport SEXP _MixtClust_up_nu(SEXP zSEXP, SEXP wSEXP, SEXP nusSEXP, SEXP psSEXP, SEXP constrSEXP, SEXP approxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w(wSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nus(nusSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ps(psSEXP);
    Rcpp::traits::input_parameter< bool >::type constr(constrSEXP);
    Rcpp::traits::input_parameter< bool >::type approx(approxSEXP);
    rcpp_result_gen = Rcpp::wrap(up_nu(z, w, nus, ps, constr, approx));
    return rcpp_result_gen;
END_RCPP
}
// SOiOEOOk
arma::cube SOiOEOOk(arma::mat sigma, arma::umat Ru);
RcppExport SEXP _MixtClust_SOiOEOOk(SEXP sigmaSEXP, SEXP RuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type Ru(RuSEXP);
    rcpp_result_gen = Rcpp::wrap(SOiOEOOk(sigma, Ru));
    return rcpp_result_gen;
END_RCPP
}
// xhatk
arma::mat xhatk(arma::mat x, arma::vec mu, arma::vec grp, int M, NumericVector SOiOEOOk);
RcppExport SEXP _MixtClust_xhatk(SEXP xSEXP, SEXP muSEXP, SEXP grpSEXP, SEXP MSEXP, SEXP SOiOEOOkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type grp(grpSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type SOiOEOOk(SOiOEOOkSEXP);
    rcpp_result_gen = Rcpp::wrap(xhatk(x, mu, grp, M, SOiOEOOk));
    return rcpp_result_gen;
END_RCPP
}
// up_mu_Lin
arma::mat up_mu_Lin(int p, arma::mat z, arma::mat w, ListOf<NumericMatrix> xhat);
RcppExport SEXP _MixtClust_up_mu_Lin(SEXP pSEXP, SEXP zSEXP, SEXP wSEXP, SEXP xhatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type w(wSEXP);
    Rcpp::traits::input_parameter< ListOf<NumericMatrix> >::type xhat(xhatSEXP);
    rcpp_result_gen = Rcpp::wrap(up_mu_Lin(p, z, w, xhat));
    return rcpp_result_gen;
END_RCPP
}
// up_Sigmak_Lin
arma::mat up_Sigmak_Lin(int M, arma::vec zk, arma::vec wk, arma::vec mu, arma::mat sigma, arma::mat xhatk, arma::vec grp, NumericVector SOiOEOOk);
RcppExport SEXP _MixtClust_up_Sigmak_Lin(SEXP MSEXP, SEXP zkSEXP, SEXP wkSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP xhatkSEXP, SEXP grpSEXP, SEXP SOiOEOOkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type zk(zkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type wk(wkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xhatk(xhatkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type grp(grpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type SOiOEOOk(SOiOEOOkSEXP);
    rcpp_result_gen = Rcpp::wrap(up_Sigmak_Lin(M, zk, wk, mu, sigma, xhatk, grp, SOiOEOOk));
    return rcpp_result_gen;
END_RCPP
}
// Q2
double Q2(arma::mat x, arma::mat z, arma::mat w, NumericVector sigmas, arma::mat mus, arma::vec grp, arma::umat Ru);
RcppExport SEXP _MixtClust_Q2(SEXP xSEXP, SEXP zSEXP, SEXP wSEXP, SEXP sigmasSEXP, SEXP musSEXP, SEXP grpSEXP, SEXP RuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type w(wSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigmas(sigmasSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mus(musSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type grp(grpSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type Ru(RuSEXP);
    rcpp_result_gen = Rcpp::wrap(Q2(x, z, w, sigmas, mus, grp, Ru));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MixtClust_to_array", (DL_FUNC) &_MixtClust_to_array, 3},
    {"_MixtClust_fix_var", (DL_FUNC) &_MixtClust_fix_var, 2},
    {"_MixtClust_mahalanobis", (DL_FUNC) &_MixtClust_mahalanobis, 4},
    {"_MixtClust_dMVT", (DL_FUNC) &_MixtClust_dMVT, 6},
    {"_MixtClust_h", (DL_FUNC) &_MixtClust_h, 6},
    {"_MixtClust_up_Z", (DL_FUNC) &_MixtClust_up_Z, 7},
    {"_MixtClust_up_W", (DL_FUNC) &_MixtClust_up_W, 6},
    {"_MixtClust_up_pi", (DL_FUNC) &_MixtClust_up_pi, 1},
    {"_MixtClust_up_mu", (DL_FUNC) &_MixtClust_up_mu, 4},
    {"_MixtClust_up_Sigma", (DL_FUNC) &_MixtClust_up_Sigma, 6},
    {"_MixtClust_up_nu", (DL_FUNC) &_MixtClust_up_nu, 6},
    {"_MixtClust_SOiOEOOk", (DL_FUNC) &_MixtClust_SOiOEOOk, 2},
    {"_MixtClust_xhatk", (DL_FUNC) &_MixtClust_xhatk, 5},
    {"_MixtClust_up_mu_Lin", (DL_FUNC) &_MixtClust_up_mu_Lin, 4},
    {"_MixtClust_up_Sigmak_Lin", (DL_FUNC) &_MixtClust_up_Sigmak_Lin, 8},
    {"_MixtClust_Q2", (DL_FUNC) &_MixtClust_Q2, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_MixtClust(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
