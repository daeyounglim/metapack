// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// BMVMR_POCov
Rcpp::List BMVMR_POCov(const arma::mat& Outcome, const arma::mat& SD, const arma::mat& XCovariate, const arma::mat& WCovariate, const arma::uvec& Treat, const arma::uvec& Trial, const arma::vec& Npt, const double& c0, const double& dj0, const double& d0, const double& s0, const double& nu0, const arma::mat& Omega0, const arma::mat& Sigma0, const int& K, const int& T, const int& fmodel, const int& ndiscard, const int& nskip, const int& nkeep, const bool& verbose);
RcppExport SEXP _metapack_BMVMR_POCov(SEXP OutcomeSEXP, SEXP SDSEXP, SEXP XCovariateSEXP, SEXP WCovariateSEXP, SEXP TreatSEXP, SEXP TrialSEXP, SEXP NptSEXP, SEXP c0SEXP, SEXP dj0SEXP, SEXP d0SEXP, SEXP s0SEXP, SEXP nu0SEXP, SEXP Omega0SEXP, SEXP Sigma0SEXP, SEXP KSEXP, SEXP TSEXP, SEXP fmodelSEXP, SEXP ndiscardSEXP, SEXP nskipSEXP, SEXP nkeepSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Outcome(OutcomeSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type SD(SDSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type XCovariate(XCovariateSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type WCovariate(WCovariateSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type Treat(TreatSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type Trial(TrialSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Npt(NptSEXP);
    Rcpp::traits::input_parameter< const double& >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const double& >::type dj0(dj0SEXP);
    Rcpp::traits::input_parameter< const double& >::type d0(d0SEXP);
    Rcpp::traits::input_parameter< const double& >::type s0(s0SEXP);
    Rcpp::traits::input_parameter< const double& >::type nu0(nu0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Omega0(Omega0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sigma0(Sigma0SEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int& >::type T(TSEXP);
    Rcpp::traits::input_parameter< const int& >::type fmodel(fmodelSEXP);
    Rcpp::traits::input_parameter< const int& >::type ndiscard(ndiscardSEXP);
    Rcpp::traits::input_parameter< const int& >::type nskip(nskipSEXP);
    Rcpp::traits::input_parameter< const int& >::type nkeep(nkeepSEXP);
    Rcpp::traits::input_parameter< const bool& >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(BMVMR_POCov(Outcome, SD, XCovariate, WCovariate, Treat, Trial, Npt, c0, dj0, d0, s0, nu0, Omega0, Sigma0, K, T, fmodel, ndiscard, nskip, nkeep, verbose));
    return rcpp_result_gen;
END_RCPP
}
// BayesNMR
Rcpp::List BayesNMR(const arma::vec& y, const arma::vec& sd, const arma::mat& x, const arma::mat& z, const arma::uvec& ids, const arma::uvec& iarm, const arma::vec& npt, const double& nu, const double& c01_inv, const double& c02_inv, const int& K, const int& nT, const int& ndiscard, const int& nskip, const int& nkeep, const bool verbose, const arma::vec& beta_init, const arma::vec& phi_init, const arma::vec& sig2_init);
RcppExport SEXP _metapack_BayesNMR(SEXP ySEXP, SEXP sdSEXP, SEXP xSEXP, SEXP zSEXP, SEXP idsSEXP, SEXP iarmSEXP, SEXP nptSEXP, SEXP nuSEXP, SEXP c01_invSEXP, SEXP c02_invSEXP, SEXP KSEXP, SEXP nTSEXP, SEXP ndiscardSEXP, SEXP nskipSEXP, SEXP nkeepSEXP, SEXP verboseSEXP, SEXP beta_initSEXP, SEXP phi_initSEXP, SEXP sig2_initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ids(idsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type iarm(iarmSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type npt(nptSEXP);
    Rcpp::traits::input_parameter< const double& >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const double& >::type c01_inv(c01_invSEXP);
    Rcpp::traits::input_parameter< const double& >::type c02_inv(c02_invSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int& >::type nT(nTSEXP);
    Rcpp::traits::input_parameter< const int& >::type ndiscard(ndiscardSEXP);
    Rcpp::traits::input_parameter< const int& >::type nskip(nskipSEXP);
    Rcpp::traits::input_parameter< const int& >::type nkeep(nkeepSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta_init(beta_initSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type phi_init(phi_initSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sig2_init(sig2_initSEXP);
    rcpp_result_gen = Rcpp::wrap(BayesNMR(y, sd, x, z, ids, iarm, npt, nu, c01_inv, c02_inv, K, nT, ndiscard, nskip, nkeep, verbose, beta_init, phi_init, sig2_init));
    return rcpp_result_gen;
END_RCPP
}
// calc_modelfit_dic
Rcpp::List calc_modelfit_dic(const arma::vec& y, const arma::mat& x, const arma::mat& z, const arma::uvec& ids, const arma::uvec& iarm, const arma::vec& npt, const double& nu, const arma::mat& betas, const arma::mat& sig2s, const arma::mat& phis, const arma::mat& lams, const arma::cube& Rhos, const int& K, const int& nT, const int& nkeep, const bool verbose);
RcppExport SEXP _metapack_calc_modelfit_dic(SEXP ySEXP, SEXP xSEXP, SEXP zSEXP, SEXP idsSEXP, SEXP iarmSEXP, SEXP nptSEXP, SEXP nuSEXP, SEXP betasSEXP, SEXP sig2sSEXP, SEXP phisSEXP, SEXP lamsSEXP, SEXP RhosSEXP, SEXP KSEXP, SEXP nTSEXP, SEXP nkeepSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ids(idsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type iarm(iarmSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type npt(nptSEXP);
    Rcpp::traits::input_parameter< const double& >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type betas(betasSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sig2s(sig2sSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type phis(phisSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lams(lamsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Rhos(RhosSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int& >::type nT(nTSEXP);
    Rcpp::traits::input_parameter< const int& >::type nkeep(nkeepSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_modelfit_dic(y, x, z, ids, iarm, npt, nu, betas, sig2s, phis, lams, Rhos, K, nT, nkeep, verbose));
    return rcpp_result_gen;
END_RCPP
}
// veclinv
arma::mat veclinv(const arma::vec& v, const int& n);
RcppExport SEXP _metapack_veclinv(SEXP vSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(veclinv(v, n));
    return rcpp_result_gen;
END_RCPP
}
// pRho_to_Rho
arma::mat pRho_to_Rho(arma::mat& pRho);
RcppExport SEXP _metapack_pRho_to_Rho(SEXP pRhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type pRho(pRhoSEXP);
    rcpp_result_gen = Rcpp::wrap(pRho_to_Rho(pRho));
    return rcpp_result_gen;
END_RCPP
}
// calc_modelfit_lpml
Rcpp::List calc_modelfit_lpml(const arma::vec& y, const arma::mat& x, const arma::mat& z, const arma::uvec& ids, const arma::uvec& iarm, const arma::vec& npt, const double& nu, const arma::mat& betas, const arma::mat& sig2s, const arma::mat& phis, const arma::mat& lams, const arma::cube& Rhos, const int& K, const int& nT, const int& nkeep, const bool verbose);
RcppExport SEXP _metapack_calc_modelfit_lpml(SEXP ySEXP, SEXP xSEXP, SEXP zSEXP, SEXP idsSEXP, SEXP iarmSEXP, SEXP nptSEXP, SEXP nuSEXP, SEXP betasSEXP, SEXP sig2sSEXP, SEXP phisSEXP, SEXP lamsSEXP, SEXP RhosSEXP, SEXP KSEXP, SEXP nTSEXP, SEXP nkeepSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ids(idsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type iarm(iarmSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type npt(nptSEXP);
    Rcpp::traits::input_parameter< const double& >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type betas(betasSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sig2s(sig2sSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type phis(phisSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lams(lamsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Rhos(RhosSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int& >::type nT(nTSEXP);
    Rcpp::traits::input_parameter< const int& >::type nkeep(nkeepSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_modelfit_lpml(y, x, z, ids, iarm, npt, nu, betas, sig2s, phis, lams, Rhos, K, nT, nkeep, verbose));
    return rcpp_result_gen;
END_RCPP
}
// rwish
arma::mat rwish(const double& v, const arma::mat& S);
RcppExport SEXP _metapack_rwish(SEXP vSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(rwish(v, S));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_metapack_BMVMR_POCov", (DL_FUNC) &_metapack_BMVMR_POCov, 21},
    {"_metapack_BayesNMR", (DL_FUNC) &_metapack_BayesNMR, 19},
    {"_metapack_calc_modelfit_dic", (DL_FUNC) &_metapack_calc_modelfit_dic, 16},
    {"_metapack_veclinv", (DL_FUNC) &_metapack_veclinv, 2},
    {"_metapack_pRho_to_Rho", (DL_FUNC) &_metapack_pRho_to_Rho, 1},
    {"_metapack_calc_modelfit_lpml", (DL_FUNC) &_metapack_calc_modelfit_lpml, 16},
    {"_metapack_rwish", (DL_FUNC) &_metapack_rwish, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_metapack(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
