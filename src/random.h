#ifndef BAYESMETA_RANDOM_H
#define BAYESMETA_RANDOM_H
#include <stdio.h>
#include <cmath>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <Rdefines.h>

namespace RNG
{
	double _gig_mode(const double& lambda, const double& omega);
	double _rgig_ROU_noshift (const double& lambda, const double& lambda_old, const double& omega, const double& alpha);
	double _rgig_newapproach1 (const double& lambda, const double& lambda_old, const double& omega, const double& alpha);
	double _rgig_ROU_shift_alt (const double& lambda, const double& lambda_old, const double& omega, const double& alpha);
	double gigrnd(const double& lambda, const double& chi, const double& psi);
	double besselM3(const double& lambda, const double& x, const bool& logvalue);
	double EGIG_x(const double& lambda, const double& chi, const double& psi);
	double EGIG_xinv(const double& lambda, const double& chi, const double& psi);
	double rtgamma(const double& a, const double& b, const double& truncpoint, const bool& up);
	arma::mat rwish(const double& v, const arma::mat& S);
	arma::mat riwish(const double& v, const arma::mat& S);
	double tnormrnd(const double& mu, const double& sigma, const double& low, const double& up);
	double rtnormrnd(const double& mu, const double& sigma, const double& up);
	double ltnormrnd(const double& mu, const double& sigma, const double& low);
}
#endif
