#ifndef BAYESMETA_RTMVN_H
#define BAYESMETA_RTMVN_H
#include <stdio.h>
#include <cmath>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <Rdefines.h>
#include "random.h"

using namespace Rcpp;
using namespace R;
using namespace arma;

arma::mat rtmvn_gibbs(const int& n, const int& p, const arma::vec& Mean, const arma::mat& Sigma_chol,
                      const arma::mat& R, const arma::vec& a, const arma::vec& b, arma::vec& z);

arma::mat rtmvn(const int& n, const arma::vec& Mean, const arma::mat& Sigma, const arma::mat& D,
                   const arma::vec& lower, const arma::vec& upper, const arma::vec& init);

#endif
