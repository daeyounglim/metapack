#ifndef BAYESMETA_LINEARALGEBRA_H
#define BAYESMETA_LINEARALGEBRA_H
#include <stdio.h>
#include <cmath>
#include <RcppArmadillo.h>
#include <Rdefines.h>
#include "linearalgebra.h"
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec vecl(const arma::mat& X);
arma::mat veclinv(const arma::vec& v, const int& n);
arma::vec vech(const arma::mat& X);
arma::mat vechinv(const arma::vec& v, const int& n);
arma::mat duplicate_matrix (const int& n);
arma::vec uppertriv(const arma::mat& A);
arma::mat blockdiag( arma::field<arma::mat>& x );
arma::mat pRho_to_Rho(arma::mat& pRho);

#endif