#ifndef BAYESMETA_MISCNMR_H
#define BAYESMETA_MISCNMR_H
#include <cmath>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <Rdefines.h>
// [[Rcpp::depends(RcppArmadillo)]]

double loglik_lam(const double& lam_k,
				  const double& nu,
				  const arma::vec& resid_k, // = y_k - X_k * beta
				  const arma::mat& ERE_k, // = Z(phi) *  E_k' * Rho * E_k * Z(phi)
				  const arma::vec& sig2_k,
				  const int& Tk);

double loglik_eta(const double& eta_k,
				  const double& nu,
				  const arma::vec& resid_k, // = y_k - X_k * beta
				  const arma::vec& Z_k,
				  const arma::mat& ERE, // = E_k' * Rho * E_k
				  const arma::vec& sig2_k);

double loglik_phi(const arma::vec& phi,
				  const arma::mat& z,
				  const double& c02inv,
				  const arma::vec& lam,
				  const arma::vec& sig2,
				  const arma::mat& Rho,
				  const arma::vec& resid,
				  const arma::field<arma::mat>& Eks,
				  const arma::field<arma::uvec>& idxks);

arma::mat pRho_to_Rho(arma::mat& pRho);

double loglik_z(const double& zprho,
				const int& index1,
				const int& index2,
				const arma::mat& pRho,
				const arma::vec& lam,
				const arma::vec& sig2,
				const arma::vec& Z,
				const arma::vec& resid,
				const arma::field<arma::mat>& Eks,
				const arma::field<arma::uvec>& idxks);


#endif
