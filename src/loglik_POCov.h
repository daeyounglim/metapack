#ifndef BAYESMETA_LOGLIK_POCOV_H
#define BAYESMETA_LOGLIK_POCOV_H
#include <cmath>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <Rdefines.h>
#include "linearalgebra.h"
#include "loglik_POCov.h"
// [[Rcpp::depends(RcppArmadillo)]]

double loglik_rik(const double& rstar,
				  const arma::rowvec& vrtk,
				  const int& kk,
				  const int& J,
				  const int& iR,
				  const int& iC,
				  const double& ntk,
				  const arma::mat& VSV);

double loglik_delta_m3(const double& logdel,
					   const arma::rowvec& delta_i,
					   const int& j,
					   const arma::mat& Rhoinv,
					   const arma::mat& qq,
					   const double& a0,
					   const double& b0,
					   const double& ntk);
double loglik_rho_m3(const double& zprho,
					 const arma::vec& vRho,
					 const arma::mat& qq,
					 const int& ii,
					 const int& iR,
					 const int& iC,
					 const int& J,
					 const double& sumNpt);
double loglik_delta_m4(const double& logdel,
					   const arma::vec& delta,
					   const int& j,
					   const arma::mat& Rho,
					   const arma::mat& vRtk,
					   const arma::mat& gamR,
					   const arma::uvec& Trial,
					   const arma::vec& Npt,
					   const arma::mat& SD,
					   const arma::mat& resid,
					   const arma::mat& WCovariate,
					   const int& N,
					   const int& J,
					   const int& K,
					   const int& T,
					   const double& d0,
					   const double& nu0,
					   const arma::mat& Sigma0inv);
double loglik_rho_m4(const double& zprho,
	        		 const arma::vec& vRho,
	        		 const int& ii,
	        		 const arma::vec& delta,
	        		 const arma::mat& WCovariate,
	        		 const arma::mat& SD,
	        		 const arma::mat& resid,
	        		 const arma::vec& Npt,
	        		 const arma::mat& vRtk,
	        		 const arma::uvec& Trial,
	        		 const arma::mat& gamR,
	        		 const int& iR,
	        		 const int& iC,
	        		 const double& d0,
	        		 const double& nu0,
	        		 const int& N,
	        		 const int& J,
	        		 const int& K,
	        		 const int& T,
	        		 const arma::mat& Sigma0inv);
#endif