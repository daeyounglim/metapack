#include <cmath>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <Rdefines.h>
#include "misc_mvmeta.h"
// [[Rcpp::depends(RcppArmadillo)]]


double lmgamma(const double& a,
			   const int& p) {
	// log multivariate gamma
	double out = 0.5 * static_cast<double>(p * (p - 1)) * M_LN_SQRT_PI;
	for (int j = 0; j < p; ++j) {
		out += R::lgammafn(a + 0.5 * static_cast<double>(1 - j));
	}
	return out;
}


double loglik_tau(const double& star,
				  const arma::mat& zz,
				  const arma::field<arma::uvec>& idxks,
				  const int& sum_nk,
				  const int& J,
				  const double& a1,
				  const double& b1) {
	using namespace arma;
	const int K = idxks.n_elem;
	const double tau = 1.0 + std::exp(star);

	// compute the gamma function
	double suma = 0.0;
	for (int j = 1; j < J+1; ++j) {
		suma += std::log(tau + static_cast<double>(J - j + 1));
	}

	double loglik = static_cast<double>(sum_nk) * (suma + (tau + 1.0) * std::log(tau))
				  + (a1 - 1.0) * std::log(tau) - b1 * tau + star;
	for (int k = 0; k < K; ++k) {
		uvec idxk = idxks(k);
		mat zz_k = zz.rows(idxk);
		loglik -= (tau + static_cast<double>(J+1)) * arma::accu(arma::log(tau + arma::sum(zz_k, 1)));
	}
	return loglik;
}

double loglik_vv(const double& star,
				 const arma::mat& resid,
				 const arma::vec& lambda,
				 const arma::field<arma::uvec>& idxks,
				 const arma::mat& Sigma,
				 const double& nu0,
				 const double& a2,
				 const double& b2,
				 const int& J) {
	using namespace arma;

	const int K = idxks.n_elem;
	const double vv = nu0 + std::exp(star);
	mat Sig_star = (vv - static_cast<double>(J+1)) * Sigma;
	double logdet_val;
	double logdet_sign;
	log_det(logdet_val, logdet_sign, Sig_star);
	double loglik = 0.5 * vv * logdet_val + (a2 - 1.0) * std::log(vv) - b2 * vv + star;
	for (int k = 0; k < K; ++k) {
		uvec idx_k = idxks(k);
		mat resid_k = resid.rows(idx_k);
		vec lambda_k = lambda(idx_k);
		mat tmpmat = Sig_star + resid_k.t() * arma::diagmat(lambda_k) * resid_k;
		log_det(logdet_val, logdet_sign, tmpmat);		
		int n_k = idx_k.n_elem;
		loglik += lmgamma(0.5 * (vv + static_cast<double>(n_k)), J) - lmgamma(0.5 * vv, J)
				  - 0.5 * (vv + static_cast<double>(n_k)) * logdet_val;
	}
	return loglik;
}


