#include <cmath>
#include <Rmath.h>
#include <algorithm>
#include <iterator>
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <Rdefines.h>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include "nelmin.h"
#include "misc_nmr.h"
// [[Rcpp::depends(RcppArmadillo, RcppProgress, BH)]]
using namespace arma;

/**************************************************
Calculate the goodness of fit measures

+ Dev(theta) = -2 * log L(theta | D_oy)
+ p_D = E[Dev(theta)] - Dev(thetabar)
+ DIC = Dev(thetabar) + 2 * p_D
***************************************************/

// [[Rcpp::export]]
Rcpp::List calc_modelfit_dic(const arma::vec& y,
						 const arma::mat& x,
						 const arma::mat& z,
						 const arma::uvec& ids,
						 const arma::uvec& iarm,
						 const arma::vec& npt,
						 const double& nu,
						 const arma::mat& betas,
						 const arma::mat& sig2s,
						 const arma::mat& phis,
						 const arma::mat& lams,
						 const arma::cube& Rhos,
						 const int& K,
						 const int& nT,
						 const int& nkeep,
						 const bool verbose) {
	using namespace arma;
	using namespace boost::math::quadrature;

	/* make a list of y_k, X_k, z_k*/
	arma::field<arma::mat> Xks(K);
	arma::field<arma::mat> Eks(K);
	arma::field<arma::uvec> idxks(K);
	for (int k = 0; k < K; ++k) {
		uvec idx = find(ids == k+1);
		idxks(k) = idx;
		Xks(k) = x.rows(idx);
		int idx_l = idx.n_elem;
		mat Ek(nT, idx_l, fill::zeros);
		
		uvec iarm_k = iarm(idx);
		for (int j = 0; j < idx_l; ++j) {
			Ek(iarm_k(j),j) = 1.0;
		}
		Eks(k) = Ek;
	}

	bool t_random_effect = false;
	if (R_FINITE(nu)) {
		t_random_effect = true;
	}

	vec beta_est = arma::mean(betas, 1);
	vec sig2_est = arma::mean(sig2s, 1);
	vec lam_est = arma::mean(lams, 1);
	vec phi_est = arma::mean(phis, 1);
	mat Rho_est(nT, nT, fill::zeros);
	for (int ikeep=0; ikeep < nkeep; ++ikeep) {
		Rho_est += Rhos.slice(ikeep);
	}
	Rho_est /= static_cast<double>(nkeep);

	/*******************************************
	Dev(thetabar) = -2 * log L(thetabar | D_oy)
	*******************************************/
	double Dev_thetabar = 0.0;
	vec maxll_est(K, fill::zeros);
	vec Z_est = arma::exp(z * phi_est);
	for (int k=0; k < K; ++k) {
		uvec idx = idxks(k);
		vec y_k = y(idx);
		mat E_k = Eks(k);
		mat X_k = arma::join_horiz(Xks(k), E_k.t());
		vec resid_k = y_k - X_k * beta_est;
		vec Z_k = Z_est(idx);
		mat ERE_k = diagmat(Z_k) * E_k.t() * Rho_est * E_k * diagmat(Z_k);
		vec sig2_k = sig2_est(idx) / npt(idx);

		int Tk = idx.n_elem;
		
    
    	if (t_random_effect) {
			auto fx_lam = [&](double eta[])->double {
				return -loglik_lam(eta[0], nu, resid_k, ERE_k, sig2_k, Tk);
			};
			double start[] = { std::log(lam_est(k)) };
			double xmin[] = { 0.0 };
			double ynewlo = 0.0;
			double reqmin = 1.0e-20;
			int konvge = 5;
			int kcount = 1000;
			double step[] = { 0.2 };
			int icount = 0;
			int numres = 0;
			int ifault = 0;
			nelmin(fx_lam, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);
			double maxll = -ynewlo;
			if (R_IsNaN(maxll)) {
				if (k != 0) {
					maxll = maxll_est(k-1);
				} else {
					maxll = 0.0;
				}
			}
			maxll_est(k) = maxll;
			

			auto fx = [&](double lam)->double {
				double loglik = (0.5 * nu - 1.0) * std::log(lam) - 0.5 * nu * lam + 0.5 * nu * (std::log(nu) - M_LN2) - R::lgammafn(0.5 * nu);
				mat ZEREZ_S = diagmat(Z_k) * E_k.t() * Rho_est * E_k * diagmat(Z_k / lam);
				ZEREZ_S.diag() += sig2_k;
				double logdet_val;
				double logdet_sign;
				log_det(logdet_val, logdet_sign, ZEREZ_S);
				loglik -= 0.5 * (logdet_val + arma::accu(resid_k % arma::solve(ZEREZ_S, resid_k))) + M_LN_SQRT_2PI * static_cast<double>(Tk);
				/***********************************
				subtract by maximum likelihood value
				for numerical stability
				***********************************/
				return std::exp(loglik - maxll);
			};

		    double error;
			double Q = gauss_kronrod<double, 15>::integrate(fx, 0.0, std::numeric_limits<double>::infinity(), 5, 1e-10, &error);
			Dev_thetabar += -2.0 * (maxll + std::log(Q));
    	} else {
    		double loglik = -M_LN_SQRT_2PI * static_cast<double>(Tk);
    		// mat ZEREZ_S = diagmat(Z_k) * E_k.t() * Rho_est * E_k * diagmat(Z_k);
			ERE_k.diag() += sig2_k;
    		double logdet_val;
			double logdet_sign;
			log_det(logdet_val, logdet_sign, ERE_k);
			loglik -= 0.5 * (logdet_val + arma::accu(resid_k % arma::solve(ERE_k, resid_k)));
    		Dev_thetabar += -2.0 * loglik;
    	}
	}


	double Dev_bar = 0;
	mat maxll_keep(K,nkeep,fill::zeros);
	{
		Progress prog(nkeep, verbose);
		for (int ikeep = 0; ikeep < nkeep; ++ikeep) {
			if (Progress::check_abort()) {
				return Rcpp::List::create(Rcpp::Named("error") = "user interrupt aborted");
			}

			vec beta_ikeep = betas.col(ikeep);
			vec sig2_ikeep = sig2s.col(ikeep);
			vec phi_ikeep = phis.col(ikeep);
			vec lam_ikeep = lams.col(ikeep);
			mat Rho_ikeep = Rhos.slice(ikeep);
			vec Z_ikeep = arma::exp(z * phi_ikeep);

			for (int k=0; k < K; ++k) {
				uvec idx = idxks(k);
				vec y_k = y(idx);
				mat E_k = Eks(k);
				mat X_k = arma::join_horiz(Xks(k), E_k.t());
				vec resid_k = y_k - X_k * beta_ikeep;
				vec Z_k = Z_ikeep(idx);
				double lam_k = lam_ikeep(k);
				mat ERE_k = diagmat(Z_k) * E_k.t() * Rho_ikeep * E_k * diagmat(Z_k);
				vec sig2_k = sig2_ikeep(idx) / npt(idx);


				int Tk = idx.n_elem;

				if (t_random_effect) {
					auto fx_lam = [&](double eta[])->double {
						return -loglik_lam(eta[0], nu, resid_k, ERE_k, sig2_k, Tk);
					};
					double start[] = { std::log(lam_k) };
					double xmin[] = { 0.0 };
					double ynewlo = 0.0;
					double reqmin = 1.0e-20;
					int konvge = 5;
					int kcount = 1000;
					double step[] = { 0.2 };
					int icount = 0;
					int numres = 0;
					int ifault = 0;
					nelmin(fx_lam, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);
					double maxll = -ynewlo;
					if (R_IsNaN(maxll)) {
						if (k != 0) {
							maxll = maxll_keep(k-1,ikeep);
						} else {
							maxll = 0.0;
						}
					}
					maxll_keep(k,ikeep) = maxll;
					

					auto fx = [&](double lam)->double {
						double loglik = (0.5 * nu - 1.0) * std::log(lam) - 0.5 * nu * lam + 0.5 * nu * (std::log(nu) - M_LN2) - R::lgammafn(0.5 * nu);
						mat ZEREZ_S = diagmat(Z_k) * E_k.t() * Rho_ikeep * E_k * diagmat(Z_k / lam);
						ZEREZ_S.diag() += sig2_k;
						double logdet_val;
						double logdet_sign;
						log_det(logdet_val, logdet_sign, ZEREZ_S);
						loglik -= 0.5 * (logdet_val + arma::accu(resid_k % arma::solve(ZEREZ_S, resid_k))) + M_LN_SQRT_2PI * static_cast<double>(Tk);
						/***********************************
						subtract by maximum likelihood value
						for numerical stability
						***********************************/
						return std::exp(loglik - maxll);
					};

					double error;
					double Q = gauss_kronrod<double, 15>::integrate(fx, 0.0, std::numeric_limits<double>::infinity(), 5, 1e-10, &error);
					Dev_bar += -2.0 * (maxll + std::log(Q));
				} else {
					double loglik = -M_LN_SQRT_2PI * static_cast<double>(Tk);
					// mat ZEREZ_S = diagmat(Z_k) * E_k.t() * Rho_ikeep * E_k * diagmat(Z_k);
					ERE_k.diag() += sig2_k;
		    		double logdet_val;
					double logdet_sign;
					log_det(logdet_val, logdet_sign, ERE_k);
					loglik -= 0.5 * (logdet_val + arma::accu(resid_k % arma::solve(ERE_k, resid_k)));
		    		Dev_bar += -2.0 * loglik;
				}
			}

			prog.increment();
		}
		Dev_bar /= static_cast<double>(nkeep);
		double p_D = Dev_bar - Dev_thetabar;
		double DIC = Dev_thetabar + 2.0 * p_D;
		return Rcpp::List::create(Rcpp::Named("dic")=DIC);
	}
}















