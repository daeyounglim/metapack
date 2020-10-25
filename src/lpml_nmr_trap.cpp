#include <cmath>
#include <Rmath.h>
#include <algorithm>
#include <iterator>
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <Rdefines.h>
#include "nelmin.h"
#include "misc_nmr.h"
#ifdef _OPENMP
  #include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]

/**************************************************
Calculate the goodness of fit measures

+ Dev(theta) = -2 * log L(theta | D_oy)
+ p_D = E[Dev(theta)] - Dev(thetabar)
+ DIC = Dev(thetabar) + 2 * p_D
***************************************************/
// [[Rcpp::export]]
Rcpp::List calc_modelfit_lpml_trap(const arma::vec& y,
						 const arma::mat& x,
						 const arma::mat& z,
						 const arma::uvec& ids,
						 const arma::uvec& iarm,
						 const arma::vec& npt,
						 const arma::vec& dfs,
						 const double& nu,
						 const arma::mat& betas,
						 const arma::mat& sig2s,
						 const arma::mat& phis,
						 const arma::mat& lams,
						 const arma::cube& Rhos,
						 const int& K,
						 const int& nT,
						 const int& nkeep,
						 const bool& sample_df,
						 const bool& verbose,
						 const int& ncores) {
	using namespace arma;

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

	mat Qs(K, nkeep, fill::zeros);
	mat g(K, nkeep, fill::zeros);
	vec gmax(K, fill::zeros);
	vec alogcpo(K, fill::zeros);
	double alpml = 0.0;


	const double h = 0.5;
	mat maxll_keep(K,nkeep,fill::zeros);
	{
		Progress prog(nkeep, verbose);
		#ifdef _OPENMP
		#pragma omp parallel for schedule(static) num_threads(ncores)
		#endif
		for (int ikeep = 0; ikeep < nkeep; ++ikeep) {
			if (!Progress::check_abort()) {
				vec beta_ikeep = betas.col(ikeep);
				vec sig2_ikeep = sig2s.col(ikeep);
				vec phi_ikeep = phis.col(ikeep);
				vec lam_ikeep = lams.col(ikeep);
				mat Rho_ikeep = Rhos.slice(ikeep);
				vec Z_ikeep = arma::exp(z * phi_ikeep);
				double df_ikeep = nu;
				if (sample_df) {
					df_ikeep = dfs(ikeep);
				}

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
						auto fx_lam = [&](double eta[]) -> double {
							double loglik = loglik_lam(eta[0], df_ikeep, resid_k, ERE_k, sig2_k, Tk);
							loglik += 0.5 * df_ikeep * (std::log(df_ikeep) - M_LN2) - R::lgammafn(0.5 * df_ikeep) - M_LN_SQRT_2PI * static_cast<double>(Tk);
							return -loglik;
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
						if (maxll > 100) {
							maxll = maxll_keep(k-1, ikeep);
						}
						maxll_keep(k,ikeep) = maxll;
						

						auto fx = [&](double lam)->double {
							double loglik = (0.5 * df_ikeep - 1.0) * std::log(lam) - 0.5 * df_ikeep * lam + 0.5 * df_ikeep * (std::log(df_ikeep) - M_LN2) - R::lgammafn(0.5 * df_ikeep);
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

						double lam_start = 0.0001;
						std::vector<double> ygrid_;
						while (true) {
							double fv = fx(lam_start);
							if (fv > 0.0) {
								ygrid_.push_back(fv);
								lam_start += h;
							} else {
								break;
							}
						}

						
						vec ygrid = arma::conv_to<arma::vec>::from(ygrid_);
						double ngrid = ygrid.n_elem;
						double Q = 0.5 * ygrid(0) + 0.5 * ygrid(ngrid-1);
						for (int ii = 1; ii < ngrid-1; ++ii) {
							Q += ygrid(ii);
						}
						Q *= h;
						Qs(k, ikeep) = Q;
						g(k,ikeep) -= maxll + std::log(Q);
					} else {
						double loglik = -M_LN_SQRT_2PI * static_cast<double>(Tk);
						ERE_k.diag() += sig2_k;
			    		double logdet_val;
						double logdet_sign;
						log_det(logdet_val, logdet_sign, ERE_k);
						loglik -= 0.5 * (logdet_val + arma::accu(resid_k % arma::solve(ERE_k, resid_k)));
			    		g(k,ikeep) -= loglik;
					}
				}
				prog.increment();
			}
		}
		for (int k = 0; k < K; ++k) {
			gmax(k) = g(k,0);
			for (int j1 = 1; j1 < nkeep; ++j1) {
				if (gmax(k) < g(k, j1)) {
					gmax(k) = g(k, j1);
				}
			}
			double sumrep = 0.0;
			for (int j1 = 1; j1 < nkeep; ++j1) {
				sumrep += std::exp(g(k,j1) - gmax(k));
			}
			alogcpo(k) -= gmax(k) + std::log(sumrep / static_cast<double>(nkeep));
			alpml += alogcpo(k);
		}
		// Dev_bar /= static_cast<double>(nkeep);
		// double p_D = Dev_bar - Dev_thetabar;
		// double DIC = Dev_thetabar + 2.0 * p_D;
		return Rcpp::List::create(Rcpp::Named("lpml")=alpml, Rcpp::Named("logcpo")=alogcpo, Rcpp::Named("g")=g, Rcpp::Named("Qs")=Qs);
	}
}




