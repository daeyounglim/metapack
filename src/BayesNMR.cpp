#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <algorithm>
#include <iterator>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "ListBuilder.h"
#include "linearalgebra.h"
#include "misc_nmr.h"
#include "nelmin.h"
// [[Rcpp::depends(RcppArmadillo,RcppProgress))]]

// [[Rcpp::export]]
Rcpp::List BayesNMR(const arma::vec& y,
					const arma::vec& sd,
					const arma::mat& x,
					const arma::mat& z,
					const arma::uvec& ids,
					const arma::uvec& iarm,
					const arma::vec& npt,
					const double& nu, // = degrees of freedom
					const double& c01_inv,
					const double& c02_inv,
					const int& K,
					const int& nT,
					const int& ndiscard,
					const int& nskip,
					const int& nkeep,
					const bool verbose,
					const arma::vec& beta_init,
					const arma::vec& phi_init,
					const arma::vec& sig2_init) {
	using namespace arma;
	using namespace Rcpp;
	using namespace R;
	using namespace std;

	const int ns = y.n_elem;
	const int nx = x.n_cols;
	const int nz = z.n_cols;
	const vec sd2 = arma::pow(sd, 2.0);

	bool t_random_effect = false;
	if (R_FINITE(nu)) {
		t_random_effect = true;
	}

	/*************************
	Parameters for adaptive MH
	*************************/
	vec lam_rates(K, fill::zeros);
	vec phi_rates(nz, fill::zeros);
	mat Rho_rates(nT,nT, fill::zeros);


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


	/********************
	Initialize parameters 
	********************/
	vec beta = beta_init;
	vec phi = phi_init;
	vec sig2 = sig2_init;
	vec xb(ns, fill::zeros);
	for (int k=0; k < K; ++k) {
		mat E_k = Eks(k);
		mat X_k = arma::join_horiz(Xks(k), E_k.t());
		uvec idx = idxks(k);
		xb(idx) = X_k * beta;
	}
	vec resid = y - xb;
	vec Rgam(ns, fill::zeros);
	mat pRho(nT,nT, fill::zeros);
	mat Rho(nT,nT, fill::zeros);
	for (int i = 0; i < nT; ++i) {
		pRho(i,i) = 1.0;
	}
	Rho = pRho_to_Rho(pRho);


	vec Z = arma::exp(z * phi);

	vec lam(K, fill::ones);

	/*******************
	Begin burn-in period
	*******************/
	if (verbose) {
		Rcout << "Warming up" << endl;
	}
	{
		Progress prog(ndiscard, verbose);
		for (int idiscard = 0; idiscard < ndiscard; ++idiscard) {
			if (Progress::check_abort()) {
				return Rcpp::List::create(Rcpp::Named("error") = "user interrupt aborted");
			}
			/**********
			Sample sig2
			**********/
			vec shape = 0.5 * npt + 0.00001;
			vec rate = 0.5 * npt % arma::pow(resid - Z % Rgam, 2.0) + 0.5 * (npt - 1.0) % sd2 + 0.00001;
			for (int i = 0; i < ns; ++i) {
				sig2(i) = rate(i) / ::Rf_rgamma(shape(i), 1.0);
			}

			/**********
			Sample beta
			**********/
			mat SigBetainv(nx+nT, nx+nT, fill::zeros);
			vec muBetaTmp(nx+nT, fill::zeros);
			for (int k=0; k < K; ++k) {
				uvec idx = idxks(k);
				mat Sig_k = arma::diagmat(sig2(idx) / npt(idx));
				vec y_k = y(idx);
				mat E_k = Eks(k);
				mat X_k = arma::join_horiz(Xks(k), E_k.t());
				mat Z_k = arma::diagmat(Z(idx));

				mat tmpmat = Z_k * E_k.t() * Rho * E_k * Z_k / lam(k) + Sig_k;
				SigBetainv += X_k.t() * arma::solve(tmpmat, X_k);
				muBetaTmp += X_k.t() * arma::solve(tmpmat, y_k);
			}
			SigBetainv.diag() += c01_inv;
			SigBetainv = 0.5 * (SigBetainv + SigBetainv.t());
			mat SigBetainvChol = chol(SigBetainv);
			vec muBeta = arma::solve(arma::trimatu(SigBetainvChol), arma::solve(trimatl(SigBetainvChol.t()), muBetaTmp));
			vec btmp(nx+nT);
			std::generate(btmp.begin(), btmp.end(), ::norm_rand);
			beta = muBeta + arma::solve(arma::trimatu(SigBetainvChol), btmp);
			for (int k=0; k < K; ++k) {
				mat E_k = Eks(k);
				mat X_k = arma::join_horiz(Xks(k), E_k.t());
				uvec idx = idxks(k);
				xb(idx) = X_k * beta;
			}
			resid = y - xb;

			/***********************
			Sample eta = log(lambda)
			***********************/
			if (t_random_effect) {
				for (int k=0; k < K; ++k) {
					uvec idx = idxks(k);
					vec sig2_k = sig2(idx) / npt(idx);
					vec resid_k = resid(idx);
					mat E_k = Eks(k);
					vec z_k = Z(idx);
					double lam_k = lam(k);
					mat ERE_k = E_k.t() * Rho * E_k;

					auto fx_lam = [&](double eta_input[])->double {
						return -loglik_eta(eta_input[0], nu, resid_k, z_k, ERE_k, sig2_k);
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

					double minll = ynewlo;
					double xmax = xmin[0];

					mat cl(5,3, fill::zeros);
					vec dl(5, fill::zeros);
					double step_size = 0.5;
					eta_burnin_block:
						for (int iii=0; iii < 5; ++iii) {
							double e1 = static_cast<double>(iii-2);
							cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
							cl(iii,1) = xmax + e1 * step_size;
							cl(iii,2) = 1.0;
							dl(iii) = -loglik_eta(xmax + e1 * step_size, nu, resid_k, z_k, ERE_k, sig2_k);
						}

					for (int ni=0; ni < 5; ++ni) {
						if ((ni+1) != 3) {
							if (dl(ni) <= minll) {
								step_size *= 1.5;
								goto eta_burnin_block;
							}
						}
					}
					vec fl = solve(cl.t() * cl,  cl.t() * dl);
					double sigmaa = std::sqrt(0.5 / fl(0));



					double eta_k = std::log(lam_k);
					double eta_prop = ::norm_rand() * sigmaa + xmax;
					
					// log-likelihood difference
					double ll_diff = loglik_eta(eta_prop, nu, resid_k, z_k, ERE_k, sig2_k) - 
									loglik_eta(eta_k, nu, resid_k, z_k, ERE_k, sig2_k) -
									0.5 * (std::pow(eta_k - xmax, 2.0) - std::pow(eta_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);

					if (std::log(::unif_rand()) < ll_diff) {
						lam(k) = std::exp(eta_prop);
						++lam_rates(k);
					}
				}
			}


			/*********
			Sample phi
			*********/
			for (int g=0; g < nz; ++g) {
				vec phi_prop = phi;
				auto fx = [&](double phig[])->double {
					phi_prop(g) = phig[0];
					vec Z = arma::exp(z * phi_prop);
					const int K = Eks.n_elem;
					double loglik = -0.5 * c02_inv * arma::accu(phi_prop % phi_prop);
					for (int k=0; k < K; ++k) {
						uvec idx_k = idxks(k);
						mat Z_k = arma::diagmat(Z(idx_k));
						vec sig2_k = sig2(idx_k) / npt(idx_k);
						mat E_k = Eks(k);
						mat ERE = E_k.t() * Rho * E_k;
						vec resid_k = resid(idx_k);


						mat tmpmat = Z_k * ERE * Z_k / lam(k);
						tmpmat.diag() += sig2_k;
						double logdet_val;
						double logdet_sign;
						log_det(logdet_val, logdet_sign, tmpmat);

						loglik += -0.5 * logdet_val - 0.5 * arma::accu(resid_k % arma::solve(tmpmat, resid_k));
					}
					return -loglik;
				};
				double start[] = { phi(g) };
				double xmin[] = { 0.0 };
				double ynewlo = 0.0;
				double reqmin = 1.0e-20;
				int konvge = 5;
				int kcount = 1000;
				double step[] = { 0.05 };
				int icount = 0;
				int numres = 0;
				int ifault = 0;
				nelmin(fx, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);

				double xmax = xmin[0];
				double minll = ynewlo;

				mat cl(5,3, fill::zeros);
				vec dl(5, fill::zeros);
				double step_size = 0.5;
				phi_burnin_block:
					for (int iii=0; iii < 5; ++iii) {
						double e1 = static_cast<double>(iii-2);
						cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
						cl(iii,1) = xmax + e1 * step_size;
						cl(iii,2) = 1.0;
						phi_prop(g) = xmax + e1 * step_size;
						dl(iii) = -loglik_phi(phi_prop, z, c02_inv, lam, sig2 / npt, Rho, resid, Eks, idxks);
					}

				for (int ni=0; ni < 5; ++ni) {
					if ((ni+1) != 3) {
						if (dl(ni) <= minll) {
							step_size *= 1.5;
							goto phi_burnin_block;
						}
					}
				}

				vec fl = solve(cl.t() * cl, cl.t() * dl);
				double sigmaa = std::sqrt(0.5 / fl(0));

				phi_prop(g) = ::norm_rand() * sigmaa + xmax;

				// log-likelihood difference
				double ll_diff = loglik_phi(phi_prop, z, c02_inv, lam, sig2 / npt, Rho, resid, Eks, idxks) -
						  loglik_phi(phi, z, c02_inv, lam, sig2 / npt, Rho, resid, Eks, idxks) -
						  0.5 * (std::pow(phi(g) - xmax, 2.0) - std::pow(phi_prop(g) - xmax, 2.0)) / std::pow(sigmaa, 2.0);
				if (std::log(::unif_rand()) < ll_diff) {
					phi(g) = phi_prop(g);
					++phi_rates(g);
				}
			}
			Z = arma::exp(z * phi);


			/*********
			Sample Rho
			*********/
			for (int iR=0; iR < nT-1; ++iR) {
				for (int iC=iR+1; iC < nT; ++iC) {
					double zprho = 0.5 * std::log((1.0 + pRho(iR,iC)) / (1.0 - pRho(iR,iC)));
					auto fx_rho = [&](double zprho_input[])->double {
						return -loglik_z(zprho_input[0], iR, iC, pRho, lam, sig2 / npt, Z, resid, Eks, idxks);
					};

					double start[] = { zprho };
					double xmin[] = { 0.0 };
					double ynewlo = 0.0;
					double reqmin = 1.0e-20;
					int konvge = 5;
					int kcount = 1000;
					double step[] = { 0.02 };
					int icount = 0;
					int numres = 0;
					int ifault = 0;
					nelmin(fx_rho, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);
					double xmax = xmin[0];
					double minll = ynewlo;

					mat cl(5,3, fill::zeros);
					vec dl(5, fill::zeros);
					double step_size = 0.2;
					Rho_burnin_block:
						for (int iii=0; iii < 5; ++iii) {
							double e1 = static_cast<double>(iii-2);
							cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
							cl(iii,1) = xmax + e1 * step_size;
							cl(iii,2) = 1.0;
							dl(iii) = -loglik_z(xmax + e1 * step_size, iR, iC, pRho, lam, sig2 / npt, Z, resid, Eks, idxks);
						}

					for (int ni=0; ni < 5; ++ni) {
						if ((ni+1) != 3) {
							if (dl(ni) <= minll) {
								step_size *= 1.2;
								goto Rho_burnin_block;
							}
						}
					}

					vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
					double sigmaa = std::sqrt(0.5 / fl(0));


					double zprho_prop = ::norm_rand() * sigmaa + xmax;

					// log-likelihood difference
					double ll_diff = loglik_z(zprho_prop, iR, iC, pRho, lam, sig2 / npt, Z, resid, Eks, idxks) -
									 loglik_z(zprho, iR, iC, pRho, lam, sig2 / npt, Z, resid, Eks, idxks) -
									 0.5 * (std::pow(zprho - xmax, 2.0) - std::pow(zprho_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
					if (std::log(::unif_rand()) < ll_diff) {
						pRho(iR,iC) = (std::exp(2.0 * zprho_prop) - 1.0) / (std::exp(2.0 * zprho_prop) + 1.0);
						pRho(iC,iR) = pRho(iR,iC);
						Rho = pRho_to_Rho(pRho);
						++Rho_rates(iR,iC);
					}
				}
			}

			/***********
			Sample Rgam
			***********/
			for (int k=0; k < K; ++k) {
				uvec idx = idxks(k);
				mat Z_k = arma::diagmat(Z(idx));
				vec sig2_k = sig2(idx) / npt(idx);
				mat E_k = Eks(k);
				mat ERE = E_k.t() * Rho * E_k;
				vec resid_k = resid(idx);

				mat SigRgami = Z_k * diagmat(1.0 / sig2_k) * Z_k + lam(k) * ERE.i();
				SigRgami = 0.5 * (SigRgami + SigRgami.t());
				
				mat SigRgamiChol = chol(SigRgami);
				vec muRgam = arma::solve(arma::trimatu(SigRgamiChol), arma::solve(arma::trimatl(SigRgamiChol.t()), Z_k * diagmat(1.0 / sig2_k) * resid_k));

				int n_k = idx.n_elem;
				vec gtmp(n_k);
				std::generate(gtmp.begin(), gtmp.end(), ::norm_rand);
				Rgam(idx) = muRgam + arma::solve(arma::trimatu(SigRgamiChol), gtmp);
			}

			prog.increment();
		}
	}

	/***********************
	Begin posterior sampling
	***********************/
	mat beta_save(nx+nT, nkeep, arma::fill::zeros);
	mat sig2_save(ns, nkeep, arma::fill::zeros);
	mat lam_save(K, nkeep, arma::fill::zeros);
	mat phi_save(nz, nkeep, arma::fill::zeros);
	cube Rho_save(nT,nT,nkeep, arma::fill::zeros);
	mat gam_save(ns, nkeep, arma::fill::zeros);

	if (verbose) {
		Rcout << "Saving posterior samples" << endl;
	}
	{
		Progress prog(nkeep, verbose);
		for (int ikeep=0; ikeep < nkeep; ++ikeep) {
			// R_CheckUserInterrupt();
			if (Progress::check_abort()) {
				return Rcpp::List::create(Rcpp::Named("error") = "user interrupt aborted");
			}
			for (int iskip=0; iskip < nskip; ++iskip) {
				/**********
				Sample sig2
				**********/
				vec shape = 0.5 * npt + 0.00001;
				vec rate = 0.5 * npt % arma::pow(resid - Z % Rgam, 2.0) + 0.5 * (npt - 1.0) % sd2 + 0.00001;
				for (int i = 0; i < ns; ++i) {
					sig2(i) = rate(i) / ::Rf_rgamma(shape(i), 1.0);
				}

				/**********
				Sample beta
				**********/
				mat SigBetainv(nx+nT, nx+nT, fill::zeros);
				vec muBetaTmp(nx+nT, fill::zeros);
				for (int k=0; k < K; ++k) {
					uvec idx = idxks(k);
					mat Sig_k = arma::diagmat(sig2(idx) / npt(idx));
					vec y_k = y(idx);
					mat E_k = Eks(k);
					mat X_k = arma::join_horiz(Xks(k), E_k.t());
					mat Z_k = arma::diagmat(Z(idx));

					mat tmpmat = Z_k * E_k.t() * Rho * E_k * Z_k / lam(k) + Sig_k;
					SigBetainv += X_k.t() * arma::solve(tmpmat, X_k);
					muBetaTmp += X_k.t() * arma::solve(tmpmat, y_k);
				}
				SigBetainv.diag() += c01_inv;
				SigBetainv = 0.5 * (SigBetainv + SigBetainv.t());
				mat SigBetainvChol = chol(SigBetainv);
				
				vec muBeta = arma::solve(arma::trimatu(SigBetainvChol), arma::solve(trimatl(SigBetainvChol.t()), muBetaTmp));
				vec btmp(nx+nT);
				std::generate(btmp.begin(), btmp.end(), ::norm_rand);
				beta = muBeta + arma::solve(arma::trimatu(SigBetainvChol), btmp);
				for (int k=0; k < K; ++k) {
					mat E_k = Eks(k);
					mat X_k = arma::join_horiz(Xks(k), E_k.t());
					uvec idx = idxks(k);
					xb(idx) = X_k * beta;
				}
				resid = y - xb;

				/***********************
				Sample eta = log(lambda)
				***********************/
				if (t_random_effect) {
					for (int k=0; k < K; ++k) {
						uvec idx = idxks(k);
						vec sig2_k = sig2(idx) / npt(idx);
						vec resid_k = resid(idx);
						mat E_k = Eks(k);
						vec z_k = Z(idx);
						double lam_k = lam(k);
						mat ERE_k = E_k.t() * Rho * E_k;

			            auto fx_lam = [&](double eta_input[])->double {
			              return -loglik_eta(eta_input[0], nu, resid_k, z_k, ERE_k, sig2_k);
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
						double minll = ynewlo;
						double xmax = xmin[0];

						mat cl(5,3, fill::zeros);
						vec dl(5, fill::zeros);
						double step_size = 0.5;
						eta_sample_block:
							for (int iii=0; iii < 5; ++iii) {
								double e1 = static_cast<double>(iii-2);
								cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
								cl(iii,1) = xmax + e1 * step_size;
								cl(iii,2) = 1.0;
								dl(iii) = -loglik_eta(xmax + e1 * step_size, nu, resid_k, z_k, ERE_k, sig2_k);
							}

						for (int ni=0; ni < 5; ++ni) {
							if ((ni+1) != 3) {
								if (dl(ni) <= minll) {
									step_size *= 1.5;
									goto eta_sample_block;
								}
							}
						}

						// mat clinv = arma::solve(cl.t() * cl);
						vec fl = solve(cl.t() * cl, cl.t() * dl);
						double sigmaa = std::sqrt(0.5 / fl(0));



						double eta_k = std::log(lam(k));
						double eta_prop = ::norm_rand() * sigmaa + xmax;
						
						// log-likelihood difference
						double ll_diff = loglik_eta(eta_prop, nu, resid_k, z_k, ERE_k, sig2_k) - 
										loglik_eta(eta_k, nu, resid_k, z_k, ERE_k, sig2_k) -
										0.5 * (std::pow(eta_k - xmax, 2.0) - std::pow(eta_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);

						if (std::log(::unif_rand()) < ll_diff) {
							lam(k) = std::exp(eta_prop);
							++lam_rates(k);
						}
					}
				}


				/*********
				Sample phi
				*********/
				for (int g=0; g < nz; ++g) {
					vec phi_prop = phi;
					auto fx = [&](double phig[])->double {
						phi_prop(g) = phig[0];
						vec Z = arma::exp(z * phi_prop);
						const int K = Eks.n_elem;
						double loglik = -0.5 * c02_inv * arma::accu(phi_prop % phi_prop);
						for (int k=0; k < K; ++k) {
							uvec idx_k = idxks(k);
							mat Z_k = arma::diagmat(Z(idx_k));
							vec sig2_k = sig2(idx_k) / npt(idx_k);
							mat E_k = Eks(k);
							mat ERE = E_k.t() * Rho * E_k;
							vec resid_k = resid(idx_k);

							mat tmpmat = Z_k * ERE * Z_k / lam(k);
							tmpmat.diag() += sig2_k;
							double logdet_val;
							double logdet_sign;
							log_det(logdet_val, logdet_sign, tmpmat);

							loglik += -0.5 * logdet_val - 0.5 * arma::accu(resid_k % arma::solve(tmpmat, resid_k));
						}
						return -loglik;
					};
					double start[] = { phi(g) };
					double xmin[] = { 0.0 };
					double ynewlo = 0.0;
					double reqmin = 1.0e-20;
					int konvge = 5;
					int kcount = 1000;
					double step[] = { 0.05 };
					int icount = 0;
					int numres = 0;
					int ifault = 0;
					nelmin(fx, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);
					double xmax = xmin[0];
					double minll = ynewlo;

					mat cl(5,3, fill::zeros);
					vec dl(5, fill::zeros);
					double step_size = 0.5;
					phi_sample_block:
						for (int iii=0; iii < 5; ++iii) {
							double e1 = static_cast<double>(iii-2);
							cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
							cl(iii,1) = xmax + e1 * step_size;
							cl(iii,2) = 1.0;
							phi_prop(g) = xmax + e1 * step_size;
							dl(iii) = -loglik_phi(phi_prop, z, c02_inv, lam, sig2 / npt, Rho, resid, Eks, idxks);
						}

					for (int ni=0; ni < 5; ++ni) {
						if ((ni+1) != 3) {
							if (dl(ni) <= minll) {
								step_size *= 1.5;
								goto phi_sample_block;
							}
						}
					}

					vec fl = solve(cl.t() * cl, cl.t() * dl);
					double sigmaa = std::sqrt(0.5 / fl(0));

					phi_prop(g) = ::norm_rand() * sigmaa + xmax;

					// log-likelihood difference
					double ll_diff = loglik_phi(phi_prop, z, c02_inv, lam, sig2 / npt, Rho, resid, Eks, idxks) -
							  loglik_phi(phi, z, c02_inv, lam, sig2 / npt, Rho, resid, Eks, idxks) -
							  0.5 * (std::pow(phi(g) - xmax, 2.0) - std::pow(phi_prop(g) - xmax, 2.0)) / std::pow(sigmaa, 2.0);
					if (std::log(::unif_rand()) < ll_diff) {
						phi(g) = phi_prop(g);
						++phi_rates(g);
					}
				}
				Z = arma::exp(z * phi);

				/*********
				Sample Rho
				*********/
				for (int iR=0; iR < nT-1; ++iR) {
					for (int iC=iR+1; iC < nT; ++iC) {
						double zprho = 0.5 * std::log((1.0 + pRho(iR,iC)) / (1.0 - pRho(iR,iC)));
			            auto fx_rho = [&](double zprho_input[])->double {
			              return -loglik_z(zprho_input[0], iR, iC, pRho, lam, sig2 / npt, Z, resid, Eks, idxks);
			            };

			            double start[] = { zprho };
			            double xmin[] = { 0.0 };
			            double ynewlo = 0.0;
			            double reqmin = 1.0e-20;
			            int konvge = 5;
			            int kcount = 1000;
			            double step[] = { 0.02 };
			            int icount = 0;
			            int numres = 0;
			            int ifault = 0;
			            nelmin(fx_rho, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);
			            double xmax = xmin[0];
			            double minll = ynewlo;

						mat cl(5,3, fill::zeros);
						vec dl(5, fill::zeros);
						double step_size = 0.5;
						Rho_sample_block:
							for (int iii=0; iii < 5; ++iii) {
								double e1 = static_cast<double>(iii-2);
								cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
								cl(iii,1) = xmax + e1 * step_size;
								cl(iii,2) = 1.0;
								dl(iii) = -loglik_z(xmax + e1 * step_size, iR, iC, pRho, lam, sig2 / npt, Z, resid, Eks, idxks);
							}

						for (int ni=0; ni < 5; ++ni) {
							if ((ni+1) != 3) {
								if (dl(ni) <= minll) {
									step_size *= 1.2;
									goto Rho_sample_block;
								}
							}
						}
						vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
						double sigmaa = std::sqrt(0.5 / fl(0));


						double zprho_prop = ::norm_rand() * sigmaa + xmax;

						// log-likelihood difference
						double ll_diff = loglik_z(zprho_prop, iR, iC, pRho, lam, sig2 / npt, Z, resid, Eks, idxks) -
										 loglik_z(zprho, iR, iC, pRho, lam, sig2 / npt, Z, resid, Eks, idxks) -
										 0.5 * (std::pow(zprho - xmax, 2.0) - std::pow(zprho_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
						if (std::log(::unif_rand()) < ll_diff) {
							pRho(iR,iC) = (std::exp(2.0 * zprho_prop) - 1.0) / (std::exp(2.0 * zprho_prop) + 1.0);
							pRho(iC,iR) = pRho(iR,iC);
							Rho = pRho_to_Rho(pRho);
							++Rho_rates(iR,iC);
						}
					}
				}

				/***********
				Sample Rgam
				***********/
				for (int k=0; k < K; ++k) {
					uvec idx = idxks(k);
					mat Z_k = arma::diagmat(Z(idx));
					vec sig2_k = sig2(idx) / npt(idx);
					mat E_k = Eks(k);
					mat ERE = E_k.t() * Rho * E_k;
					vec resid_k = resid(idx);

					mat SigRgami = Z_k * diagmat(1.0 / sig2_k) * Z_k + lam(k) * ERE.i();
					SigRgami = 0.5 * (SigRgami + SigRgami.t());
					mat SigRgamiChol = chol(SigRgami);
					vec muRgam = arma::solve(arma::trimatu(SigRgamiChol), arma::solve(arma::trimatl(SigRgamiChol.t()), Z_k * diagmat(1.0 / sig2_k) * resid_k));

					int n_k = idx.n_elem;
					vec gtmp(n_k);
					std::generate(gtmp.begin(), gtmp.end(), ::norm_rand);
					Rgam(idx) = muRgam + arma::solve(arma::trimatu(SigRgamiChol), gtmp);
				}
			}
			beta_save.col(ikeep) = beta;
			phi_save.col(ikeep) = phi;
			lam_save.col(ikeep) = lam;
			sig2_save.col(ikeep) = sig2;
			Rho_save.slice(ikeep) = Rho;
			gam_save.col(ikeep) = Rgam;
			

			prog.increment();
		}
	}

	Rho_rates /= static_cast<double>(ndiscard+nkeep*nskip);
	phi_rates /= static_cast<double>(ndiscard+nkeep*nskip);
	lam_rates /= static_cast<double>(ndiscard+nkeep*nskip);

	return ListBuilder()
	.add("beta", beta_save)
	.add("phi", phi_save)
	.add("lam", lam_save)
	.add("sig2", sig2_save)
	.add("Rho", Rho_save)
	.add("gam", gam_save)
	.add("Rho_acceptance", Rho_rates)
	.add("lam_acceptance", lam_rates)
	.add("phi_acceptance", phi_rates);
}
