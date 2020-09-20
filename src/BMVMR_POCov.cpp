#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <algorithm>
#include <iterator>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "misc.h"
#include "random.h"
#include "linearalgebra.h"
#include "loglik_POCov.h"
#include "nelmin.h"
#include "ListBuilder.h"
// [[Rcpp::depends(RcppArmadillo,RcppProgress))]]

// [[Rcpp::export]]
Rcpp::List BMVMR_POCov(const arma::mat& Outcome,
					   const arma::mat& SD,
					   const arma::mat& XCovariate,
					   const arma::mat& WCovariate,
					   const arma::uvec& Treat,
					   const arma::uvec& Trial,
					   const arma::vec& Npt,
					   const double& c0,
					   const double& dj0, // hyperparameter for Omega
					   const double& d0, // hyperparameter for Sigma
					   const double& s0, // hyperparameter for Sigma
					   const double& nu0, // hyperparameter for Sigma
					   const double& a0,
					   const double& b0,
					   const arma::mat& Omega0,
					   const arma::mat& Sigma0,
					   const int& K, // # of Trials
					   const int& T, // # of Treatments
					   const int& fmodel,
					   const int& ndiscard,
					   const int& nskip,
					   const int& nkeep,
					   const double& R_stepsize,
					   const double& Rho_stepsize,
					   const double& delta_stepsize,
					   const int& delta_rep,
					   const int& rho_rep,
					   const bool& verbose) {
	using namespace arma;
	using namespace std;
	using namespace Rcpp;
	using namespace R;

	const int N = Outcome.n_rows;
	const int J = Outcome.n_cols;
	const int xcols = XCovariate.n_cols;
	const int nw = WCovariate.n_cols;
	const int nt = (xcols + nw) * J;

	arma::field<arma::uvec> idxks(K);
	for (int k = 0; k < K; ++k) {
		uvec idx = find(Trial == k);
		idxks(k) = idx;
	}

	/***********************
	Parameter Initialization
	***********************/

	// vec theta(nt, fill::zeros);
	mat gamR(nw*J, K, fill::zeros);
	mat Omegainv(nw*J, nw*J, fill::eye);
	mat Omega(arma::size(Omegainv), fill::eye);
	mat Sig_lt(N, (J * (J + 1)) /2, fill::zeros); // Sigma lower triangular := vech(Sigma_{tk})
	mat Siginv_lt(N, (J * (J + 1)) /2, fill::zeros); // Sigma^{-1} lower triangular := vech(Sigma_{tk}^{-1})
	mat vRtk(N, J * (J - 1) / 2, fill::zeros);
	vRtk.fill(0.1);
	mat delta = SD;
	vec vRho(J*(J-1)/2);
	std::generate(vRho.begin(), vRho.end(), ::norm_rand);
	mat pRho = arma::tanh(vecrinv(vRho, J));
	pRho.diag().fill(1.0);
	mat Rho = pRho_to_Rho(pRho);
	mat Rhoinv = arma::inv_sympd(Rho);
	mat xx(nt, nt, fill::zeros);
	vec xy(nt, fill::zeros);
	for (int i = 0; i < N; ++i) {
		if (fmodel == 3) {
			mat Sigma = arma::diagmat(delta.row(i)) * Rho * arma::diagmat(delta.row(i));
			mat Sigmainv = arma::diagmat(1.0 / delta.row(i)) * Rhoinv * arma::diagmat(1.0 / delta.row(i));
			Sig_lt.row(i) = arma::trans(vech(Sigma));
			Siginv_lt.row(i) = arma::trans(vech(Sigmainv));
		} else if (fmodel == 2) {
			mat Sigma = arma::diagmat(delta.row(0)) * Rho * arma::diagmat(delta.row(0));
			mat Sigmainv = arma::diagmat(1.0 / delta.row(0)) * Rhoinv * arma::diagmat(1.0 / delta.row(0));
			Sig_lt.row(i) = arma::trans(vech(Sigma));
			Siginv_lt.row(i) = arma::trans(vech(Sigmainv));
		}
		vec y_i = arma::trans(Outcome.row(i));
		rowvec x_i = XCovariate.row(i);
		rowvec w_i = WCovariate.row(i);
		mat X(J, xcols*J, fill::zeros);
		mat W(J, nw*J, fill::zeros);
		for (int j = 0; j < J; ++j) {
			X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
			W(j, span(j*nw, (j+1)*nw-1)) = w_i;
		}
		mat Xstar = join_horiz(X, W);
		xx += Xstar.t() * Xstar;
		xy += Xstar.t() * y_i;
	}
	vec theta = arma::solve(xx, xy);
	mat Delta(J, J, fill::eye);

	
	mat Sig0 = Delta * Rho * Delta;

	/************
	Miscellaneous
	************/
	const mat Omega0inv = arma::inv_sympd(Omega0);
	const mat Sigma0inv = arma::inv_sympd(Sigma0);
	const double sumNpt = arma::accu(Npt);
	const double df = d0 + sumNpt;
	const double shape_omega = static_cast<double>(K) + dj0;
	mat resid = Outcome;
	mat vR_rates(arma::size(vRtk), fill::zeros);
	mat delta_rates(arma::size(delta), fill::zeros);
	mat Delta_rates(arma::size(Delta), fill::zeros);
	vec vRho_rates(arma::size(vRho), fill::zeros);

	/*********
	Containers
	*********/
	mat theta_save(nt, nkeep, fill::zeros);
	mat Delta_save(J, nkeep, fill::zeros);
	cube Sig_save(N, (J * (J + 1)) / 2, nkeep, fill::zeros);
	cube Omegainv_save(nw*J, nw*J, nkeep, fill::zeros);
	cube Omega_save(nw*J, nw*J, nkeep, fill::zeros);
	cube Sigma_save(J, J, nkeep, fill::zeros);
	cube delta_save(N, J, nkeep, fill::zeros);
	cube Rho_save(J, J, nkeep, fill::zeros);
	cube Rtk_save(N, J * (J - 1) / 2, nkeep, fill::zeros);
	cube ypred_save(N, J, nkeep, fill::zeros);
	cube pRtk_save(N, J*(J-1)/2, nkeep, fill::zeros);
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
			
			/***********
			Update theta
			***********/
			mat Sig_theta(nt, nt, fill::zeros);
			Sig_theta.diag().fill(1.0 / c0);
			vec mu_theta(nt, fill::zeros);
			for (int k = 0; k < K; ++k) {
				uvec idxk = idxks(k);
				int n_k = idxk.n_elem;
				mat XSX(nt, nt, fill::zeros);
				mat WSX(nw*J, nt, fill::zeros);
				vec WSy(nw*J, fill::zeros);
				vec XSy(nt, fill::zeros);
				mat Sig_gamk = Omegainv;
				for (int i = 0; i < n_k; ++i) {
					int i_k = idxk(i);
					rowvec x_i = XCovariate.row(i_k);
					rowvec w_i = WCovariate.row(i_k);
					rowvec y_i = Outcome.row(i_k);
					mat Sigmainv = vechinv(arma::trans(Siginv_lt.row(i_k)), J);
					double ntk = Npt(i_k);
					mat X(J, xcols * J, fill::zeros);
					mat W(J, nw*J, fill::zeros);
					for (int j = 0; j < J; ++j) {
						X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
						W(j, span(j*nw, (j+1)*nw-1)) = w_i;
					}
					mat Xstar = join_horiz(X, W);
					XSX += ntk * Xstar.t() * Sigmainv * Xstar;
					XSy += ntk * Xstar.t() * Sigmainv * y_i.t();
					mat WS = ntk * W.t() * Sigmainv;
					Sig_gamk += WS * W;
					WSX += WS * Xstar;
					WSy += WS * y_i.t();
				}
				mat Sig_gamk_inv = arma::inv_sympd(Sig_gamk);
				Sig_theta += XSX - WSX.t() * Sig_gamk_inv * WSX;
				mu_theta += XSy - WSX.t() * Sig_gamk_inv * WSy;
			}
			Sig_theta = 0.5 * (Sig_theta + Sig_theta.t());
			mat Sig_thetaChol = chol(Sig_theta);
			vec atmp(nt);
			std::generate(atmp.begin(), atmp.end(), ::norm_rand);
			vec muTheta = arma::solve(arma::trimatu(Sig_thetaChol), arma::solve(trimatl(Sig_thetaChol.t()), mu_theta));
			theta = muTheta + arma::solve(arma::trimatu(Sig_thetaChol), atmp);
			for (int i = 0; i < N; ++i) {
				rowvec x_i = XCovariate.row(i);
				rowvec w_i = WCovariate.row(i);
				rowvec y_i = Outcome.row(i);
				mat X(J, xcols * J, fill::zeros);
				mat W(J, nw*J, fill::zeros);
				for (int j = 0; j < J; ++j) {
					X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
					W(j, span(j*nw, (j+1)*nw-1)) = w_i;
				}
				mat Xstar = join_horiz(X, W);
				vec resid_i = y_i.t() - Xstar * theta;
				resid.row(i) = resid_i.t();
			}


			/**********
			Update gamR
			**********/
			for (int k = 0; k < K; ++k) {
				mat Siggam = Omegainv;
				vec mugam(nw*J, fill::zeros);
				uvec idxk = idxks(k);
				int n_k = idxk.n_elem;
				for (int i = 0; i < n_k; ++i) {
					int i_k = idxk(i);
					rowvec x_i = XCovariate.row(i_k);
					rowvec w_i = WCovariate.row(i_k);
					rowvec resid_i = resid.row(i_k);
					double ntk = Npt(i_k);
					mat Sigmainv = vechinv(arma::trans(Siginv_lt.row(i_k)), J);
					
					mat W(J, nw*J, fill::zeros);
					for (int j = 0; j < J; ++j) {
						W(j, span(j*nw, (j+1)*nw-1)) = w_i;
					}
					mat WS = W.t() * Sigmainv;
					Siggam += ntk * WS * W;
					mugam += ntk * WS * resid_i.t();
				}
				Siggam = 0.5 * (Siggam + Siggam.t());
				mat SiggamChol = chol(Siggam);
				vec gtmp(nw*J);
				std::generate(gtmp.begin(), gtmp.end(), ::norm_rand);
				vec muGam = arma::solve(arma::trimatu(SiggamChol), arma::solve(trimatl(SiggamChol.t()), mugam));
				gamR.col(k) = muGam + arma::solve(arma::trimatu(SiggamChol), gtmp);
			}
			/***********
			Update Omega
			***********/
			for (int jj = 0; jj < J; ++jj) {
				mat gamstar = gamR.rows(jj*nw, (jj+1)*nw-1);
				mat qq = Omega0inv + gamstar * gamstar.t();

				mat ominv = rwish(shape_omega, arma::inv_sympd(qq));
				mat om = arma::inv_sympd(ominv);
				Omegainv(span(jj*nw, (jj+1)*nw-1), span(jj*nw, (jj+1)*nw-1)) = ominv;
				Omega(span(jj*nw, (jj+1)*nw-1), span(jj*nw, (jj+1)*nw-1)) = om;
			}

			/***********
			Update Sigma
			***********/
			if (fmodel == 1) {
				for (int i = 0; i < N; ++i) {
					int k = Trial(i);
					double ntk = Npt(i);
					double shape = a0 + 0.5 * ntk;
					rowvec sd2_i = arma::square(SD.row(i));
					vec gam_k = gamR.col(k);
					rowvec w_i = WCovariate.row(i);
					mat W(J, nw*J, fill::zeros);
					for (int j = 0; j < J; ++j) {
						W(j, span(j*nw, (j+1)*nw-1)) = w_i;
					}
					vec resid_i2 = arma::square(arma::trans(resid.row(i)) - W * gam_k);
					vec sig2inv(J);
					for (int j = 0; j < J; ++j) {
						double rate = b0 + ntk * 0.5 * resid_i2(j) + (ntk - 1.0) * 0.5 * sd2_i(j);
						sig2inv(j) = ::Rf_rgamma(shape, 1.0) / rate;
					}
					Siginv_lt.row(i) = arma::trans(vech(arma::diagmat(sig2inv)));
					Sig_lt.row(i) = arma::trans(vech(arma::diagmat(1.0 / sig2inv)));
				}
			} else {
				if (fmodel == 2) {
					mat qq = Sigma0inv;
					for (int i = 0; i < N; ++i) {
						int k = Trial(i);
						double ntk = Npt(i);
						vec gam_k = gamR.col(k);

						mat pRR = vecrinv(trans(arma::tanh(vRtk.row(i))), J);
						pRR.diag().fill(1.0);
						mat R = pRho_to_Rho(pRR);

						mat V = arma::diagmat(SD.row(i));
						rowvec w_i = WCovariate.row(i);
						mat W(J, nw*J, fill::zeros);
						for (int j = 0; j < J; ++j) {
							W(j, span(j*nw, (j+1)*nw-1)) = w_i;
						}
						vec resid_i = arma::trans(resid.row(i)) - W * gam_k;
						qq += ntk * resid_i * resid_i.t() + (ntk - 1.0) * V * R * V;
					}
					mat Siginv_new = rwish(df, arma::inv_sympd(qq));
					rowvec vechsiginv_new = arma::trans(vech(Siginv_new));
					rowvec vechsig_new = arma::trans(vech(arma::inv_sympd(Siginv_new)));
					for (int i = 0; i < N; ++i) {
						Sig_lt.row(i) = vechsig_new;
						Siginv_lt.row(i) = vechsiginv_new;
					}
				} else if (fmodel == 3) {
					// Update sig2
					for (int i = 0; i < N; ++i) {
						int k = Trial(i);
						double ntk = Npt(i);
						rowvec w_i = WCovariate.row(i);
						mat pRR = arma::tanh(vecrinv(trans(vRtk.row(i)), J));
						pRR.diag().fill(1.0);
						mat R = pRho_to_Rho(pRR);
						// mat pRR = constructR(trans(vRtk.row(i)), J);
						// mat R = pRR.t() * pRR;
						vec gam_k = gamR.col(k);
						mat V = arma::diagmat(SD.row(i));
						mat W(J, nw*J, fill::zeros);
						for (int j = 0; j < J; ++j) {
							W(j, span(j*nw, (j+1)*nw-1)) = w_i;
						}
						vec resid_i = arma::trans(resid.row(i)) - W * gam_k;
						mat qq = ntk * resid_i * resid_i.t() + (ntk - 1.0) * V * R * V;

						rowvec delta_i = delta.row(i);

						for (int j = 0; j < J; ++j) {
							auto fx_delta = [&](double delta_input[])->double {
								return -loglik_delta_m3(delta_input[0], delta_i, j, Rhoinv, qq, a0, b0, ntk);
							};

							double dstar = std::log(delta_i(j));
							double start[] = { dstar };
							double xmin[] = { 0.0 };
							double ynewlo = 0.0;
							double reqmin = 1.0e-20;
							int konvge = 5;
							int kcount = 1000;
							double step[] = { 0.2 };
							int icount = 0;
							int numres = 0;
							int ifault = 0;
							nelmin(fx_delta, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);					
							double xmax = xmin[0];
							double minll = ynewlo;

							mat cl(5,3, fill::zeros);
							vec dl(5, fill::zeros);
							double step_size = delta_stepsize;
							delta_burnin_block:
								for (int iii=0; iii < 5; ++iii) {
									double e1 = static_cast<double>(iii-2);
									cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
									cl(iii,1) = xmax + e1 * step_size;
									cl(iii,2) = 1.0;
									dl(iii) = -loglik_delta_m3(xmax + e1 * step_size, delta_i, j, Rhoinv, qq, a0, b0, ntk);
								}

							for (int ni=0; ni < 5; ++ni) {
								if ((ni+1) != 3) {
									if (dl(ni) <= minll) {
										step_size *= 1.2;
										goto delta_burnin_block;
									}
								}
							}

							vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
							double sigmaa = std::sqrt(0.5 / fl(0));

							// log-likelihood difference
							icount = 0;
							bool cont_flag = true; // continue flag
							while (icount < delta_rep && cont_flag) {
								double dstar_prop = ::norm_rand() * sigmaa + xmax;
								double ll_diff = loglik_delta_m3(dstar_prop, delta_i, j, Rhoinv, qq, a0, b0, ntk) - loglik_delta_m3(dstar, delta_i, j, Rhoinv, qq, a0, b0, ntk)
										    - 0.5 * (std::pow(dstar - xmax, 2.0) - std::pow(dstar_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
								if (std::log(::unif_rand()) < ll_diff) {
									delta_i(j) = std::exp(dstar_prop);
									delta(i, j) = delta_i(j);
									mat siginvm = arma::diagmat(1.0 / delta_i);
									mat Siginv_new = siginvm * Rhoinv * siginvm;
									mat Sig_new = arma::diagmat(delta_i) * Rho * arma::diagmat(delta_i);
									Sig_lt.row(i) = arma::trans(vech(Sig_new));
									Siginv_lt.row(i) = arma::trans(vech(Siginv_new));

									++delta_rates(i,j);
									cont_flag = false;
								} else {
									++icount;
								}
							}
						}
					}

					// Update Rho
					mat qq(J, J, fill::zeros);
					for (int i = 0; i < N; ++i) {
						int k = Trial(i);
						double ntk = Npt(i);
						rowvec w_i = WCovariate.row(i);
						mat pRR = arma::tanh(vecrinv(trans(vRtk.row(i)), J));
						pRR.diag().fill(1.0);
						mat R = pRho_to_Rho(pRR);

						// mat pRR = constructR(trans(vRtk.row(i)), J);
						// mat R = pRR.t() * pRR;
						vec gam_k = gamR.col(k);
						mat V = arma::diagmat(SD.row(i));
						mat W(J, nw*J, fill::zeros);
						for (int j = 0; j < J; ++j) {
							W(j, span(j*nw, (j+1)*nw-1)) = w_i;
						}
						vec resid_i = arma::trans(resid.row(i)) - W * gam_k;
						mat dAd = ntk * resid_i * resid_i.t() + (ntk - 1.0) * V * R * V;
						mat siginvm = arma::diagmat(1.0 / delta.row(i));
						qq += siginvm * dAd * siginvm;
					}
					for (int ii = 0; ii < (J*(J-1))/2; ++ii) {
						int iR = J - 2 - static_cast<int>(std::sqrt(-8.0*static_cast<double>(ii) + 4.0*static_cast<double>(J*(J-1))-7.0)/2.0 - 0.5); // row index
						int iC = ii + iR + 1 - (J*(J-1))/2 + ((J-iR)*((J-iR)-1))/2; // column index
						double zprho = vRho(ii);
						auto fx_zrho = [&](double zprho_input[])->double {
							return -loglik_rho_m3(zprho_input[0], vRho, qq, ii, iR, iC, J, sumNpt);
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
						nelmin(fx_zrho, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);					
						double xmax = xmin[0];
						double minll = ynewlo;

						mat cl(5,3, fill::zeros);
						vec dl(5, fill::zeros);
						double step_size = Rho_stepsize;
						pRho_burnin_block:
							for (int iii=0; iii < 5; ++iii) {
								double e1 = static_cast<double>(iii-2);
								cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
								cl(iii,1) = xmax + e1 * step_size;
								cl(iii,2) = 1.0;
								dl(iii) = -loglik_rho_m3(xmax + e1 * step_size, vRho, qq, ii, iR, iC, J, sumNpt);
							}

						for (int ni=0; ni < 5; ++ni) {
							if ((ni+1) != 3) {
								if (dl(ni) <= minll) {
									step_size *= 1.2;
									goto pRho_burnin_block;
								}
							}
						}

						vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
						double sigmaa = std::sqrt(0.5 / fl(0));

						icount = 0;
						bool cont_flag = true; // continue flag
						while (icount < rho_rep && cont_flag) {
							// log-likelihood difference
							double zprho_prop = ::norm_rand() * sigmaa + xmax;
							double ll_diff = loglik_rho_m3(zprho_prop, vRho, qq, ii, iR, iC, J, sumNpt) - loglik_rho_m3(zprho, vRho, qq, ii, iR, iC, J, sumNpt)
									    		- 0.5 * (std::pow(zprho - xmax, 2.0) - std::pow(zprho_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
							if (std::log(::unif_rand()) < ll_diff) {
								vRho(ii) = zprho_prop;
								mat pRR = arma::tanh(vecrinv(vRho, J));
								pRR.diag().fill(1.0);
								Rho = pRho_to_Rho(pRR);
								Rhoinv = arma::inv_sympd(Rho);

								// Update Sigmainvs
								for (int i = 0; i < N; ++i) {
									mat siginvm = arma::diagmat(1.0 / delta.row(i));
									mat Siginv_new = siginvm * Rhoinv * siginvm;
									mat Sig_new = arma::diagmat(delta.row(i)) * Rho * arma::diagmat(delta.row(i));
									Sig_lt.row(i) = arma::trans(vech(Sig_new));
									Siginv_lt.row(i) = arma::trans(vech(Siginv_new));
								}
								++vRho_rates(ii);
								cont_flag = false;
							} else {
								++icount;
							}
						}
					}
				}  else if (fmodel == 4) {
					// Update Delta
					for (int j = 0; j < J; ++j) {
						mat Delta_prop = Delta;
						auto fx_delta = [&](double delta_input[])->double {
							double logdeltap = delta_input[0];
							double deltap = std::exp(logdeltap);
							Delta_prop(j,j) = deltap;
							mat Sig_prop = Delta_prop * Rho * Delta_prop;
							double loglik = (static_cast<double>(T*K) * nu0 + d0 - static_cast<double>(J)) * logdeltap - 0.5 * arma::dot(Sigma0inv, Sig_prop);
							for (int i = 0; i < N; ++i) {
								rowvec w_i = WCovariate.row(i);
								// mat RR = constructR(trans(vRtk.row(i)), J);
								// mat R = RR.t() * RR;

								mat pRR = arma::tanh(vecrinv(trans(vRtk.row(i)), J));
								pRR.diag().fill(1.0);
								mat R = pRho_to_Rho(pRR);

								int k = Trial(i);
								vec gam_k = gamR.col(k);
								mat V = arma::diagmat(SD.row(i));
								mat W(J, nw*J, fill::zeros);
								for (int j = 0; j < J; ++j) {
									W(j, span(j*nw, (j+1)*nw-1)) = w_i;
								}
								double ntk = Npt(i);
								vec resid_i = arma::trans(resid.row(i)) - W * gam_k;
								mat qq = ntk * resid_i * resid_i.t() + (ntk - 1.0) * V * R * V + (nu0 - static_cast<double>(J) - 1.0) * Sig_prop;
								double logdet_val;
								double logdet_sign;
								log_det(logdet_val, logdet_sign, qq);
								loglik -= 0.5 * (ntk + nu0) * logdet_val;
							}
						    return -loglik;
						};

						double dstar = std::log(Delta(j,j));
						double start[] = { dstar };
						double xmin[] = { 0.0 };
						double ynewlo = 0.0;
						double reqmin = 1.0e-20;
						int konvge = 5;
						int kcount = 1000;
						double step[] = { 0.2 };
						int icount = 0;
						int numres = 0;
						int ifault = 0;
						nelmin(fx_delta, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);					
						double xmax = xmin[0];
						double minll = ynewlo;

						mat cl(5,3, fill::zeros);
						vec dl(5, fill::zeros);
						double step_size = delta_stepsize;
						delta2_burnin_block:
							for (int iii=0; iii < 5; ++iii) {
								double e1 = static_cast<double>(iii-2);
								cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
								cl(iii,1) = xmax + e1 * step_size;
								cl(iii,2) = 1.0;
								start[0] = xmax + e1 * step_size;
								dl(iii) = fx_delta(start);
							}

						for (int ni=0; ni < 5; ++ni) {
							if ((ni+1) != 3) {
								if (dl(ni) <= minll) {
									step_size *= 1.2;
									goto delta2_burnin_block;
								}
							}
						}

						vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
						double sigmaa = std::sqrt(0.5 / fl(0));


						double dstar_prop = ::norm_rand() * sigmaa + xmax;
						// log-likelihood difference
						start[0] = dstar;
						double ll_diff = fx_delta(start);
						start[0] = dstar_prop;
						ll_diff += -fx_delta(start)
								    - 0.5 * (std::pow(dstar - xmax, 2.0) - std::pow(dstar_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
						if (std::log(::unif_rand()) < ll_diff) {
							Delta(j, j) = std::exp(dstar_prop);
							++Delta_rates(j,j);
						}
					}

					// Update Rho
					for (int ii = 0; ii < (J*(J-1))/2; ++ii) {
						int iR = J - 2 - static_cast<int>(std::sqrt(-8.0*static_cast<double>(ii) + 4.0*static_cast<double>(J*(J-1))-7.0)/2.0 - 0.5); // row index
						int iC = ii + iR + 1 - (J*(J-1))/2 + ((J-iR)*((J-iR)-1))/2; // column index
						double zprho = vRho(ii);
						auto fx_zrho = [&](double zprho_input[])->double {
							double zprhop = zprho_input[0];
							double z = std::tanh(zprhop);
							vec vRhop = vRho;
							vRhop(ii) = zprhop;
							// mat RRp = constructR(vRhop, J);
							// mat Rhop = RRp.t() * RRp;
							// mat pRRp = arma::tanh(vecrinv(vRhop, J));
							// pRRp.diag().fill(1.0);
							// mat Rhop = pRho_to_Rho(pRRp);

							mat pRRp = constructR(vRhop, J);
							mat Rhop = pRRp.t() * pRRp;
							mat RhopChol = arma::chol(Rhop);
							double logdet_val = 2.0 * arma::accu(arma::log(RhopChol.diag()));
							// double logdet_val = 2.0 * arma::prod(RRp.diag());

							mat Sig_prop = Delta * Rhop * Delta;
							double loglik = 0.5 * (static_cast<double>(T*K) * nu0 + d0 - static_cast<double>(J) - 1.0) * logdet_val - 0.5 * arma::dot(Sigma0inv, Sig_prop);
							for (int i = 0; i < N; ++i) {
								rowvec w_i = WCovariate.row(i);

								// mat pRR = arma::tanh(vecrinv(trans(vRtk.row(i)), J));
								// pRR.diag().fill(1.0);
								// mat R = pRho_to_Rho(pRR);
								mat pRR = constructR(trans(vRtk.row(i)), J);
								mat R = pRR.t() * pRR;

								// mat RR = constructR(trans(vRtk.row(i)), J);
								// mat R = RR.t() * RR;
								int k = Trial(i);
								vec gam_k = gamR.col(k);
								mat V = arma::diagmat(SD.row(i));
								mat W(J, nw*J, fill::zeros);
								for (int j = 0; j < J; ++j) {
									W(j, span(j*nw, (j+1)*nw-1)) = w_i;
								}
								double ntk = Npt(i);
								vec resid_i = arma::trans(resid.row(i)) - W * gam_k;
								mat qq = ntk * resid_i * resid_i.t() + (ntk - 1.0) * V * R * V + (nu0 - static_cast<double>(J) - 1.0) * Sig_prop;
								double logdet_sign;
								log_det(logdet_val, logdet_sign, qq);
								loglik -= 0.5 * (ntk + nu0) * logdet_val;
							}
							loglik += 0.5 * (static_cast<double>(J+1-std::abs(iC-iR))) * std::log1p(1.0 - z*z);
						    return -loglik;
						};
						double start[] = { zprho };
						double xmin[] = { 0.0 };
						double ynewlo = 0.0;
						double reqmin = 1.0e-20;
						int konvge = 5;
						int kcount = 1000;
						double step[] = { 0.2 };
						int icount = 0;
						int numres = 0;
						int ifault = 0;
						nelmin(fx_zrho, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);					
						double xmax = xmin[0];
						double minll = ynewlo;

						mat cl(5,3, fill::zeros);
						vec dl(5, fill::zeros);
						double step_size = Rho_stepsize;
						pRho2_burnin_block:
							for (int iii=0; iii < 5; ++iii) {
								double e1 = static_cast<double>(iii-2);
								cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
								cl(iii,1) = xmax + e1 * step_size;
								cl(iii,2) = 1.0;
								start[0] = xmax + e1 * step_size;
								dl(iii) = fx_zrho(start);
							}

						for (int ni=0; ni < 5; ++ni) {
							if ((ni+1) != 3) {
								if (dl(ni) <= minll) {
									step_size *= 1.2;
									goto pRho2_burnin_block;
								}
							}
						}

						vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
						double sigmaa = std::sqrt(0.5 / fl(0));


						double zprho_prop = ::norm_rand() * sigmaa + xmax;
						// log-likelihood difference
						start[0] = zprho;
						double ll_diff = fx_zrho(start);
						start[0] = zprho_prop;
						ll_diff += -fx_zrho(start)
								    - 0.5 * (std::pow(zprho - xmax, 2.0) - std::pow(zprho_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
						if (std::log(::unif_rand()) < ll_diff) {
							vRho(ii) = zprho_prop;
							++vRho_rates(ii);
						}
					}
				}

				// Update R
				for (int i = 0; i < N; ++i) {
					rowvec y_i = Outcome.row(i);
					rowvec vrtk = vRtk.row(i);
					// rowvec phi_tk = R_angles.row(i);
					mat V = arma::diagmat(SD.row(i));
					double ntk = Npt(i);
					mat VSV = V * vechinv(arma::trans(Siginv_lt.row(i)), J) * V;

					for (int kk = 0; kk < (J*(J-1))/2; ++kk) {
						int iR = J - 2 - static_cast<int>(std::sqrt(-8.0*static_cast<double>(kk) + 4.0*static_cast<double>(J*(J-1))-7.0)/2.0 - 0.5); // row index
						int iC = kk + iR + 1 - (J*(J-1))/2 + ((J-iR)*((J-iR)-1))/2; // column index
						auto fx_rstar = [&](double rstar_input[])->double {
							return -loglik_rik(rstar_input[0], vrtk, kk, J, iR, iC, ntk, VSV);
						};
						// double zstar = std::log(phi_tk(kk)) - std::log(M_PI - phi_tk(kk));
						double zstar = vrtk(kk);
						double start[] = { zstar };
						double xmin[] = { 0.0 };
						double ynewlo = 0.0;
						double reqmin = 1.0e-20;
						int konvge = 5;
						int kcount = 1000;
						double step[] = { 0.02 };
						int icount = 0;
						int numres = 0;
						int ifault = 0;
						nelmin(fx_rstar, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);					
						double xmax = xmin[0];
						double minll = ynewlo;

						mat cl(5,3, fill::zeros);
						vec dl(5, fill::zeros);
						double step_size = R_stepsize;

						bool cont_flag = true;
						while (cont_flag) {
							for (int iii=0; iii < 5; ++iii) {
								double e1 = static_cast<double>(iii-2);
								cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
								cl(iii,1) = xmax + e1 * step_size;
								cl(iii,2) = 1.0;
								dl(iii) = -loglik_rik(xmax + e1 * step_size, vrtk, kk, J, iR, iC, ntk, VSV);
							}
							if (any(dl < minll)) {
								step_size *= 1.2;
							} else {
								cont_flag = false;
							}
						}

						vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
						double sigmaa = std::sqrt(0.5 / fl(0));

						double zstar_prop = ::norm_rand() * sigmaa + xmax;
						// log-likelihood difference
						double ll_diff = loglik_rik(zstar_prop, vrtk, kk, J, iR, iC, ntk, VSV) - loglik_rik(zstar, vrtk, kk, J, iR, iC, ntk, VSV)
											- 0.5 * (std::pow(zstar - xmax, 2.0) - std::pow(zstar_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
								    
						if (std::log(::unif_rand()) < ll_diff) {
							vrtk(kk) = zstar_prop;
							vRtk(i,kk) = zstar_prop;
							++vR_rates(i,kk);
							// phi_tk(kk) = std::exp(2.0 * M_LN_SQRT_PI - log1exp(zstar_prop));
							// R_angles.row(i) = phi_tk;
						}
					}
				}
			}
			prog.increment();
		}
	}

	mat ypred(arma::size(Outcome), fill::zeros);
	if (verbose) {
		Rcout << "Sampling" << endl;
	}
	{
		Progress prog(nkeep, verbose);
		int incr = 0;
		for (int ikeep = 0; ikeep < nkeep; ++ikeep) {
			if (Progress::check_abort()) {
				return Rcpp::List::create(Rcpp::Named("error") = "user interrupt aborted");
			}
			for (int iskip = 0; iskip < nskip; ++iskip) {
				/***********
				Update theta
				***********/
				mat Sig_theta(nt, nt, fill::zeros);
				Sig_theta.diag().fill(1.0 / c0);
				vec mu_theta(nt, fill::zeros);
				for (int k = 0; k < K; ++k) {
					uvec idxk = idxks(k);
					int n_k = idxk.n_elem;
					mat XSX(nt, nt, fill::zeros);
					mat WSX(nw*J, nt, fill::zeros);
					vec WSy(nw*J, fill::zeros);
					vec XSy(nt, fill::zeros);
					mat Sig_gamk = Omegainv;
					for (int i = 0; i < n_k; ++i) {
						int i_k = idxk(i);
						rowvec x_i = XCovariate.row(i_k);
						rowvec w_i = WCovariate.row(i_k);
						rowvec y_i = Outcome.row(i_k);
						mat Sigmainv = vechinv(arma::trans(Siginv_lt.row(i_k)), J);
						double ntk = Npt(i_k);
						mat X(J, xcols * J, fill::zeros);
						mat W(J, nw*J, fill::zeros);
						for (int j = 0; j < J; ++j) {
							X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
							W(j, span(j*nw, (j+1)*nw-1)) = w_i;
						}
						mat Xstar = join_horiz(X, W);
						XSX += ntk * Xstar.t() * Sigmainv * Xstar;
						XSy += ntk * Xstar.t() * Sigmainv * y_i.t();
						mat WS = ntk * W.t() * Sigmainv;
						Sig_gamk += WS * W;
						WSX += WS * Xstar;
						WSy += WS * y_i.t();
					}
					mat Sig_gamk_inv = arma::inv_sympd(Sig_gamk);
					Sig_theta += XSX - WSX.t() * Sig_gamk_inv * WSX;
					mu_theta += XSy - WSX.t() * Sig_gamk_inv * WSy;
				}
				Sig_theta = 0.5 * (Sig_theta + Sig_theta.t());
				mat Sig_thetaChol = chol(Sig_theta);
				vec atmp(nt);
				std::generate(atmp.begin(), atmp.end(), ::norm_rand);
				vec muTheta = arma::solve(arma::trimatu(Sig_thetaChol), arma::solve(trimatl(Sig_thetaChol.t()), mu_theta));
				theta = muTheta + arma::solve(arma::trimatu(Sig_thetaChol), atmp);
				for (int i = 0; i < N; ++i) {
					rowvec x_i = XCovariate.row(i);
					rowvec w_i = WCovariate.row(i);
					rowvec y_i = Outcome.row(i);
					mat X(J, xcols * J, fill::zeros);
					mat W(J, nw*J, fill::zeros);
					for (int j = 0; j < J; ++j) {
						X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
						W(j, span(j*nw, (j+1)*nw-1)) = w_i;
					}
					vec resid_i = y_i.t() - join_horiz(X, W) * theta;
					resid.row(i) = resid_i.t();
				}

				/***********
				Update Omega
				***********/
				for (int jj = 0; jj < J; ++jj) {
					mat gamstar = gamR.rows(jj*nw, (jj+1)*nw-1);
					mat qq = Omega0inv + gamstar * gamstar.t();
					mat ominv = rwish(shape_omega, arma::inv_sympd(qq));
					mat om = ominv.i();
					Omegainv(span(jj*nw, (jj+1)*nw-1), span(jj*nw, (jj+1)*nw-1)) = ominv;
					Omega(span(jj*nw, (jj+1)*nw-1), span(jj*nw, (jj+1)*nw-1)) = om;
				}

				/***********
				Update Sigma
				***********/
				if (fmodel == 1) {
					for (int i = 0; i < N; ++i) {
						int k = Trial(i);
						double ntk = Npt(i);
						double shape = a0 + 0.5 * ntk;
						rowvec sd2_i = arma::square(SD.row(i));
						vec gam_k = gamR.col(k);
						rowvec w_i = WCovariate.row(i);
						mat W(J, nw*J, fill::zeros);
						for (int j = 0; j < J; ++j) {
							W(j, span(j*nw, (j+1)*nw-1)) = w_i;
						}
						vec resid_i2 = arma::square(arma::trans(resid.row(i)) - W * gam_k);
						vec sig2inv(J);
						for (int j = 0; j < J; ++j) {
							double rate = b0 + ntk * 0.5 * resid_i2(j) + (ntk - 1.0) * 0.5 * sd2_i(j);
							sig2inv(j) = ::Rf_rgamma(shape, 1.0) / rate;
						}
						Siginv_lt.row(i) = arma::trans(vech(arma::diagmat(sig2inv)));
						Sig_lt.row(i) = arma::trans(vech(arma::diagmat(1.0 / sig2inv)));
					}
				} else {
					if (fmodel == 2) {
						mat qq = Sigma0inv;
						for (int i = 0; i < N; ++i) {
							int k = Trial(i);
							double ntk = Npt(i);
							vec gam_k = gamR.col(k);

							// rowvec phi_tk = R_angles.row(i);
							// mat phimat = veclinv(phi_tk, J);
							// mat Bp = constructB(phimat);
							// mat R = Bp * Bp.t();

							// mat RR = constructR(trans(vRtk.row(i)), J);
							// mat R = RR.t() * RR;

							mat pRR = arma::tanh(vecrinv(trans(vRtk.row(i)), J));
							pRR.diag().fill(1.0);
							mat R = pRho_to_Rho(pRR);
							mat V = arma::diagmat(SD.row(i));
							rowvec w_i = WCovariate.row(i);
							mat W(J, nw*J, fill::zeros);
							for (int j = 0; j < J; ++j) {
								W(j, span(j*nw, (j+1)*nw-1)) = w_i;
							}
							vec resid_i = arma::trans(resid.row(i)) - W * gam_k;
							qq += ntk * resid_i * resid_i.t() + (ntk - 1.0) * V * R * V;
						}
						mat Siginv_new = rwish(df, arma::inv_sympd(qq));
						rowvec vechsiginv_new = arma::trans(vech(Siginv_new));
						rowvec vechsig_new = arma::trans(vech(arma::inv_sympd(Siginv_new)));
						for (int i = 0; i < N; ++i) {
							Sig_lt.row(i) = vechsig_new;
							Siginv_lt.row(i) = vechsiginv_new;
						}
					} else if (fmodel == 3) {
						// Update sig2
						for (int i = 0; i < N; ++i) {
							int k = Trial(i);
							double ntk = Npt(i);
							rowvec w_i = WCovariate.row(i);
							mat pRR = arma::tanh(vecrinv(trans(vRtk.row(i)), J));
							pRR.diag().fill(1.0);
							mat R = pRho_to_Rho(pRR);
							// mat pRR = constructR(trans(vRtk.row(i)), J);
							// mat R = pRR.t() * pRR;
							vec gam_k = gamR.col(k);
							mat V = arma::diagmat(SD.row(i));
							mat W(J, nw*J, fill::zeros);
							for (int j = 0; j < J; ++j) {
								W(j, span(j*nw, (j+1)*nw-1)) = w_i;
							}
							vec resid_i = arma::trans(resid.row(i)) - W * gam_k;
							mat qq = ntk * resid_i * resid_i.t() + (ntk - 1.0) * V * R * V;

							rowvec delta_i = delta.row(i);

							for (int j = 0; j < J; ++j) {
								auto fx_delta = [&](double delta_input[])->double {
									return -loglik_delta_m3(delta_input[0], delta_i, j, Rhoinv, qq, a0, b0, ntk);
								};

								double dstar = std::log(delta_i(j));
								double start[] = { dstar };
								double xmin[] = { 0.0 };
								double ynewlo = 0.0;
								double reqmin = 1.0e-20;
								int konvge = 5;
								int kcount = 1000;
								double step[] = { 0.2 };
								int icount = 0;
								int numres = 0;
								int ifault = 0;
								nelmin(fx_delta, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);					
								double xmax = xmin[0];
								double minll = ynewlo;

								mat cl(5,3, fill::zeros);
								vec dl(5, fill::zeros);
								double step_size = delta_stepsize;
								delta_sampling_block:
									for (int iii=0; iii < 5; ++iii) {
										double e1 = static_cast<double>(iii-2);
										cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
										cl(iii,1) = xmax + e1 * step_size;
										cl(iii,2) = 1.0;
										dl(iii) = -loglik_delta_m3(xmax + e1 * step_size, delta_i, j, Rhoinv, qq, a0, b0, ntk);
									}

								for (int ni=0; ni < 5; ++ni) {
									if ((ni+1) != 3) {
										if (dl(ni) <= minll) {
											step_size *= 1.2;
											goto delta_sampling_block;
										}
									}
								}

								vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
								double sigmaa = std::sqrt(0.5 / fl(0));

								// log-likelihood difference
								icount = 0;
								bool cont_flag = true; // continue flag
								while (icount < delta_rep && cont_flag) {
									double dstar_prop = ::norm_rand() * sigmaa + xmax;
									double ll_diff = loglik_delta_m3(dstar_prop, delta_i, j, Rhoinv, qq, a0, b0, ntk) - loglik_delta_m3(dstar, delta_i, j, Rhoinv, qq, a0, b0, ntk)
											    - 0.5 * (std::pow(dstar - xmax, 2.0) - std::pow(dstar_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
									if (std::log(::unif_rand()) < ll_diff) {
										delta_i(j) = std::exp(dstar_prop);
										delta(i, j) = delta_i(j);
										mat siginvm = arma::diagmat(1.0 / delta_i);
										mat Siginv_new = siginvm * Rhoinv * siginvm;
										mat Sig_new = arma::diagmat(delta_i) * Rho * arma::diagmat(delta_i);
										Sig_lt.row(i) = arma::trans(vech(Sig_new));
										Siginv_lt.row(i) = arma::trans(vech(Siginv_new));

										++delta_rates(i,j);
										cont_flag = false;
									} else {
										++icount;
									}
								}
							}
						}

						// Update Rho
						mat qq(J, J, fill::zeros);
						for (int i = 0; i < N; ++i) {
							int k = Trial(i);
							double ntk = Npt(i);
							rowvec w_i = WCovariate.row(i);
							mat pRR = arma::tanh(vecrinv(trans(vRtk.row(i)), J));
							pRR.diag().fill(1.0);
							mat R = pRho_to_Rho(pRR);

							// mat pRR = constructR(trans(vRtk.row(i)), J);
							// mat R = pRR.t() * pRR;
							vec gam_k = gamR.col(k);
							mat V = arma::diagmat(SD.row(i));
							mat W(J, nw*J, fill::zeros);
							for (int j = 0; j < J; ++j) {
								W(j, span(j*nw, (j+1)*nw-1)) = w_i;
							}
							vec resid_i = arma::trans(resid.row(i)) - W * gam_k;
							mat dAd = ntk * resid_i * resid_i.t() + (ntk - 1.0) * V * R * V;
							mat siginvm = arma::diagmat(1.0 / delta.row(i));
							qq += siginvm * dAd * siginvm;
						}
						for (int ii = 0; ii < (J*(J-1))/2; ++ii) {
							int iR = J - 2 - static_cast<int>(std::sqrt(-8.0*static_cast<double>(ii) + 4.0*static_cast<double>(J*(J-1))-7.0)/2.0 - 0.5); // row index
							int iC = ii + iR + 1 - (J*(J-1))/2 + ((J-iR)*((J-iR)-1))/2; // column index
							auto fx_zrho = [&](double zprho_input[])->double {
								return -loglik_rho_m3(zprho_input[0], vRho, qq, ii, iR, iC, J, sumNpt);
							};
							double zprho = vRho(ii);
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
							nelmin(fx_zrho, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);					
							double xmax = xmin[0];
							double minll = ynewlo;

							mat cl(5,3, fill::zeros);
							vec dl(5, fill::zeros);
							double step_size = Rho_stepsize;
							pRho_sampling_block:
								for (int iii=0; iii < 5; ++iii) {
									double e1 = static_cast<double>(iii-2);
									cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
									cl(iii,1) = xmax + e1 * step_size;
									cl(iii,2) = 1.0;
									dl(iii) = -loglik_rho_m3(xmax + e1 * step_size, vRho, qq, ii, iR, iC, J, sumNpt);
								}

							for (int ni=0; ni < 5; ++ni) {
								if ((ni+1) != 3) {
									if (dl(ni) <= minll) {
										step_size *= 1.2;
										goto pRho_sampling_block;
									}
								}
							}

							vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
							double sigmaa = std::sqrt(0.5 / fl(0));

							icount = 0;
							bool cont_flag = true; // continue flag
							while (icount < rho_rep && cont_flag) {
								// log-likelihood difference
								double zprho_prop = ::norm_rand() * sigmaa + xmax;
								double ll_diff = loglik_rho_m3(zprho_prop, vRho, qq, ii, iR, iC, J, sumNpt) - loglik_rho_m3(zprho, vRho, qq, ii, iR, iC, J, sumNpt)
										    		- 0.5 * (std::pow(zprho - xmax, 2.0) - std::pow(zprho_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
								if (std::log(::unif_rand()) < ll_diff) {
									vRho(ii) = zprho_prop;
									mat pRR = arma::tanh(vecrinv(vRho, J));
									pRR.diag().fill(1.0);
									Rho = pRho_to_Rho(pRR);
									Rhoinv = arma::inv_sympd(Rho);

									// Update Sigmainvs
									for (int i = 0; i < N; ++i) {
										mat siginvm = arma::diagmat(1.0 / delta.row(i));
										mat Siginv_new = siginvm * Rhoinv * siginvm;
										mat Sig_new = arma::diagmat(delta.row(i)) * Rho * arma::diagmat(delta.row(i));
										Sig_lt.row(i) = arma::trans(vech(Sig_new));
										Siginv_lt.row(i) = arma::trans(vech(Siginv_new));
									}
									++vRho_rates(ii);
									cont_flag = false;
								} else {
									++icount;
								}
							}
						}
					}  else if (fmodel == 4) {
						// Update Delta
						for (int j = 0; j < J; ++j) {
							mat Delta_prop = Delta;
							auto fx_delta = [&](double delta_input[])->double {
								double logdeltap = delta_input[0];
								double deltap = std::exp(logdeltap);
								Delta_prop(j,j) = deltap;
								mat Sig_prop = Delta_prop * Rho * Delta_prop;
								double loglik = (static_cast<double>(T*K) * nu0 + d0 - static_cast<double>(J)) * logdeltap - 0.5 * arma::dot(Sigma0inv, Sig_prop);
								for (int i = 0; i < N; ++i) {
									rowvec w_i = WCovariate.row(i);
									mat RR = constructR(trans(vRtk.row(i)), J);
									mat R = RR.t() * RR;

									// mat pRR = arma::tanh(vecrinv(trans(vRtk.row(i)), J));
									// pRR.diag().fill(1.0);
									// mat R = pRho_to_Rho(pRR);

									int k = Trial(i);
									vec gam_k = gamR.col(k);
									mat V = arma::diagmat(SD.row(i));
									mat W(J, nw*J, fill::zeros);
									for (int j = 0; j < J; ++j) {
										W(j, span(j*nw, (j+1)*nw-1)) = w_i;
									}
									double ntk = Npt(i);
									vec resid_i = arma::trans(resid.row(i)) - W * gam_k;
									mat qq = ntk * resid_i * resid_i.t() + (ntk - 1.0) * V * R * V + (nu0 - static_cast<double>(J) - 1.0) * Sig_prop;
									double logdet_val;
									double logdet_sign;
									log_det(logdet_val, logdet_sign, qq);
									loglik -= 0.5 * (ntk + nu0) * logdet_val;
								}
							    return -loglik;
							};

							double dstar = std::log(Delta(j,j));
							double start[] = { dstar };
							double xmin[] = { 0.0 };
							double ynewlo = 0.0;
							double reqmin = 1.0e-20;
							int konvge = 5;
							int kcount = 1000;
							double step[] = { 0.2 };
							int icount = 0;
							int numres = 0;
							int ifault = 0;
							nelmin(fx_delta, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);					
							double xmax = xmin[0];
							double minll = ynewlo;

							mat cl(5,3, fill::zeros);
							vec dl(5, fill::zeros);
							double step_size = delta_stepsize;
							delta2_sampling_block:
								for (int iii=0; iii < 5; ++iii) {
									double e1 = static_cast<double>(iii-2);
									cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
									cl(iii,1) = xmax + e1 * step_size;
									cl(iii,2) = 1.0;
									start[0] = xmax + e1 * step_size;
									dl(iii) = fx_delta(start);
								}

							for (int ni=0; ni < 5; ++ni) {
								if ((ni+1) != 3) {
									if (dl(ni) <= minll) {
										step_size *= 1.2;
										goto delta2_sampling_block;
									}
								}
							}

							vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
							double sigmaa = std::sqrt(0.5 / fl(0));


							double dstar_prop = ::norm_rand() * sigmaa + xmax;
							// log-likelihood difference
							start[0] = dstar;
							double ll_diff = fx_delta(start);
							start[0] = dstar_prop;
							ll_diff += -fx_delta(start)
									    - 0.5 * (std::pow(dstar - xmax, 2.0) - std::pow(dstar_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
							if (std::log(::unif_rand()) < ll_diff) {
								Delta(j, j) = std::exp(dstar_prop);
								++Delta_rates(j,j);
							}
						}

						// Update Rho
						for (int ii = 0; ii < (J*(J-1))/2; ++ii) {
							int iR = J - 2 - static_cast<int>(std::sqrt(-8.0*static_cast<double>(ii) + 4.0*static_cast<double>(J*(J-1))-7.0)/2.0 - 0.5); // row index
							int iC = ii + iR + 1 - (J*(J-1))/2 + ((J-iR)*((J-iR)-1))/2; // column index
							double zprho = vRho(ii);
							auto fx_zrho = [&](double zprho_input[])->double {
								double zprhop = zprho_input[0];
								double z = std::tanh(zprhop);
								vec vRhop = vRho;
								vRhop(ii) = zprhop;
								mat RRp = constructR(vRhop, J);
								mat Rhop = RRp.t() * RRp;
								// mat pRRp = arma::tanh(vecrinv(vRhop, J));
								// pRRp.diag().fill(1.0);
								// mat Rhop = pRho_to_Rho(pRRp);
								// mat RhopChol = arma::chol(Rhop);
								// double logdet_val = 2.0 * arma::accu(arma::log(RhopChol.diag()));
								double logdet_val = 2.0 * arma::accu(arma::log(RRp.diag()));

								mat Sig_prop = Delta * Rhop * Delta;
								double loglik = 0.5 * (static_cast<double>(T*K) * nu0 + d0 - static_cast<double>(J) - 1.0) * logdet_val - 0.5 * arma::dot(Sigma0inv, Sig_prop);
								for (int i = 0; i < N; ++i) {
									rowvec w_i = WCovariate.row(i);

									// mat pRR = arma::tanh(vecrinv(trans(vRtk.row(i)), J));
									// pRR.diag().fill(1.0);
									// mat R = pRho_to_Rho(pRR);

									mat RR = constructR(trans(vRtk.row(i)), J);
									mat R = RR.t() * RR;
									int k = Trial(i);
									vec gam_k = gamR.col(k);
									mat V = arma::diagmat(SD.row(i));
									mat W(J, nw*J, fill::zeros);
									for (int j = 0; j < J; ++j) {
										W(j, span(j*nw, (j+1)*nw-1)) = w_i;
									}
									double ntk = Npt(i);
									vec resid_i = arma::trans(resid.row(i)) - W * gam_k;
									mat qq = ntk * resid_i * resid_i.t() + (ntk - 1.0) * V * R * V + (nu0 - static_cast<double>(J) - 1.0) * Sig_prop;
									double logdet_sign;
									log_det(logdet_val, logdet_sign, qq);
									loglik -= 0.5 * (ntk + nu0) * logdet_val;
								}
								loglik += 0.5 * (static_cast<double>(J+1-std::abs(iC-iR))) * std::log1p(1.0 - z*z);
							    return -loglik;
							};
							double start[] = { zprho };
							double xmin[] = { 0.0 };
							double ynewlo = 0.0;
							double reqmin = 1.0e-20;
							int konvge = 5;
							int kcount = 1000;
							double step[] = { 0.2 };
							int icount = 0;
							int numres = 0;
							int ifault = 0;
							nelmin(fx_zrho, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);					
							double xmax = xmin[0];
							double minll = ynewlo;

							mat cl(5,3, fill::zeros);
							vec dl(5, fill::zeros);
							double step_size = Rho_stepsize;
							pRho2_sampling_block:
								for (int iii=0; iii < 5; ++iii) {
									double e1 = static_cast<double>(iii-2);
									cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
									cl(iii,1) = xmax + e1 * step_size;
									cl(iii,2) = 1.0;
									start[0] = xmax + e1 * step_size;
									dl(iii) = fx_zrho(start);
								}

							for (int ni=0; ni < 5; ++ni) {
								if ((ni+1) != 3) {
									if (dl(ni) <= minll) {
										step_size *= 1.2;
										goto pRho2_sampling_block;
									}
								}
							}

							vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
							double sigmaa = std::sqrt(0.5 / fl(0));


							double zprho_prop = ::norm_rand() * sigmaa + xmax;
							// log-likelihood difference
							start[0] = zprho;
							double ll_diff = fx_zrho(start);
							start[0] = zprho_prop;
							ll_diff += -fx_zrho(start)
									    - 0.5 * (std::pow(zprho - xmax, 2.0) - std::pow(zprho_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
							if (std::log(::unif_rand()) < ll_diff) {
								vRho(ii) = zprho_prop;
								++vRho_rates(ii);
							}
						}
					}

					// Update R
					for (int i = 0; i < N; ++i) {
						rowvec y_i = Outcome.row(i);
						rowvec vrtk = vRtk.row(i);
						// rowvec phi_tk = R_angles.row(i);
						mat V = arma::diagmat(SD.row(i));
						double ntk = Npt(i);
						mat VSV = V * vechinv(arma::trans(Siginv_lt.row(i)), J) * V;

						for (int kk = 0; kk < (J*(J-1))/2; ++kk) {
							int iR = J - 2 - static_cast<int>(std::sqrt(-8.0*static_cast<double>(kk) + 4.0*static_cast<double>(J*(J-1))-7.0)/2.0 - 0.5); // row index
							int iC = kk + iR + 1 - (J*(J-1))/2 + ((J-iR)*((J-iR)-1))/2; // column index
							auto fx_rstar = [&](double rstar_input[])->double {
								return -loglik_rik(rstar_input[0], vrtk, kk, J, iR, iC, ntk, VSV);
							};
							// double zstar = std::log(phi_tk(kk)) - std::log(M_PI - phi_tk(kk));
							double zstar = vrtk(kk);
							double start[] = { zstar };
							double xmin[] = { 0.0 };
							double ynewlo = 0.0;
							double reqmin = 1.0e-20;
							int konvge = 5;
							int kcount = 1000;
							double step[] = { 0.02 };
							int icount = 0;
							int numres = 0;
							int ifault = 0;
							nelmin(fx_rstar, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);					
							double xmax = xmin[0];
							double minll = ynewlo;

							mat cl(5,3, fill::zeros);
							vec dl(5, fill::zeros);
							double step_size = R_stepsize;
							R_sampling_block:
								for (int iii=0; iii < 5; ++iii) {
									double e1 = static_cast<double>(iii-2);
									cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
									cl(iii,1) = xmax + e1 * step_size;
									cl(iii,2) = 1.0;
									dl(iii) = -loglik_rik(xmax + e1 * step_size, vrtk, kk, J, iR, iC, ntk, VSV);
								}

							for (int ni=0; ni < 5; ++ni) {
								if ((ni+1) != 3) {
									if (dl(ni) <= minll) {
										step_size *= 1.2;
										goto R_sampling_block;
									}
								}
							}

							vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
							double sigmaa = std::sqrt(0.5 / fl(0));

							double zstar_prop = ::norm_rand() * sigmaa + xmax;
							// log-likelihood difference
							double ll_diff = loglik_rik(zstar_prop, vrtk, kk, J, iR, iC, ntk, VSV) - loglik_rik(zstar, vrtk, kk, J, iR, iC, ntk, VSV)
												- 0.5 * (std::pow(zstar - xmax, 2.0) - std::pow(zstar_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
									    
							if (std::log(::unif_rand()) < ll_diff) {
								vrtk(kk) = zstar_prop;
								vRtk(i,kk) = zstar_prop;
								++vR_rates(i,kk);
								// phi_tk(kk) = std::exp(2.0 * M_LN_SQRT_PI - log1exp(zstar_prop));
								// R_angles.row(i) = phi_tk;
							}
						}
					}
				}
				++incr;
			}
			theta_save.col(ikeep) = theta;
			Omega_save.slice(ikeep) = Omega;
			ypred_save.slice(ikeep) = ypred;


			if (fmodel == 1) {
				Sig_save.slice(ikeep) = Sig_lt;
			} else {
				if (fmodel == 2) {
					Sigma_save.slice(ikeep) = vechinv(arma::trans(Sig_lt.row(0)), J);
				} else if (fmodel == 3) {
					delta_save.slice(ikeep) = delta;
					Rho_save.slice(ikeep) = Rho;
					Sig_save.slice(ikeep) = Sig_lt;
				} else if (fmodel == 4) {
					Sigma_save.slice(ikeep) = Sig0;
					Sig_save.slice(ikeep) = Sig_lt;
					vec del = arma::sqrt(Sig0.diag());
					Delta_save.col(ikeep) = del;
					Rho_save.slice(ikeep) = arma::diagmat(1.0 / del) * Sig0 * arma::diagmat(1.0 / del);
				}

				mat Rtk(arma::size(vRtk), fill::zeros);
				for (int i = 0; i < N; ++i) {
					mat RR = arma::tanh(vecrinv(trans(vRtk.row(i)), J));
					RR.diag().fill(1.0);
					mat R = pRho_to_Rho(RR);
					// mat RR = constructR(trans(vRtk.row(i)), J);
					// mat R = RR.t() * RR;
					Rtk.row(i) = arma::trans(vecl(R));
					// rowvec phi_tk = R_angles.row(i);
					// mat phimat = veclinv(phi_tk, J);
					// mat Bp = constructB(phimat);
					// mat R = Bp * Bp.t();
					// Rtk.row(i) = arma::trans(vecl(R));
				}
				Rtk_save.slice(ikeep) = Rtk;
				pRtk_save.slice(ikeep) = arma::tanh(vRtk);
			}
			prog.increment();
		}
	}

	ListBuilder out;
	if (fmodel == 1) {
		out = ListBuilder()
			.add("ypred", ypred_save)
			.add("theta", theta_save)
			.add("Omega", Omega_save)
			.add("Sigma", Sig_save);
	} else {
		if (fmodel == 2) {
			out = ListBuilder()
				.add("ypred", ypred_save)
				.add("theta", theta_save)
				.add("Omega", Omega_save)
			    .add("Sigma", Sigma_save)
			    .add("R", Rtk_save)
			    .add("pR", pRtk_save)
			    .add("vR_acceptance", vR_rates / static_cast<double>(ndiscard + nskip*nkeep));
		} else if (fmodel == 3) {
			out = ListBuilder()
				.add("ypred", ypred_save)
				.add("theta", theta_save)
				.add("Omega", Omega_save)
				.add("Sigma", Sig_save)
				.add("delta", delta_save)
				.add("Rho", Rho_save)
				.add("R", Rtk_save)
			    .add("pR", pRtk_save)
				.add("delta_acceptance", delta_rates / static_cast<double>(ndiscard + nskip*nkeep))
				.add("vRho_acceptance", vRho_rates / static_cast<double>(ndiscard + nskip*nkeep))
			    .add("vR_acceptance", vR_rates / static_cast<double>(ndiscard + nskip*nkeep));
		} else if (fmodel == 4) {
			out = ListBuilder()
				.add("ypred", ypred_save)
				.add("theta", theta_save)
				.add("Omega", Omega_save)
				.add("Delta", Delta_save)
				.add("Rho", Rho_save)
				.add("Sigma0", Sigma_save)
				.add("Sigma", Sig_save)
				.add("R", Rtk_save)
			    .add("pR", pRtk_save)
			    .add("vR_acceptance", vR_rates / static_cast<double>(ndiscard + nskip*nkeep));
		} else {
			out = ListBuilder().add("warning", "Returned without a result. Please check your fmodel.");
		}
	}
	return out;
}