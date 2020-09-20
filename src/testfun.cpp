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
#include "nelmin.h"
#include "ListBuilder.h"
// [[Rcpp::depends(RcppArmadillo,RcppProgress))]]

// [[Rcpp::export]]
Rcpp::List testfun(const arma::mat& VSV,
				   const arma::vec& vrtk,
				   const int& j,
				   const int& J,
				   const int& iR,
				   const int& iC,
				   const double& ntk) {
	using namespace arma;
	using namespace std;
	using namespace Rcpp;
	using namespace R;

	auto fx_rstar = [&](double rho_input[])->double {
		double rstar = rho_input[0];
		double z = std::tanh(rstar);
		vec tmp_vrtk = vrtk;
		tmp_vrtk(j-1) = rstar;

		mat pRR_prop = arma::tanh(vecrinv(tmp_vrtk, J));
		pRR_prop.diag().fill(1.0);
		mat R_prop = pRho_to_Rho(pRR_prop);
		mat R_prop_chol = arma::chol(R_prop);
		double logdet_val = 2.0 * arma::accu(arma::log(R_prop_chol.diag()));
		// Rcpp::Rcout << "tanh(tmp_vrtk) = " << arma::tanh(tmp_vrtk) << std::endl;
		// Rcpp::Rcout << "logdet_val = " << logdet_val << std::endl;
		// Rcpp::Rcout << "ntk = " << ntk << std::endl;
		// Rcpp::Rcout << "dot(R_prop, VSV) = " << arma::dot(R_prop, VSV) << std::endl;
		// Rcpp::Rcout << "(J+1-iC+iR)/2*log1p(-z*z) = " << 0.5 * static_cast<double>(J + 1 - iC + iR) * std::log1p(-z*z);
		double loglik = 0.5 * (ntk - static_cast<double>(J) - 2.0) * logdet_val - 0.5 * (ntk - 1.0) * arma::dot(R_prop, VSV);
		loglik += 0.5 * static_cast<double>(J + 1 - iC + iR) * std::log1p(-z*z);
	    return -loglik;
	};

	// double xeval[] = { a };
	// double fval = fx_rstar(xeval);
	// return Rcpp::List::create(Rcpp::Named("fval")=fval);

	double start[] = { vrtk(j-1) };
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
	start[0] = vrtk(j-1);
	double fval = fx_rstar(start);
	return Rcpp::List::create(Rcpp::Named("xmax")=xmax,
					   Rcpp::Named("ymin")=minll,
					   Rcpp::Named("fval")=fval);
}
// Rcpp::List testfun(const double& a) {
// 	using namespace arma;
// 	using namespace std;
// 	using namespace Rcpp;
// 	using namespace R;
// 	auto fx_rstar = [&](double rho_input[])->double {
// 		double r = rho_input[0];
// 		return 3.0 * r*r - 2.0 * r + 5.0;
// 	};

// 	double start[] = { a };
// 	double xmin[] = { 0.0 };
// 	double ynewlo = 0.0;
// 	double reqmin = 1.0e-20;
// 	int konvge = 5;
// 	int kcount = 1000;
// 	double step[] = { 0.02 };
// 	int icount = 0;
// 	int numres = 0;
// 	int ifault = 0;
// 	nelmin(fx_rstar, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);					
// 	double xmax = xmin[0];
// 	double minll = ynewlo;
// 	start[0] = a;
// 	double fval = fx_rstar(start);
// 	return Rcpp::List::create(Rcpp::Named("xmax")=xmax,
// 					   Rcpp::Named("ymin")=minll,
// 					   Rcpp::Named("fval")=fval);
// }



// Rcpp::List testfun(const arma::mat& X,
// 				   const int& ndiscard, const int& nskip, const int& nkeep) {
// 	using namespace arma;

// 	const int n = X.n_cols;
// 	const int p = X.n_rows;

// 	vec delta(p, fill::ones);
// 	vec pRho = vecr(arma::eye(p,p));
// 	mat Rho = arma::eye(p,p);
// 	mat S = X * X.t();

// 	vec pR_rates(p*(p-1)/2, fill::zeros);
// 	vec delta_rates(p, fill::zeros);
// 	Rcpp::Rcout << "Warming up" << std::endl;

// 	{
// 		Progress prog(ndiscard, true);
// 		for (int idiscard = 0; idiscard < ndiscard; ++idiscard) {
// 			mat DSD = arma::diagmat(1.0 / delta) * S * arma::diagmat(1.0 / delta);
// 			for (int ii = 0; ii < (p*(p-1))/2; ++ii) {
// 				int iR = p - 2 - static_cast<int>(std::sqrt(-8.0*static_cast<double>(ii) + 4.0*static_cast<double>(p*(p-1))-7.0)/2.0 - 0.5); // row index
// 				int iC = ii + iR + 1 - (p*(p-1))/2 + ((p-iR)*((p-iR)-1))/2; // column index
// 				auto fx_rstar = [&](double zstarpinput[])->double {
// 					double zstarp = zstarpinput[0];
// 					double z = std::tanh(zstarp);
// 					vec prhop = pRho;
// 					prhop(ii) = z;
// 					mat pRhop = vecrinv(prhop, p);
// 					pRhop.diag().fill(1.0);
// 					mat Rhop = pRho_to_Rho(pRhop);
// 					mat Rhopinv = arma::inv_sympd(Rhop);
// 					mat RhopChol = arma::chol(Rhop);
// 					double logdet_val = 2.0 * arma::accu(arma::log(RhopChol.diag()));
// 					double loglik = -0.5 * static_cast<double>(n) * logdet_val - 0.5 * arma::dot(Rhopinv, DSD)
// 									+ 0.5 * static_cast<double>(p+1-std::abs(iC-iR)) * std::log1p(-z*z);
// 					return -loglik;
// 				};

// 				double zstar = std::atanh(pRho(ii));
// 				double start[] = { zstar };
// 				double xmin[] = { 0.0 };
// 				double ynewlo = 0.0;
// 				double reqmin = 1.0e-20;
// 				int konvge = 5;
// 				int kcount = 1000;
// 				double step[] = { 0.02 };
// 				int icount = 0;
// 				int numres = 0;
// 				int ifault = 0;
// 				nelmin(fx_rstar, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);					
// 				double xmax = xmin[0];
// 				double minll = ynewlo;

// 				mat cl(5,3, fill::zeros);
// 				vec dl(5, fill::zeros);
// 				double step_size = 0.02;
// 				R_burnin_block:
// 					for (int iii=0; iii < 5; ++iii) {
// 						double e1 = static_cast<double>(iii-2);
// 						cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
// 						cl(iii,1) = xmax + e1 * step_size;
// 						cl(iii,2) = 1.0;
// 						start[0] = xmax + e1 * step_size;
// 						dl(iii) = fx_rstar(start);
// 					}

// 				for (int ni=0; ni < 5; ++ni) {
// 					if ((ni+1) != 3) {
// 						if (dl(ni) <= minll) {
// 							step_size *= 1.2;
// 							goto R_burnin_block;
// 						}
// 					}
// 				}

// 				vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
// 				double sigmaa = std::sqrt(0.5 / fl(0));

// 				double zstar_prop = ::norm_rand() * sigmaa + xmax;
// 				// log-likelihood difference
// 				start[0] = zstar;
// 				double ll_diff = fx_rstar(start);
// 				start[0] = zstar_prop;
// 				ll_diff += -fx_rstar(start)
// 						    - 0.5 * (std::pow(zstar - xmax, 2.0) - std::pow(zstar_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
// 				if (std::log(::unif_rand()) < ll_diff) {
// 					pRho(ii) = std::tanh(zstar_prop);
// 					++pR_rates(ii);
// 				}
// 			}

// 			mat Rhoinv = arma::inv_sympd(Rho);
// 			for (int j = 0; j < p; ++j) {
// 				auto fx_delta = [&](double xistar[])->double {
// 					double xip = xistar[0];
// 					vec delinvp = 1.0 / delta;
// 					delinvp(j) = exp(-xip);
// 					mat DRDinv = arma::diagmat(delinvp) * Rhoinv * arma::diagmat(delinvp);
// 					double loglik = -(static_cast<double>(n) + 0.1) * xip - 0.5 * arma::dot(DRDinv, S) - 0.1 * std::exp(-xip);
// 					return -loglik;
// 				};
// 				double dstar = std::log(delta(j));
// 				double start[] = { dstar };
// 				double xmin[] = { 0.0 };
// 				double ynewlo = 0.0;
// 				double reqmin = 1.0e-20;
// 				int konvge = 5;
// 				int kcount = 1000;
// 				double step[] = { 0.2 };
// 				int icount = 0;
// 				int numres = 0;
// 				int ifault = 0;
// 				nelmin(fx_delta, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);					
// 				double xmax = xmin[0];
// 				double minll = ynewlo;

// 				mat cl(5,3, fill::zeros);
// 				vec dl(5, fill::zeros);
// 				double step_size = 0.2;
// 				delta_burning_block:
// 					for (int iii=0; iii < 5; ++iii) {
// 						double e1 = static_cast<double>(iii-2);
// 						cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
// 						cl(iii,1) = xmax + e1 * step_size;
// 						cl(iii,2) = 1.0;
// 						start[0] = xmax + e1 * step_size;
// 						dl(iii) = fx_delta(start);
// 					}

// 				for (int ni=0; ni < 5; ++ni) {
// 					if ((ni+1) != 3) {
// 						if (dl(ni) <= minll) {
// 							step_size *= 1.2;
// 							goto delta_burning_block;
// 						}
// 					}
// 				}

// 				vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
// 				double sigmaa = std::sqrt(0.5 / fl(0));

// 				// log-likelihood difference
// 				double dstar_prop = ::norm_rand() * sigmaa + xmax;
// 				start[0] = dstar;
// 				double ll_diff = fx_delta(start);
// 				start[0] = dstar_prop;
// 				ll_diff += -fx_delta(start)
// 						    - 0.5 * (std::pow(dstar - xmax, 2.0) - std::pow(dstar_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
// 				if (std::log(::unif_rand()) < ll_diff) {
// 					delta(j) = std::exp(dstar_prop);	
// 					++delta_rates(j);
// 				}
// 			}
// 			prog.increment();
// 		}
// 	}

// 	mat delta_save(p, nkeep, fill::zeros);
// 	cube rho_save(p,p, nkeep, fill::zeros);
// 	cube prho_save(p,p, nkeep, fill::zeros);
// 	Rcpp::Rcout << "Sampling" << std::endl;
// 	{
// 		Progress prog(nkeep, true);
// 		for (int ikeep = 0; ikeep < nkeep; ++ikeep) {
// 			if (Progress::check_abort()) {
// 				return Rcpp::List::create(Rcpp::Named("error") = "user interrupt aborted");
// 			}
// 			for (int iskip = 0; iskip < nskip; ++iskip) {
// 				mat DSD = arma::diagmat(1.0 / delta) * S * arma::diagmat(1.0 / delta);
// 				for (int ii = 0; ii < (p*(p-1))/2; ++ii) {
// 					int iR = p - 2 - static_cast<int>(std::sqrt(-8.0*static_cast<double>(ii) + 4.0*static_cast<double>(p*(p-1))-7.0)/2.0 - 0.5); // row index
// 					int iC = ii + iR + 1 - (p*(p-1))/2 + ((p-iR)*((p-iR)-1))/2; // column index
// 					auto fx_rstar = [&](double zstarpinput[])->double {
// 						double zstarp = zstarpinput[0];
// 						double z = std::tanh(zstarp);
// 						vec prhop = pRho;
// 						prhop(ii) = z;
// 						mat pRhop = vecrinv(prhop, p);
// 						pRhop.diag().fill(1.0);
// 						mat Rhop = pRho_to_Rho(pRhop);
// 						mat Rhopinv = arma::inv_sympd(Rhop);
// 						mat RhopChol = arma::chol(Rhop);
// 						double logdet_val = 2.0 * arma::accu(arma::log(RhopChol.diag()));
// 						double loglik = -0.5 * static_cast<double>(n) * logdet_val - 0.5 * arma::dot(Rhopinv, DSD)
// 										+ 0.5 * static_cast<double>(p+1-std::abs(iC-iR)) * std::log1p(-z*z);
// 						return -loglik;
// 					};

// 					double zstar = std::atanh(pRho(ii));
// 					double start[] = { zstar };
// 					double xmin[] = { 0.0 };
// 					double ynewlo = 0.0;
// 					double reqmin = 1.0e-20;
// 					int konvge = 5;
// 					int kcount = 1000;
// 					double step[] = { 0.02 };
// 					int icount = 0;
// 					int numres = 0;
// 					int ifault = 0;
// 					nelmin(fx_rstar, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);					
// 					double xmax = xmin[0];
// 					double minll = ynewlo;

// 					mat cl(5,3, fill::zeros);
// 					vec dl(5, fill::zeros);
// 					double step_size = 0.02;
// 					R_sampling_block:
// 						for (int iii=0; iii < 5; ++iii) {
// 							double e1 = static_cast<double>(iii-2);
// 							cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
// 							cl(iii,1) = xmax + e1 * step_size;
// 							cl(iii,2) = 1.0;
// 							start[0] = xmax + e1 * step_size;
// 							dl(iii) = fx_rstar(start);
// 						}

// 					for (int ni=0; ni < 5; ++ni) {
// 						if ((ni+1) != 3) {
// 							if (dl(ni) <= minll) {
// 								step_size *= 1.2;
// 								goto R_sampling_block;
// 							}
// 						}
// 					}

// 					vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
// 					double sigmaa = std::sqrt(0.5 / fl(0));

// 					double zstar_prop = ::norm_rand() * sigmaa + xmax;
// 					// log-likelihood difference
// 					start[0] = zstar;
// 					double ll_diff = fx_rstar(start);
// 					start[0] = zstar_prop;
// 					ll_diff += -fx_rstar(start)
// 							    - 0.5 * (std::pow(zstar - xmax, 2.0) - std::pow(zstar_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
// 					if (std::log(::unif_rand()) < ll_diff) {
// 						pRho(ii) = std::tanh(zstar_prop);
// 						++pR_rates(ii);
// 					}
// 				}

// 				mat Rhoinv = arma::inv_sympd(Rho);
// 				for (int j = 0; j < p; ++j) {
// 					auto fx_delta = [&](double xistar[])->double {
// 						double xip = xistar[0];
// 						vec delinvp = 1.0 / delta;
// 						delinvp(j) = exp(-xip);
// 						mat DRDinv = arma::diagmat(delinvp) * Rhoinv * arma::diagmat(delinvp);
// 						double loglik = -(static_cast<double>(n) + 0.1) * xip - 0.5 * arma::dot(DRDinv, S) - 0.1 * std::exp(-xip);
// 						return -loglik;
// 					};
// 					double dstar = std::log(delta(j));
// 					double start[] = { dstar };
// 					double xmin[] = { 0.0 };
// 					double ynewlo = 0.0;
// 					double reqmin = 1.0e-20;
// 					int konvge = 5;
// 					int kcount = 1000;
// 					double step[] = { 0.2 };
// 					int icount = 0;
// 					int numres = 0;
// 					int ifault = 0;
// 					nelmin(fx_delta, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);					
// 					double xmax = xmin[0];
// 					double minll = ynewlo;

// 					mat cl(5,3, fill::zeros);
// 					vec dl(5, fill::zeros);
// 					double step_size = 0.2;
// 					delta_sampling_block:
// 						for (int iii=0; iii < 5; ++iii) {
// 							double e1 = static_cast<double>(iii-2);
// 							cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
// 							cl(iii,1) = xmax + e1 * step_size;
// 							cl(iii,2) = 1.0;
// 							start[0] = xmax + e1 * step_size;
// 							dl(iii) = fx_delta(start);
// 						}

// 					for (int ni=0; ni < 5; ++ni) {
// 						if ((ni+1) != 3) {
// 							if (dl(ni) <= minll) {
// 								step_size *= 1.2;
// 								goto delta_sampling_block;
// 							}
// 						}
// 					}

// 					vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
// 					double sigmaa = std::sqrt(0.5 / fl(0));

// 					// log-likelihood difference
// 					double dstar_prop = ::norm_rand() * sigmaa + xmax;
// 					start[0] = dstar;
// 					double ll_diff = fx_delta(start);
// 					start[0] = dstar_prop;
// 					ll_diff += -fx_delta(start)
// 							    - 0.5 * (std::pow(dstar - xmax, 2.0) - std::pow(dstar_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
// 					if (std::log(::unif_rand()) < ll_diff) {
// 						delta(j) = std::exp(dstar_prop);	
// 						++delta_rates(j);
// 					}
// 				}
// 			}
// 			delta_save.col(ikeep) = delta;
// 			mat prho = vecrinv(pRho, p);
// 			prho.diag().fill(1.0);
// 			prho_save.slice(ikeep) = prho;
// 			rho_save.slice(ikeep) = pRho_to_Rho(prho);
// 			prog.increment();
// 		}
// 	}

// 	return ListBuilder()
// 			.add("delta", delta_save)
// 			.add("rho", rho_save)
// 			.add("prho", prho_save)
// 			.add("delta_rates", delta_rates / static_cast<double>(ndiscard + nskip*nkeep))
// 			.add("rho_rates", pR_rates / static_cast<double>(ndiscard + nskip*nkeep));
// }
