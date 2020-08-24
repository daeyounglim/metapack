#include <cmath>
#include <Rmath.h>
#include <algorithm>
#include <iterator>
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <Rdefines.h>
#include "linearalgebra.h"
// [[Rcpp::depends(RcppArmadillo, RcppProgress, BH)]]


/**************************************************
Calculate the goodness of fit measures

+ Dev(theta) = -2 * log L(theta | D_oy)
+ p_D = E[Dev(theta)] - Dev(thetabar)
+ DIC = Dev(thetabar) + 2 * p_D
***************************************************/

double mvnpdf(const arma::vec& x, const arma::vec& mu, const arma::mat& Sig, const bool logp) {
	using namespace arma;

	int k = x.n_elem;
	double logdet_val;
	double logdet_sign;
	log_det(logdet_val, logdet_sign, Sig);
	double lpdf = -static_cast<double>(k) * M_LN_SQRT_2PI	- 0.5 * logdet_val - 0.5 * arma::dot(x - mu, arma::solve(Sig, x - mu));
	if (logp) {
		return lpdf;
	} else {
		return std::exp(lpdf);
	}
}

// [[Rcpp::export]]
Rcpp::List lpml_parcov(const arma::mat& Outcome,
			  		   const arma::mat& XCovariate,
			  		   const arma::mat& WCovariate,
			  		   const arma::vec& Npt,
			  		   const arma::cube& Sigma,
			  		   const arma::cube& Omega,
			  		   const arma::mat& theta,
			  		   const arma::vec& thetahat,
			  		   const arma::mat& Sigmahat,
			  		   const arma::mat& Omegahat,
			  		   const int& fmodel,
			  		   const int& nkeep,
			  		   const bool verbose) {
	using namespace arma;
	const int N = Outcome.n_rows;
	const int J = Outcome.n_cols;
	const int xcols = XCovariate.n_cols;
	const int nw = WCovariate.n_cols;
	mat g(N, nkeep, fill::zeros);
	double alpml = 0.0;
	{
		Progress prog(nkeep, verbose);
		for (int ikeep = 0; ikeep < nkeep; ++ikeep) {
			if (Progress::check_abort()) {
				return Rcpp::List::create(Rcpp::Named("error") = "user interrupt aborted");
			}

			mat Sig_ikeep = Sigma.slice(ikeep);
			mat Omega_ikeep = Omega.slice(ikeep);
			vec theta_ikeep = theta.col(ikeep);

			for (int i = 0; i < N; ++i) {
				double ntk = Npt(i);
				rowvec x_i = XCovariate.row(i);
				rowvec w_i = WCovariate.row(i);
				vec y_i = arma::trans(Outcome.row(i));
				mat X(J, xcols * J, fill::zeros);
				mat W(J, nw*J, fill::zeros);
				for (int j = 0; j < J; ++j) {
					X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
					W(j, span(j*nw, (j+1)*nw-1)) = w_i;
				}
				mat Xstar = arma::join_horiz(X, W);
				if (fmodel >= 3) {
					mat Sig_i = vechinv(arma::trans(Sig_ikeep.row(i)), J);
					mat Q = Sig_i / ntk + W * Omega_ikeep * W.t();
					g(i,ikeep) -= mvnpdf(y_i, Xstar * theta_ikeep, Q, true);
				} else {
					mat Q = Sig_ikeep / ntk + W * Omega_ikeep * W.t();
					g(i,ikeep) -= mvnpdf(y_i, Xstar * theta_ikeep, Q, true);
				}
			}
			prog.increment();
		}

		vec gmax(N, fill::zeros);
		vec alogcpo(N, fill::zeros);
		for (int i = 0; i < N; ++i) {
			gmax(i) = g(i,0);
			for (int j1 = 1; j1 < nkeep; ++j1) {
				if (gmax(i) < g(i, j1)) {
					gmax(i) = g(i, j1);
				}
			}
			double sumrep = 0.0;
			for (int j1 = 1; j1 < nkeep; ++j1) {
				sumrep += std::exp(g(i,j1) - gmax(i));
			}
			alogcpo(i) -= gmax(i) + std::log(sumrep / static_cast<double>(nkeep));
			alpml += alogcpo(i);
		}
	}
	return Rcpp::List::create(Rcpp::Named("lpml")=alpml);
}


// [[Rcpp::export]]
Rcpp::List dic_parcov(const arma::mat& Outcome,
			  		  const arma::mat& XCovariate,
			  		  const arma::mat& WCovariate,
			  		  const arma::vec& Npt,
			  		  const arma::cube& Sigma,
			  		  const arma::cube& Omega,
			  		  const arma::mat& theta,
			  		  const arma::vec& thetahat,
			  		  const arma::mat& Sigmahat,
			  		  const arma::mat& Omegahat,
			  		  const int& fmodel,
			  		  const int& nkeep,
			  		  const bool verbose) {
	using namespace arma;
	double Dev_thetabar = 0.0;
	const int N = Outcome.n_rows;
	const int J = Outcome.n_cols;
	const int xcols = XCovariate.n_cols;
	const int nw = WCovariate.n_cols;
	for (int i = 0; i < N; ++i) {
		double ntk = Npt(i);
		rowvec x_i = XCovariate.row(i);
		rowvec w_i = WCovariate.row(i);
		vec y_i = arma::trans(Outcome.row(i));
		mat X(J, xcols * J, fill::zeros);
		mat W(J, nw*J, fill::zeros);
		for (int j = 0; j < J; ++j) {
			X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
			W(j, span(j*nw, (j+1)*nw-1)) = w_i;
		}
		mat Xstar = arma::join_horiz(X, W);
		if (fmodel >= 3) {
			mat Sighat = vechinv(arma::trans(Sigmahat.row(i)), J);
			mat Q = Sighat / ntk + W * Omegahat * W.t();
			Dev_thetabar -= 2.0 * mvnpdf(y_i, Xstar * thetahat, Q, true);
		} else {
			mat Q = Sigmahat / ntk + W * Omegahat * W.t();
			Dev_thetabar -= 2.0 * mvnpdf(y_i, Xstar * thetahat, Q, true);
		}
	}


	double Dev_bar = 0;
	{
		Progress prog(nkeep, verbose);
		for (int ikeep = 0; ikeep < nkeep; ++ikeep) {
			if (Progress::check_abort()) {
				return Rcpp::List::create(Rcpp::Named("error") = "user interrupt aborted");
			}

			mat Sig_ikeep = Sigma.slice(ikeep);
			mat Omega_ikeep = Omega.slice(ikeep);
			vec theta_ikeep = theta.col(ikeep);

			for (int i = 0; i < N; ++i) {
				double ntk = Npt(i);
				rowvec x_i = XCovariate.row(i);
				rowvec w_i = WCovariate.row(i);
				vec y_i = arma::trans(Outcome.row(i));
				mat X(J, xcols * J, fill::zeros);
				mat W(J, nw*J, fill::zeros);
				for (int j = 0; j < J; ++j) {
					X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
					W(j, span(j*nw, (j+1)*nw-1)) = w_i;
				}
				mat Xstar = arma::join_horiz(X, W);
				if (fmodel >= 3) {
					mat Sig_i = vechinv(arma::trans(Sig_ikeep.row(i)), J);
					mat Q = Sig_i / ntk + W * Omega_ikeep * W.t();
					Dev_bar -= 2.0 * mvnpdf(y_i, Xstar * theta_ikeep, Q, true);
				} else {
					mat Q = Sig_ikeep / ntk + W * Omega_ikeep * W.t();
					Dev_bar -= 2.0 * mvnpdf(y_i, Xstar * theta_ikeep, Q, true);
				}
			}
			prog.increment();
		}
	}
	Dev_bar /= static_cast<double>(nkeep);
	double p_D = Dev_bar - Dev_thetabar;
	double DIC = Dev_thetabar + 2.0 * p_D;
	return Rcpp::List::create(Rcpp::Named("dic")=DIC);
}



