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
				  const arma::mat& VSV) {
	using namespace arma;
	using namespace std;
	using namespace Rcpp;
	using namespace R;

	double z = std::tanh(rstar);
	rowvec tmp_vrtk = vrtk;
	tmp_vrtk(kk) = rstar;

	mat pRR_prop = vecrinv(trans(arma::tanh(tmp_vrtk)), J);
	pRR_prop.diag().fill(1.0);
	mat R_prop = pRho_to_Rho(pRR_prop);
	double logdet_val, logdet_sign;
	arma::log_det(logdet_val, logdet_sign, R_prop);
	double loglik = 0.5 * (ntk - static_cast<double>(J) - 2.0) * logdet_val - 0.5 * (ntk - 1.0) * arma::dot(R_prop, VSV);
	loglik += 0.5 * static_cast<double>(J + 1 - std::abs(iC - iR)) * std::log1p(-z*z);
    return loglik;
}


double loglik_delta_m3(const double& logdel,
					   const arma::rowvec& delta_i,
					   const int& j,
					   const arma::mat& Rhoinv,
					   const arma::mat& qq,
					   const double& a0,
					   const double& b0,
					   const double& ntk) {
	using namespace arma;
	using namespace std;
	using namespace Rcpp;
	using namespace R;

	rowvec siginv = 1.0 / delta_i;
	siginv(j) = std::exp(-logdel);
	mat Vinvp = arma::diagmat(siginv);
	mat VRV = Vinvp * Rhoinv * Vinvp;

	return -0.5 * arma::dot(qq, VRV) - (a0 + ntk) * logdel - b0 * std::exp(-logdel);
}

double loglik_rho_m3(const double& zprho,
					 const arma::vec& vRho,
					 const arma::mat& qq,
					 const int& ii,
					 const int& iR,
					 const int& iC,
					 const int& J,
					 const double& sumNpt) {
	using namespace arma;
	using namespace std;
	using namespace Rcpp;
	using namespace R;

	vec vRhop = vRho;
	vRhop(ii) = zprho;
	double z = std::tanh(zprho);

	mat pRRho = vecrinv(arma::tanh(vRhop), J);
	pRRho.diag().fill(1.0);
	mat Rhop = pRho_to_Rho(pRRho);
	mat Rhopinv = arma::inv(Rhop);
	double logdet_val, logdet_sign;
	arma::log_det(logdet_val, logdet_sign, Rhop);

	double loglik = -0.5 * arma::dot(qq, Rhopinv) - 0.5 * sumNpt * logdet_val;
	loglik += 0.5 * static_cast<double>(J + 1 - std::abs(iC - iR)) * std::log1p(-z*z);
	return loglik;
}

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
					   const arma::mat& Sigma0inv) {
	using namespace arma;
	using namespace std;
	using namespace Rcpp;
	using namespace R;

	vec delta_prop = delta;
	delta_prop(j) = std::exp(logdel);
	int nw = WCovariate.n_cols;
	mat Delta_prop = arma::diagmat(delta_prop);
	mat Sig_prop = Delta_prop * Rho * Delta_prop;
	double loglik = (static_cast<double>(T*K) * nu0 + d0 - static_cast<double>(J)) * logdel - 0.5 * arma::dot(Sigma0inv, Sig_prop);
	for (int i = 0; i < N; ++i) {
		rowvec w_i = WCovariate.row(i);
		mat pRR = vecrinv(trans(arma::tanh(vRtk.row(i))), J);
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
    return loglik;
}

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
	        		 const arma::mat& Sigma0inv) {
	using namespace arma;
	using namespace std;
	using namespace Rcpp;
	using namespace R;

	double z = std::tanh(zprho);
	vec vRhop = vRho;
	vRhop(ii) = zprho;
	mat pRRp = vecrinv(arma::tanh(vRhop), J);
	pRRp.diag().fill(1.0);
	mat Rhop = pRho_to_Rho(pRRp);

	double logdet_val;
	double logdet_sign;
	log_det(logdet_val, logdet_sign, Rhop);

	mat Delta = arma::diagmat(delta);

	int nw = WCovariate.n_cols;
	mat Sig_prop = Delta * Rhop * Delta;
	double loglik = 0.5 * (static_cast<double>(T*K) * nu0 + d0 - static_cast<double>(J) - 1.0) * logdet_val - 0.5 * arma::dot(Sigma0inv, Sig_prop);
	for (int i = 0; i < N; ++i) {
		rowvec w_i = WCovariate.row(i);

		mat pRR = vecrinv(trans(arma::tanh(vRtk.row(i))), J);
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
		log_det(logdet_val, logdet_sign, qq);
		loglik -= 0.5 * (ntk + nu0) * logdet_val;
	}
	loglik += 0.5 * (static_cast<double>(J+1-std::abs(iC-iR))) * std::log1p(-z*z);
    return loglik;
}

