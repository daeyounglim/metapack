#ifndef BAYESMETA_LINEARALGEBRA_H
#define BAYESMETA_LINEARALGEBRA_H
#include <stdio.h>
#include <cmath>
#include <RcppArmadillo.h>
#include <Rdefines.h>
#include "linearalgebra.h"
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec vecl(const arma::mat& X);

template <typename TN>
arma::mat veclinv(TN& v, const int& J) {
	using namespace arma;
	const int vdim = v.n_elem;
	mat out(J, J, fill::zeros);
	for (int kk = 0; kk < vdim; ++kk) {
		int iR = J - 2 - static_cast<int>(std::sqrt(-8.0*static_cast<double>(kk) + 4.0*static_cast<double>(J*(J-1))-7.0)/2.0 - 0.5); // row index
		int iC = kk + iR + 1 - (J*(J-1))/2 + ((J-iR)*((J-iR)-1))/2; // column index
		out(iC,iR) = v(kk);
	}
	return out;
};
arma::vec vech(const arma::mat& X);


arma::mat vechinv(const arma::vec& v, const int& n);
arma::mat duplicate_matrix (const int& n);
arma::vec uppertriv(const arma::mat& A);
arma::mat blockdiag( arma::field<arma::mat>& x );
arma::mat pRho_to_Rho(arma::mat& pRho);
arma::mat constructR(const arma::vec& vphi, const int& J);
arma::mat vecrinv(const arma::vec& X, const int& J);
// template <typename TN>
// arma::mat vecrinv(TN& X, const int& J) {
// 	const int vdim = X.n_elem;
// 	arma::mat R(J, J, arma::fill::zeros);
// 	for (int kk = 0; kk < vdim; ++kk) {
// 		int iR = J - 2 - static_cast<int>(std::sqrt(-8.0*static_cast<double>(kk) + 4.0*static_cast<double>(J*(J-1))-7.0)/2.0 - 0.5); // row index
// 		int iC = kk + iR + 1 - (J*(J-1))/2 + ((J-iR)*((J-iR)-1))/2; // column index
// 		R(iR,iC) = X(kk);
// 		R(iC,iR) = X(kk);
// 	}
// 	return R;
// }
arma::mat vecr(const arma::mat& X);
arma::mat constructB(const arma::mat& Rangle);
#endif