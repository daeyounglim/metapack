#include <stdio.h>
#include <cmath>
#include <RcppArmadillo.h>
#include <Rdefines.h>
#include "linearalgebra.h"
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec vecl(const arma::mat& X) {
	using namespace arma;
	int n = X.n_rows;
	arma::vec out(n*(n-1)/2);
	for (int j = 0; j < n-1; ++j) {
		for (int i = j+1; i < n; ++i) {
			out((n-1)*j - (j-1)*j/2 + i - 1 - j) = X(i, j);
		}
	}
	return out;
}

// [[Rcpp::export]]
arma::mat veclinv(const arma::vec& v, const int& n) {
	using namespace arma;
	mat out(n, n, fill::zeros);
	int count1 = 0;
	int count2 = n-2;
	for (int i = 0; i < n-1; ++i) {
		vec vv = v(span(count1, count2));
		out(span(i+1, n-1), i) = vv;
		out(i, span(i+1,n-1)) = vv.t();
		count1 = count2+1;
		count2 += n-2-i;
	}
	return out;
}


arma::vec vech(const arma::mat& X) {
	int n = X.n_rows;
	arma::vec out(n*(n+1)/2);
	for (int j = 0; j < n; ++j) {
		for (int i = j; i < n; ++i) {
			out(n*j - (j-1)*j/2 + i - j) = X(i, j);
		}
	}
	return out;
}

arma::mat vechinv(const arma::vec& v, const int& n) {
	using namespace arma;
	mat out(n, n, fill::zeros);
	int count1 = 0;
	int count2 = n-1;
	for (int i = 0; i < n-1; ++i) {
		vec vv = v(span(count1, count2));
		out(span(i+1, n-1),i) = vv.tail(n-1-i);
		out(i, span(i, n-1)) = vv.t();
		count1 = count2 + 1;
		count2 += n-1-i;
	}
	out(n-1,n-1) = v(n*(n+1)/2-1);
	return out;
}

arma::mat duplicate_matrix (const int& n) {
  arma::mat mat1 = arma::eye<arma::mat>(n, n);
  arma::vec index = arma::vectorise(arma::linspace<arma::vec>(1, n*(n+1)/2, n*(n+1)/2));
  for (int j = 0; j < n; ++j) {
    for (int i = j; i < n; ++i) {
      if (j == 0) mat1(i, j) = index(i);
      else {
        mat1(i, j) = index(n*j - (j-1)*j/2 + i - j, 0);
      }
    }
  }
  arma::mat mat2 = arma::symmatl(mat1);
  arma::vec temp_vec = arma::vectorise(mat2);
  int t = temp_vec.size();
  int s = index.size();
  arma::mat result(t, s);
  for (int k = 0; k < t; ++k) {
    for (int u = 0; u < s; ++u) {
      if (temp_vec(k) == index(u)) result(k, u) = 1;
      else result(k, u) = 0;
    }
  }
  
  return result;
}

arma::vec uppertriv(const arma::mat& A) {
	using namespace arma;
	int p = A.n_cols;
	umat indmx(p, p, fill::ones);
	for (int i = 0; i < p; ++i) {
		indmx(i,i) = 0;
	}
	indmx = arma::trimatu(indmx);
	return A.elem(find(indmx > 0));
}


arma::mat blockdiag( arma::field<arma::mat>& x ) {
	unsigned int n = x.n_rows;
	int dimen = 0;
	arma::ivec dimvec(n);

	for (unsigned int i = 0; i < n; ++i) {
		dimvec(i) = x(i,0).n_rows; 
		dimen += dimvec(i);
	}

	arma::mat X(dimen, dimen, arma::fill::zeros);
	int idx = 0;

	for (unsigned int i = 0; i < n; ++i) {
		X.submat( idx, idx, idx + dimvec(i) - 1, idx + dimvec(i) - 1 ) = x(i, 0);
		idx += dimvec(i);
	}

	return X;
}

/**************************
Convert partial correlation
to correlation
**************************/
// [[Rcpp::export]]
arma::mat pRho_to_Rho(arma::mat& pRho) {
	using namespace arma;
	using namespace Rcpp;
	using namespace R;
	using namespace std;

	const int nT = pRho.n_rows;
	mat Rho = pRho;

	for (int iL=2; iL < nT; ++iL) {
		mat rmat(iL-1,iL-1, fill::zeros);

		for (int iR=1; iR < nT-iL+1; ++iR) {
			vec rvec1(nT-2, fill::zeros);
			vec rvec3(nT-2, fill::zeros);

			for (int i1=1; i1 < iL; ++i1) {
				rvec1(i1-1) = Rho(iR-1, iR+i1-1);
				rvec3(i1-1) = Rho(iR+i1-1, iR+iL-1);
			}
			double rr11 = 0.0;
			double rr13 = 0.0;
			double rr33 = 0.0;

			for (int i1=1; i1 < iL; ++i1) {
				for (int j1=1; j1 < iL; ++j1) {
					rmat(i1-1,j1-1) = Rho(iR+i1-1, iR+j1-1);
				}
			}
			mat rmatinv = rmat.i();
			for (int i1=1; i1 < iL; ++i1) {
				for (int j1=1; j1 < iL; ++j1) {
					rr11 += rvec1(i1-1) * rmatinv(i1-1,j1-1) * rvec1(j1-1);
					rr13 += rvec1(i1-1) * rmatinv(i1-1,j1-1) * rvec3(j1-1);
					rr33 += rvec3(i1-1) * rmatinv(i1-1,j1-1) * rvec3(j1-1);
				}
			}
			Rho(iR-1, iR+iL-1) = rr13 + pRho(iR-1,iR+iL-1) * std::sqrt((1.0 - rr11) * (1.0 - rr33));
			Rho(iR+iL-1,iR-1) = Rho(iR-1, iR+iL-1);
		}
	}
	return Rho;
}
