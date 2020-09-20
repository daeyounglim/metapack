#include <stdio.h>
#include <cmath>
#include <RcppArmadillo.h>
#include <Rdefines.h>
#include "misc.h"
// [[Rcpp::depends(RcppArmadillo)]]

double logSumExp(const double& a, const double& b) {
	double a_ = std::log(a);
	double b_ = std::log(b);
  if (a_ > b_) {
  	return a_ + std::log1p(std::exp(b_ - a_));
  } else 
  	return b_ + std::log1p(std::exp(a_ - b_));
}


double log1exp(const double& z) {
	// computes log(1 + exp(-z))
	if (-z > 0) {
		return -z + std::log1p(std::exp(z));
	} else {
		return std::log1p(std::exp(-z));
	}
}