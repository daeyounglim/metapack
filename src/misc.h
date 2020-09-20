#ifndef BAYESMETA_MISC_H
#define BAYESMETA_MISC_H

#include <stdio.h>
#include <cmath>
#include <RcppArmadillo.h>
#include <Rdefines.h>
// [[Rcpp::depends(RcppArmadillo)]]

double logSumExp(const double& a, const double& b);
double log1exp(const double& z);

#endif