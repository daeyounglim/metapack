#include <iostream>
#include <cmath>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <Rdefines.h>
#include "random.h"
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat rtmvn_gibbs(const int& n, const int& p, const arma::vec& Mean, const arma::mat& Sigma_chol,
                      const arma::mat& R, const arma::vec& a, const arma::vec& b, arma::vec& z)
{
  arma::mat keep_x(p, n);
  arma::vec temp;
  
  for(int i = 0; i < n; ++i)
  {
    for(int j = 0; j < p; ++j)
    {
      // loop to create indicies for subsetting
      // 1 to p, excluding j
      arma::uvec idx(p - 1);
      for(int k = 0; k < p; ++k)
      {
        if(k == j)
          continue;
        else if (k > j)
          idx(k - 1) = k;
        else
          idx(k) = k;
      }
      
      // subsetting the matricies
      arma::vec rj = R.col(j);
      arma::mat Rj = R.cols(idx);
      arma::vec zj = z.rows(idx);
      arma::vec Rj_zj = Rj * zj;
      
      arma::vec a_temp = a - Rj_zj;
      arma::vec b_temp = b - Rj_zj;
      
      arma::uvec pos = rj > 0;
      arma::uvec neg = rj < 0;
      
      // which are negtive and positive
      arma::uvec neg_idx = find(neg == 1);
      arma::uvec pos_idx = find(pos == 1);
      
      double lower_pos, upper_pos, lower_neg, upper_neg;
      if(sum(pos) == 0)
      {
        lower_pos = R_NegInf;
        upper_pos = R_PosInf;
      }
      else
      {
        temp = a_temp.rows(pos_idx) / rj.rows(pos_idx);
        lower_pos = max(temp);
        
        temp = b_temp.rows(pos_idx) / rj.rows(pos_idx);
        upper_pos = min(temp);
      }
      
      if(sum(neg) == 0)
      {
        upper_neg = R_PosInf;
        lower_neg = R_NegInf;
      }
      else
      {
        temp = a_temp.rows(neg_idx) / rj.rows(neg_idx);
        upper_neg = min(temp);
        
        temp = b_temp.rows(neg_idx) / rj.rows(neg_idx);
        lower_neg = max(temp);
      }
      
      double lower_j = std::max(lower_pos, lower_neg);
      double upper_j = std::min(upper_pos, upper_neg);
      
      z(j) = RNG::tnormrnd(0., 1., lower_j, upper_j);
      
      // cout << "Z_j: " << z(j) << endl;
    }
    keep_x.col(i) = Sigma_chol * z + Mean;
  }
  return(keep_x);
}

arma::mat rtmvn(const int& n, const arma::vec& Mean, const arma::mat& Sigma, const arma::mat& D,
                   const arma::vec& lower, const arma::vec& upper, const arma::vec& init) {
  /******************************
  Test whether the initial value
  satisfies the linear constraint
  ******************************/
  arma::vec inits_test = D * init;
  if (arma::any(inits_test <= lower) || arma::any(inits_test >= upper)) {
    Rf_error("initial value outside bounds");
  }
  
  arma::mat Rtilde = D;
  arma::vec a = lower - Rtilde * Mean;
  arma::vec b = upper - Rtilde * Mean;
  arma::mat Sigma_chol = arma::chol(Sigma, "lower");
  arma::mat R = Rtilde * Sigma_chol;
  int p = R.n_cols;
  arma::vec z = arma::solve(arma::trimatl(Sigma_chol), init - Mean);
  return rtmvn_gibbs(n, p, Mean, Sigma_chol, R, a, b, z);
}