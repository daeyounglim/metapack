#include <stdio.h>
#include <cmath>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <Rdefines.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "random.h"
using namespace std;
using namespace Rcpp;
using namespace R;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]

using namespace std;
using namespace RNG;

/* 
 * rgig.c:
 *
 * Ester Pantaleo and Robert B. Gramacy, 2010
 *
 * adapted from the C code in the monomvm package for R.
 */
#define ZTOL (DOUBLE_EPS*10.0)


double RNG::_gig_mode(const double& lambda, const double& omega)
/*---------------------------------------------------------------------------*/
/* Compute mode of GIG distribution.                                         */
/*                                                                           */
/* Parameters:                                                               */
/*   lambda .. parameter for distribution                                    */
/*   omega ... parameter for distribution                                    */
/*                                                                           */
/* Return:                                                                   */
/*   mode                                                                    */
/*---------------------------------------------------------------------------*/
{
  if (lambda >= 1.)
    /* mode of fgig(x) */
    return (sqrt((lambda-1.)*(lambda-1.) + omega*omega)+(lambda-1.))/omega;
  else
    /* 0 <= lambda < 1: use mode of f(1/x) */
    return omega / (sqrt((1.-lambda)*(1.-lambda) + omega*omega)+(1.-lambda));
} /* end of _gig_mode() */

/*---------------------------------------------------------------------------*/

double RNG::_rgig_ROU_noshift (const double& lambda, const double& lambda_old, const double& omega, const double& alpha)
/*---------------------------------------------------------------------------*/
/* Tpye 1:                                                                   */
/* Ratio-of-uniforms without shift.                                          */
/*   Dagpunar (1988), Sect.~4.6.2                                            */
/*   Lehner (1989)                                                           */
/*---------------------------------------------------------------------------*/
{
  double xm, nc;     /* location of mode; c=log(f(xm)) normalization constant */
  double ym, um;     /* location of maximum of x*sqrt(f(x)); umax of MBR */
  double s, t;       /* auxiliary variables */
  double U, V, X;    /* random variables */

  int count = 0;     /* counter for total number of iterations */

  /* -- Setup -------------------------------------------------------------- */

  /* shortcuts */
  t = 0.5 * (lambda-1.);
  s = 0.25 * omega;
  
  /* mode = location of maximum of sqrt(f(x)) */
  xm = RNG::_gig_mode(lambda, omega);

  /* normalization constant: c = log(sqrt(f(xm))) */
  nc = t*log(xm) - s*(xm + 1./xm);

  /* location of maximum of x*sqrt(f(x)):           */
  /* we need the positive root of                   */
  /*    omega/2*y^2 - (lambda+1)*y - omega/2 = 0    */
  ym = ((lambda+1.) + sqrt((lambda+1.)*(lambda+1.) + omega*omega))/omega;

  /* boundaries of minmal bounding rectangle:                   */
  /* we us the "normalized" density f(x) / f(xm). hence         */
  /* upper boundary: vmax = 1.                                  */
  /* left hand boundary: umin = 0.                              */
  /* right hand boundary: umax = ym * sqrt(f(ym)) / sqrt(f(xm)) */
  um = exp(0.5*(lambda+1.)*log(ym) - s*(ym + 1./ym) - nc);

  /* -- Generate sample ---------------------------------------------------- */


  do {
    ++count;
    U = um * ::unif_rand();        /* U(0,umax) */
    V = ::unif_rand();             /* U(0,vmax) */
    X = U/V;
  }                              /* Acceptance/Rejection */
  while (((log(V)) > (t*log(X) - s*(X + 1./X) - nc)));

  /* store random point */
  return(lambda_old < 0.) ? (alpha / X) : (alpha * X);

}


/*---------------------------------------------------------------------------*/

double RNG::_rgig_newapproach1 (const double& lambda, const double& lambda_old, const double& omega, const double& alpha)
/*---------------------------------------------------------------------------*/
/* Type 4:                                                                   */
/* New approach, constant hat in log-concave part.                           */
/* Draw sample from GIG distribution.                                        */
/*                                                                           */
/* Case: 0 < lambda < 1, 0 < omega < 1                                       */
/*                                                                           */
/* Parameters:                                                               */
/*   n ....... sample size (positive integer)                                */
/*   lambda .. parameter for distribution                                    */
/*   omega ... parameter for distribution                                    */
/*                                                                           */
/* Return:                                                                   */
/*   random sample of size 'n'                                               */
/*---------------------------------------------------------------------------*/
{
  /* parameters for hat function */
  double A[3], Atot;  /* area below hat */
  double k0;          /* maximum of PDF */
  double k1, k2;      /* multiplicative constant */

  double xm;          /* location of mode */
  double x0;          /* splitting point T-concave / T-convex */
  double a;           /* auxiliary variable */

  double U, V, X;     /* random numbers */
  double hx;          /* hat at X */

  int count = 0;      /* counter for total number of iterations */

  /* -- Check arguments ---------------------------------------------------- */

  if (lambda >= 1. || omega >1.)
    throw Rcpp::exception ("invalid parameters");

  /* -- Setup -------------------------------------------------------------- */

  /* mode = location of maximum of sqrt(f(x)) */
  xm = RNG::_gig_mode(lambda, omega);

  /* splitting point */
  x0 = omega/(1.-lambda);

  /* domain [0, x_0] */
  k0 = exp((lambda-1.)*log(xm) - 0.5*omega*(xm + 1./xm));     /* = f(xm) */
  A[0] = k0 * x0;

  /* domain [x_0, Infinity] */
  if (x0 >= 2./omega) {
    k1 = 0.;
    A[1] = 0.;
    k2 = pow(x0, lambda-1.);
    A[2] = k2 * 2. * exp(-omega*x0/2.)/omega;
  }
  
  else {
    /* domain [x_0, 2/omega] */
    k1 = exp(-omega);
    A[1] = (lambda == 0.) 
      ? k1 * log(2./(omega*omega))
      : k1 / lambda * ( pow(2./omega, lambda) - pow(x0, lambda) );

    /* domain [2/omega, Infinity] */
    k2 = pow(2/omega, lambda-1.);
    A[2] = k2 * 2 * exp(-1.)/omega;
  }

  /* total area */
  Atot = A[0] + A[1] + A[2];

  /* -- Generate sample ---------------------------------------------------- */

  do {
    ++count;

    /* get uniform random number */
    V = Atot * ::unif_rand();
    
    do {

      /* domain [0, x_0] */
      if (V <= A[0]) {
        X = x0 * V / A[0];
        hx = k0;
        break;
      }

      /* domain [x_0, 2/omega] */
      V -= A[0];
      if (V <= A[1]) {
        if (lambda == 0.) {
          X = omega * exp(exp(omega)*V);
          hx = k1 / X;
        }
        else {
          X = pow(pow(x0, lambda) + (lambda / k1 * V), 1./lambda);
          hx = k1 * pow(X, lambda-1.);
        }
        break;
      }

      /* domain [max(x0,2/omega), Infinity] */
      V -= A[1];
      a = (x0 > 2./omega) ? x0 : 2./omega;
      X = -2./omega * log(exp(-omega/2. * a) - omega/(2.*k2) * V);
      hx = k2 * exp(-omega/2. * X);
      break;

    } while(false);
    
    /* accept or reject */
    U = ::unif_rand() * hx;

    if (log(U) <= (lambda-1.) * log(X) - omega/2. * (X+1./X)) {
      return (lambda_old < 0.) ? (alpha / X) : (alpha * X);
    }
  } while(true);
    
}


/*---------------------------------------------------------------------------*/

double RNG::_rgig_ROU_shift_alt (const double& lambda, const double& lambda_old, const double& omega, const double& alpha)
/*---------------------------------------------------------------------------*/
/* Type 8:                                                                   */
/* Ratio-of-uniforms with shift by 'mode', alternative implementation.       */
/*   Dagpunar (1989)                                                         */
/*   Lehner (1989)                                                           */
/*---------------------------------------------------------------------------*/
{
  double xm, nc;     /* location of mode; c=log(f(xm)) normalization constant */
  double s, t;       /* auxiliary variables */
  double U, V, X;    /* random variables */

  int count = 0;     /* counter for total number of iterations */

  double a, b, c;    /* coefficent of cubic */
  double p, q;       /* coefficents of depressed cubic */
  double fi, fak;    /* auxiliary results for Cardano's rule */

  double y1, y2;     /* roots of (1/x)*sqrt(f((1/x)+m)) */

  double uplus, uminus;  /* maximum and minimum of x*sqrt(f(x+m)) */

  /* -- Setup -------------------------------------------------------------- */

  /* shortcuts */
  t = 0.5 * (lambda-1.);
  s = 0.25 * omega;

  /* mode = location of maximum of sqrt(f(x)) */
  xm = RNG::_gig_mode(lambda, omega);

  /* normalization constant: c = log(sqrt(f(xm))) */
  nc = t*log(xm) - s*(xm + 1./xm);

  /* location of minimum and maximum of (1/x)*sqrt(f(1/x+m)):  */

  /* compute coeffients of cubic equation y^3+a*y^2+b*y+c=0 */
  a = -(2.*(lambda+1.)/omega + xm);       /* < 0 */
  b = (2.*(lambda-1.)*xm/omega - 1.);
  c = xm;

  /* we need the roots in (0,xm) and (xm,inf) */

  /* substitute y=z-a/3 for depressed cubic equation z^3+p*z+q=0 */
  p = b - a*a/3.;
  q = (2.*a*a*a)/27. - (a*b)/3. + c;

  /* use Cardano's rule */
  fi = acos(-q/(2.*sqrt(-(p*p*p)/27.)));
  fak = 2.*sqrt(-p/3.);
  y1 = fak * cos(fi/3.) - a/3.;
  y2 = fak * cos(fi/3. + 4./3.*M_PI) - a/3.;

  /* boundaries of minmal bounding rectangle:                  */
  /* we us the "normalized" density f(x) / f(xm). hence        */
  /* upper boundary: vmax = 1.                                 */
  /* left hand boundary: uminus = (y2-xm) * sqrt(f(y2)) / sqrt(f(xm)) */
  /* right hand boundary: uplus = (y1-xm) * sqrt(f(y1)) / sqrt(f(xm)) */
  uplus  = (y1-xm) * exp(t*log(y1) - s*(y1 + 1./y1) - nc);
  uminus = (y2-xm) * exp(t*log(y2) - s*(y2 + 1./y2) - nc);

  /* -- Generate sample ---------------------------------------------------- */

  do {
    ++count;
    U = uminus + ::unif_rand() * (uplus - uminus);    /* U(u-,u+)  */
    V = ::unif_rand();                                /* U(0,vmax) */
    X = U/V + xm;
  }                                         /* Acceptance/Rejection */
  while ((X <= 0.) || ((log(V)) > (t*log(X) - s*(X + 1./X) - nc)));

    /* store random point */
  return (lambda_old < 0.) ? (alpha / X) : (alpha * X);
}


double RNG::gigrnd(const double& lambda, const double& chi, const double& psi) {
  double omega, alpha;
  if ( !(R_FINITE(lambda) && R_FINITE(chi) && R_FINITE(psi)) ||
       (chi <  0. || psi < 0)      || 
       (chi == 0. && lambda <= 0.) ||
       (psi == 0. && lambda >= 0.) ) {
    throw Rcpp::exception("invalid parameters for GIG distribution");
  }

  if (chi < ZTOL) { 
    /* special cases which are basically Gamma and Inverse Gamma distribution */
    if (lambda > 0.0) {
      return R::rgamma(lambda, 2.0/psi); 
    }
    else {
      return 1.0/R::rgamma(-lambda, 2.0/psi); 
    }    
  }

  else if (psi < ZTOL) {
    /* special cases which are basically Gamma and Inverse Gamma distribution */
    if (lambda > 0.0) {
      return 1.0/R::rgamma(lambda, 2.0/chi); 
    }
    else {
      return R::rgamma(-lambda, 2.0/chi); 
    }    

  }

  else {
    double lam = lambda;
    double lambda_old = lam;
    if (lam < 0.) {lam = -lam;}
    alpha = sqrt(chi/psi);
    omega = sqrt(psi*chi);

    do {
      if (lam > 2. || omega > 3.) {
        return RNG::_rgig_ROU_shift_alt(lam, lambda_old, omega, alpha);
      }

      if (lam >= 1.-2.25*omega*omega || omega > 0.2) {
        return RNG::_rgig_ROU_noshift(lam, lambda_old, omega, alpha);
      }

      if (lam >= 0. && omega > 0.) {
        return RNG::_rgig_newapproach1(lam, lambda_old, omega, alpha);
      }
      throw Rcpp::exception("parameters must satisfy lambda>=0 and omega>0.");
      
    } while (true);
  }
}

double RNG::besselM3(const double& lambda, const double& x, const bool& logvalue)
{
	double res = 0.0;
	if (abs(lambda) == 0.5) {
		if (!logvalue) {
			res = sqrt(M_PI / (2.0 * x)) * exp(-x);
		} else {
			res = 0.5 * log(M_PI / (2.0 * x)) - x;
		}
	} else {
		if (!logvalue) {
			res = ::Rf_bessel_k(x, lambda, 1.0);
		} else {
			res = log( ::Rf_bessel_k(x, lambda, 2.0)) - x;
			// res = ::Rf_bessel_k(x, lambda, );
		}
	}
	return res;
}

/*******************************
Calculate E(X)
where X ~ GIG(lambda, chi, psi)
*******************************/
double RNG::EGIG_x(const double& lambda, const double& chi, const double& psi)
{
	if (psi == 0.0) {
		// Inv Gamma -> Student-t
		double beta = 0.5 * chi;
		double alpha = -lambda;
		if (alpha <= 1.0) {alpha = NA_REAL;}
		return beta / (alpha - 1.0);
	} else if (chi == 0.0) {
		// Gamma -> VG
		double beta = 0.5 * psi;
		double alpha = lambda;
		return alpha / beta;
	} else {
		double chi_ = abs(chi) < datum::eps ? sign(chi) * datum::eps : chi;
		if (chi_ == 0.0) {chi_ = datum::eps;}

		double alpha_bar = sqrt(chi_ * psi);
		double term1 = 0.5 * log(chi_ / psi);
		double term2 = besselM3(lambda + 1.0, alpha_bar, true);
		double term3 = besselM3(lambda, alpha_bar, true);
		return exp(term1 + term2 - term3);
	}
}

/*******************************
Calculate E(1/X)
where X ~ GIG(lambda, chi, psi)
*******************************/
double RNG::EGIG_xinv(const double& lambda, const double& chi, const double& psi)
{
	if (psi == 0.0) {
		double beta = 0.5 * chi;
		double alpha = -lambda;
		return alpha / beta;
	} else if (chi == 0.0) {
		Rcpp::warning("Case 'chi == 0' and 'func = 1/x' is not implemented");
		return NA_REAL;
	} else {
		// GIG -> ghyp, hyp, NIG
		double chi_ = abs(chi) < datum::eps ? sign(chi) * datum::eps : chi;
		if (chi_ == 0.0) {chi_ = datum::eps;}

		double alpha_bar = sqrt(chi_ * psi);
		double term1 = -0.5 * log(chi_ / psi);
		double term2 = besselM3(lambda - 1.0, alpha_bar, true);
		double term3 = besselM3(lambda, alpha_bar, true);
		return exp(term1 + term2 - term3);
	}
}


double RNG::rtgamma(const double& a, const double& b, const double& truncpoint, const bool& up)
{
  double x = 0.;
  if (truncpoint < 0.) {
    if (up) {
      x = R::rgamma(a,1./b);
    } else {
      Rf_error("truncation point is negative");
    }
  }
  if (up) {
    double a2 = a;
    double gup = R::pgamma(b*truncpoint,a,1.,1,0);
    if (gup >= 1.) {
      gup = 0.99999;
    }
    double p = R::runif(0.,1.)*gup;
    x = R::qgamma(p,a2,1.,1,0)/b;
    x = (truncpoint-x)*(x>truncpoint)+x;
  } else {
    double a2 = a;
    double glow = R::pgamma(b*truncpoint,a,1.,1,0);
    if (glow >= 0.9999) {
      x = truncpoint-log1p(-R::runif(0.,1.))/b;
    } else {
      double p = R::runif(0.,1.)*(1.-glow)+glow;
      x = R::qgamma(p,a2,1.,1,0)/b;
    }
    x = (truncpoint-x)*(x<truncpoint)+x;
  }
  return x;
}

arma::mat RNG::rwish(const double& v, const arma::mat& S) {
  int p = S.n_rows;
  arma::mat R = arma::chol(S);
  arma::mat A(p,p);
  // std::generate(A.begin(), A.end(), ::norm_rand);
  for (int i = 0; i < p-1; ++i) {
    for (int j = i+1; j < p; ++j) {
      A(j,i) = ::norm_rand();
    }
  }
  // A = arma::trimatl(A);
  
  for (int i = 0; i < p; ++i) {
    A(i,i) = std::sqrt(::Rf_rchisq(v-(double)i));
  }
  return R.t() * A * A.t() * R;
}

arma::mat RNG::riwish(const double& v, const arma::mat& S) {
  // int p = S.n_rows;
  // arma::mat R = arma::chol(S.i());
  // arma::mat A(p,p);
  // for (int i = 0; i < p-1; ++i) {
  //   for (int j = i+1; j < p; ++j) {
  //     A(j,i) = ::norm_rand();
  //   }
  // }
  // // std::generate(A.begin(), A.end(), ::norm_rand);
  // // A = arma::trimatl(A);
  
  // for (int i = 0; i < p; ++i) {
  //   A(i,i) = std::sqrt(::Rf_rchisq(v-(double)i));
  // }
  // return arma::inv_sympd(R.t() * A * A.t() * R);
  arma::mat Sinv = arma::inv_sympd(S);
  arma::mat out = rwish(v, Sinv);
  return arma::inv_sympd(out);
}

double RNG::tnormrnd(const double& mu, const double& sigma, const double& low, const double& up) {
    double u, pleft, pright, y;

    pleft=pnorm(low,mu,sigma,1,0);
    pright=pnorm(up,mu,sigma,1,0);
    if (pleft>0.9999)
        return (low)+0.0001*fmax2((up)-(low),sigma);
    if (pright<0.0001)
        return (up)-0.0001*fmax2((up)-(low),sigma);
    u=::unif_rand();
    u=pleft+(pright-pleft)*u;
    y=R::qnorm5(u,mu,sigma,1,0);

    return y;
}

double RNG::rtnormrnd(const double& mu, const double& sigma, const double& up)
{
    double u,pcut,x;

    if ((sigma)==0.0){ /* degenerate sigma=0 */
        return ((mu)<(up)) ? (mu) : (up);
    }
    pcut=R::pnorm5(up,mu,sigma,1,0);
    if (pcut < 0.0001)
        return (up)-0.0001*(sigma);
    u=::unif_rand();
    u=u*pcut;
    x=R::qnorm5(u,mu,sigma,1,0);

    return x;
}

double RNG::ltnormrnd(const double& mu, const double& sigma, const double& low)
{
    double u,pcut,x;

    if ((sigma)==0.0){ /* degenerate sigma=0 */
        return ((mu)>(low)) ? (mu) : (low);
    }

    pcut=R::pnorm5(low,mu,sigma,1,0);
    if (pcut>0.9999)
        return (low)+0.0001*(sigma);
    u=::unif_rand();
    u=pcut+(1.0-pcut)*u;
    x=R::qnorm5(u,mu,sigma,1,0);

    return x;
}
