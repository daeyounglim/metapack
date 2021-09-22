#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <algorithm>
#include <iterator>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "linearalgebra.h"
#include "loglik_POCov.h"
#include "ramcmc.h"
#include "nelmin.h"
// [[Rcpp::depends(RcppArmadillo,RcppProgress)]]

double loglik_delta_mcorr(const double &logdel,
                    const arma::rowvec &delta_i,
                    const int &j,
                    const arma::mat &Rhoinv_tk,
                    const arma::mat &qq, // = ntk * eeT + (ntk - 1.0) * Stk;
                    const double &a0,
                    const double &b0,
                    const double &ntk)
{
    using namespace arma;
    using namespace std;
    using namespace Rcpp;
    using namespace R;

    rowvec siginv = 1.0 / delta_i;
    siginv(j) = std::exp(-logdel);
    mat Vinvp = arma::diagmat(siginv);
    mat VRV = Vinvp * Rhoinv_tk * Vinvp; // sig_{tk}^{-1} * rho_{tk}^{-1} * sig_{tk}^{-1} (not the sample covariance matrix)

    return -0.5 * arma::dot(qq, VRV) - (a0 + ntk) * logdel - b0 * std::exp(-logdel);
}

double loglik_eta_mcorr(const arma::vec &eta,
                  const arma::mat &WCovariate,
                  const arma::mat &ZCovariate,
                  const arma::uvec &Trial,
                  const arma::uvec &Second,
                  const arma::mat &gamR,
                  const arma::mat &delta,
                  const arma::mat &resid,
                  const arma::vec &Npt,
                  const arma::mat &S,
                  const arma::vec &phi_all,
                  const int &transform_type,
                  const arma::mat &I_JJm12,
                  const double &TOL)
{
    using namespace arma;
    using namespace std;
    using namespace Rcpp;
    using namespace R;

    int N = Npt.n_elem;
    int J = delta.n_cols;
    int nn = WCovariate.n_cols;
    int nw = nn * 2;

    double loglik = -0.5 * arma::accu(arma::square(eta) / phi_all);
    for (int i = 0; i < N; ++i)
    {
        int k = Trial(i);
        vec gam_k = gamR.col(k);
        rowvec w_i = WCovariate.row(i);
        rowvec wstar_i(nw, fill::zeros);
        if (Second(i) == 0)
        {
            wstar_i.head(nn) = w_i;
        }
        else
        {
            wstar_i.tail(nn) = w_i;
        }
        mat W(J, nw * J, fill::zeros);
        for (int j = 0; j < J; ++j)
        {
            W(j, arma::span(j * nw, (j + 1) * nw - 1)) = wstar_i;
        }
        vec resid_i = arma::trans(resid.row(i)) - W * gam_k;

        rowvec z_tk = ZCovariate.row(i);
        mat Zmat = arma::kron(I_JJm12, z_tk);
        vec xi_tk = Zmat * eta;

        // Calculate rho_tk from xi_tk
        // Store the off-diagonal elements as a vector in vRho
        mat rho_tk(J, J, fill::zeros);
        if (transform_type == 1)
        {
            // rho is modeled as a matrix exponential of vecrinv(xi_tk)
            mat varrho_tk = find_diag(xi_tk, TOL);
            rho_tk = arma::expmat_sym(varrho_tk);
            rho_tk.diag().fill(1.0);
        }
        else if (transform_type == 2)
        {
            // rho is modeled as recovered from the partial correlation matrix
            // partial correlation matrix is further modeled as tanh(vecrinv(xi_tk))
            vec tanh_xi_tk = arma::tanh(xi_tk);
            mat pRho_tk = vecrinv(tanh_xi_tk, J);
            rho_tk = pRho_to_Rho(pRho_tk);
        }
        else if (transform_type == 3)
        {
            // rho is modeled by its Cholesky factor
            // The Cholesky factor is further modeled as hyperspherical coordinates (or hyperangles)
            mat tmp = M_PI / (1.0 + arma::exp(-xi_tk));
            mat varrho_tk = vecrinv(tmp, J);
            mat B_tk = constructB(varrho_tk);
            rho_tk = arma::trimatl(B_tk) * arma::trimatu(B_tk.t());
        }

        double logdet_val, logdet_sign;
        arma::log_det(logdet_val, logdet_sign, rho_tk);

        mat siginv = arma::diagmat(1.0 / delta.row(i));
        mat Siginv = siginv * rho_tk * siginv;
        rowvec S_i = S.row(i);
        mat qq = Npt(i) * resid_i * resid_i.t() + (Npt(i) - 1.0) * vechinv(S_i, J);
        loglik -= 0.5 * (Npt(i) * logdet_val + arma::dot(Siginv, qq));
    }

    return loglik;
}

// [[Rcpp::export]]
Rcpp::List fmodel_corr_modeling(const arma::mat &Outcome,
                                const arma::mat &SD,
                                const arma::mat &XCovariate,
                                const arma::mat &WCovariate,
                                const arma::mat &ZCovariate,
                                const arma::uvec &Treat,
                                const arma::uvec &Trial,
                                const arma::uvec &Second,
                                const arma::vec &Npt,
                                const double &c0,
                                const double &dj0, // hyperparameter for Omega
                                const double &a0,  // hyperparameter for delta
                                const double &b0,  // hyperparameter for delta
                                const double &a1,  // hyperparameter for phi
                                const double &b1,  // hyperparameter for phi
                                const double &a2,  // hyperparameter for phi_{ell(i)}
                                const double &b2,  // hyperparameter for phi_{ell(i)}
                                const double &a3,  // hyperparameter for phi_i
                                const double &b3,  // hyperparameter for phi_i
                                const arma::mat &Omega0,
                                const int &K, // # of Trials
                                const int &T, // # of Treatments
                                const int &ndiscard,
                                const int &nskip,
                                const int &nkeep,
                                const int &transform_type, // 1: matrix exponential; 2: partial correlation; 3: hyperspherical coordinates
                                const double &delta_stepsize,
                                const double &R_stepsize,
                                const double &TOL,
                                const arma::vec &theta_init,
                                const arma::mat &gamR_init,
                                const arma::mat &Omega_init,
                                const bool &verbose)
{
    using namespace arma;
    using namespace std;
    using namespace Rcpp;
    using namespace R;

    const int N = Outcome.n_rows;
    const int J = Outcome.n_cols;
    const int xcols = XCovariate.n_cols;
    const int zcols = ZCovariate.n_cols;
    const int nn = WCovariate.n_cols;
    const int nw = nn * 2;
    const int nt = (xcols + nw) * J;
    const int JJm12 = (J * (J - 1)) / 2; // J * (J minus 1) / 2 => JJm12
    const int JJp12 = (J * (J + 1)) / 2; // J * (J plus 1) / 2 => JJp12
    const int neta = JJm12 * zcols;

    arma::field<arma::uvec> idxks(K);
    vec onstat(K, fill::zeros);
    for (int k = 0; k < K; ++k)
    {
        uvec idx = find(Trial == k);
        idxks(k) = idx;
        int i_k = idx(0);
        onstat(k) = static_cast<double>(Second(i_k));
    }

    /***********************
	Parameter Initialization
	***********************/
    vec theta = theta_init;
    mat gamR = gamR_init;
    mat Omega = Omega_init;
    mat Omegainv = Omega.i();
    mat Sig_lt(N, JJp12, fill::zeros);    // store the diagonal-inclusive lower triangular (lt) elements of Sig
    mat Siginv_lt(N, JJp12, fill::zeros); // store the diagonal-inclusive lower triangular (lt) elements of Siginv
    mat vRtk(N, JJm12, fill::zeros);      // store the off-diagonal lower triangular elements of normal variates for Rtk
    vRtk.fill(0.5);
    mat S(N, JJp12, fill::zeros);
    vec phi_group(JJm12, fill::ones);
    vec phi_local(neta, fill::ones);
    double phi_global = 1.0;

    vec phi_all = phi_local * phi_global;
    for (int i = 0; i < JJm12; ++i) {
        phi_all(arma::span(i*zcols,(i+1)*zcols-1)) *= phi_group(i);
    }
    vec eta(neta, fill::ones);
    mat I_JJm12 = arma::eye(JJm12, JJm12);
    mat xi(N, JJm12, fill::zeros);
    mat Rho(N, JJm12, fill::zeros);
    mat Rhoinv(N, JJm12, fill::zeros);
    mat delta = SD;
    for (int i = 0; i < N; ++i) {
        rowvec z_tk = ZCovariate.row(i);
        mat Zmat = arma::kron(I_JJm12, z_tk);
        vec xi_tk = Zmat * eta;
        xi.row(i) = xi_tk.t();

        // Calculate rho_tk from xi_tk
        // Store the off-diagonal elements as a vector in vRho
        mat rho_tk(J, J, fill::zeros);
        if (transform_type == 1) {
            // rho is modeled as a matrix exponential of vecrinv(xi_tk)
            mat varrho_tk = find_diag(xi_tk, TOL);
            rho_tk = arma::expmat_sym(varrho_tk);
            rho_tk.diag().fill(1.0);
        } else if (transform_type == 2) {
            // rho is modeled as recovered from the partial correlation matrix
            // partial correlation matrix is further modeled as tanh(vecrinv(xi_tk))
            vec tanh_xi_tk = arma::tanh(xi_tk);
            mat pRho_tk = vecrinv(tanh_xi_tk, J);
            rho_tk = pRho_to_Rho(pRho_tk);
        } else if (transform_type == 3) {
            // rho is modeled by its Cholesky factor
            // The Cholesky factor is further modeled as hyperspherical coordinates (or hyperangles)
            mat tmp = M_PI / (1.0 + arma::exp(-xi_tk));
            mat varrho_tk = vecrinv(tmp, J);
            mat B_tk = constructB(varrho_tk);
            rho_tk = arma::trimatl(B_tk) * arma::trimatu(B_tk.t());
        }
        Rho.row(i) = arma::trans(vecr(rho_tk));
        Rhoinv.row(i) = arma::trans(vecr(rho_tk.i()));

        // Use rho_tk to initialize Sig_lt
        mat Sigma = arma::diagmat(delta.row(i)) * rho_tk * arma::diagmat(delta.row(i));
        mat Sigmainv = arma::diagmat(1.0 / delta.row(i)) * rho_tk.i() * arma::diagmat(1.0 / delta.row(i));
        Sig_lt.row(i) = arma::trans(vech(Sigma));
        Siginv_lt.row(i) = arma::trans(vech(Sigmainv));

        // Store S matrix
        rowvec s_tk = SD.row(i);
        mat Vtk = arma::diagmat(s_tk);
        mat pRtk = vecrinv(arma::trans(arma::tanh(vRtk.row(i))), J);
        pRtk.diag().fill(1.0);
        mat Rtk = pRho_to_Rho(pRtk);
        S.row(i) = arma::trans(vech(Vtk * Rtk * Vtk));
    }


    const mat Omega0inv = arma::inv(Omega0);
    const double K2 = arma::accu(onstat);
    const double K1 = static_cast<double>(K) - K2;
    const double shape_omega1 = K1 + dj0;
    const double shape_omega2 = K2 + dj0;
    const double a2_group = a2 + 0.5 * static_cast<double>(zcols);
    const double a3_local = a3 + 0.5;
    mat resid = Outcome;
    mat delta_rates(arma::size(delta), fill::zeros);
    double eta_rates = 0.0;
    mat vR_rates(N, JJm12, fill::zeros);

    const double alpha_star = 0.234;
    const double gam_exp = 2.0 / 3.0;
    mat SS(neta, neta, fill::eye);
    /*********
	Containers
	*********/
    // cube Omega_save(nw * J, nw * J, nkeep, fill::zeros);
    // cube Rho_save(N, JJm12, nkeep, fill::zeros);
    mat theta_save(nt, nkeep, fill::zeros);
    cube Rtk_save(N, JJm12, nkeep, fill::zeros);
    cube resid_save(N, J, nkeep, fill::zeros);
    cube delta_save(N, J, nkeep, fill::zeros);
    mat eta_save(neta, nkeep, fill::zeros);
    int icount_mh = 0;
    /*******************
	Begin burn-in period
	*******************/
    if (verbose)
    {
        Rcout << "Warming up" << endl;
    }
    {
        Progress prog(ndiscard, verbose);
        for (int idiscard = 0; idiscard < ndiscard; ++idiscard)
        {
            if (Progress::check_abort())
            {
                return Rcpp::List::create(Rcpp::Named("error") = "user interrupt aborted");
            }
            ++icount_mh;

            // Update theta
            Rcpp::Rcout << "theta" << std::endl;
            mat Sig_theta(nt, nt, fill::zeros);
            Sig_theta.diag().fill(1.0 / c0);
            vec mu_theta(nt, fill::zeros);
            for (int k = 0; k < K; ++k)
            {
                uvec idxk = idxks(k);
                int n_k = idxk.n_elem;
                mat XSX(nt, nt, fill::zeros);
                mat WSX(nw * J, nt, fill::zeros);
                vec WSy(nw * J, fill::zeros);
                vec XSy(nt, fill::zeros);
                mat Sig_gamk = Omegainv;
                for (int i = 0; i < n_k; ++i)
                {
                    int i_k = idxk(i);
                    rowvec x_i = XCovariate.row(i_k);
                    rowvec w_i = WCovariate.row(i_k);
                    rowvec wstar_i(nw, fill::zeros);
                    if (Second(i_k) == 0)
                    {
                        wstar_i.head(nn) = w_i;
                    }
                    else
                    {
                        wstar_i.tail(nn) = w_i;
                    }
                    rowvec y_i = Outcome.row(i_k);
                    double ntk = Npt(i_k);
                    mat Siginv = vechinv(arma::trans(Siginv_lt.row(i_k)), J);
                    mat X(J, xcols * J, fill::zeros);
                    mat W(J, nw * J, fill::zeros);
                    for (int j = 0; j < J; ++j)
                    {
                        X(j, arma::span(j * xcols, (j + 1) * xcols - 1)) = x_i;
                        W(j, arma::span(j * nw, (j + 1) * nw - 1)) = wstar_i;
                    }
                    mat Xstar = arma::join_horiz(X, W);
                    XSX += ntk * (Xstar.t() * Siginv * Xstar);
                    XSy += ntk * (Xstar.t() * Siginv * y_i.t());
                    mat WS = ntk * (W.t() * Siginv);
                    Sig_gamk += WS * W;
                    WSX += WS * Xstar;
                    WSy += WS * y_i.t();
                }
                mat Sig_gamk_inv = arma::inv(Sig_gamk);
                Sig_theta += XSX - WSX.t() * Sig_gamk_inv * WSX;
                mu_theta += XSy - WSX.t() * Sig_gamk_inv * WSy;
            }
            Sig_theta = 0.5 * (Sig_theta + Sig_theta.t());
            mat Sig_theta_chol = arma::chol(Sig_theta);
            vec ttmp(nt);
            std::generate(ttmp.begin(), ttmp.end(), ::norm_rand);
            theta = arma::solve(arma::trimatu(Sig_theta_chol), arma::solve(arma::trimatl(Sig_theta_chol.t()), mu_theta) + ttmp);

            for (int i = 0; i < N; ++i)
            {
                rowvec x_i = XCovariate.row(i);
                rowvec w_i = WCovariate.row(i);
                rowvec y_i = Outcome.row(i);
                rowvec wstar_i(nw, fill::zeros);
                if (Second(i) == 0)
                {
                    wstar_i.head(nn) = w_i;
                }
                else
                {
                    wstar_i.tail(nn) = w_i;
                }
                mat X(J, xcols * J, fill::zeros);
                mat W(J, nw * J, fill::zeros);
                for (int j = 0; j < J; ++j)
                {
                    X(j, arma::span(j * xcols, (j + 1) * xcols - 1)) = x_i;
                    W(j, arma::span(j * nw, (j + 1) * nw - 1)) = wstar_i;
                }
                mat Xstar = arma::join_horiz(X, W);
                resid.row(i) = arma::trans(y_i.t() - Xstar * theta);
            }

            // Update gamR
            Rcpp::Rcout << "gamR" << std::endl;
            for (int k = 0; k < K; ++k)
            {
                mat Siggam = Omegainv;
                vec mugam(nw * J, fill::zeros);
                uvec idxk = idxks(k);
                int n_k = idxk.n_elem;
                for (int i = 0; i < n_k; ++i)
                {
                    int i_k = idxk(i);
                    rowvec w_i = WCovariate.row(i_k);
                    rowvec wstar_i(nw, fill::zeros);
                    if (Second(i_k) == 0)
                    {
                        wstar_i.head(nn) = w_i;
                    }
                    else
                    {
                        wstar_i.tail(nn) = w_i;
                    }
                    rowvec resid_i = resid.row(i_k);
                    double ntk = Npt(i_k);
                    mat Siginv = vechinv(arma::trans(Siginv_lt.row(i_k)), J);
                    mat W(J, nw * J, fill::zeros);
                    for (int j = 0; j < J; ++j)
                    {
                        W(j, arma::span(j * nw, (j + 1) * nw - 1)) = wstar_i;
                    }
                    mat WS = W.t() * Siginv;
                    Siggam += ntk * (WS * W);
                    mugam += ntk * (WS * resid_i.t());
                }
                Siggam = 0.5 * (Siggam + Siggam.t());
                mat SiggamChol = arma::chol(Siggam);
                vec gtmp(nw * J);
                std::generate(gtmp.begin(), gtmp.end(), ::norm_rand);
                gamR.col(k) = arma::solve(arma::trimatu(SiggamChol), arma::solve(arma::trimatl(SiggamChol.t()), mugam) + gtmp);
                for (int j = 0; j < J; ++j)
                {
                    for (int j2 = 0; j2 < nn; ++j2)
                    {
                        gamR(nw * j + j2, k) = (1.0 - onstat(k)) * gamR(nw * j + j2, k);
                        gamR(nw * j + nn + j2, k) = onstat(k) * gamR(nw * j + nn + j2, k);
                    }
                }
            }

            // Update Omega
            Rcpp::Rcout << "Omega" << std::endl;
            for (int jj = 0; jj < J; ++jj)
            {
                mat gamstar = gamR.rows(nw * jj, nw * jj + nn - 1);
                mat qq1 = Omega0inv + (gamstar * gamstar.t());
                gamstar = gamR.rows(nw * jj + nn, nw * (jj + 1) - 1);
                mat qq2 = Omega0inv + (gamstar * gamstar.t());
                mat ominv1 = arma::wishrnd(arma::inv(qq1), shape_omega1);
                mat ominv2 = arma::wishrnd(arma::inv(qq2), shape_omega2);
                mat om1 = arma::inv(ominv1);
                mat om2 = arma::inv(ominv2);
                Omegainv(arma::span(nw * jj, nw * jj + nn - 1), arma::span(nw * jj, nw * jj + nn - 1)) = ominv1;
                Omegainv(arma::span(nw * jj + nn, nw * (jj + 1) - 1), arma::span(nw * jj + nn, nw * (jj + 1) - 1)) = ominv2;
                Omega(arma::span(nw * jj, nw * jj + nn - 1), arma::span(nw * jj, nw * jj + nn - 1)) = om1;
                Omega(arma::span(nw * jj + nn, nw * (jj + 1) - 1), arma::span(nw * jj + nn, nw * (jj + 1) - 1)) = om2;
            }

            // Update Sigma
            // Update sig2
            Rcpp::Rcout << "sig2" << std::endl;
            for (int i = 0; i < N; ++i)
            {
                int k = Trial(i);
                double ntk = Npt(i);
                rowvec w_i = WCovariate.row(i);
                rowvec wstar_i(nw, fill::zeros);
                if (Second(i) == 0)
                {
                    wstar_i.head(nn) = w_i;
                }
                else
                {
                    wstar_i.tail(nn) = w_i;
                }
                rowvec S_i = S.row(i);
                mat Stk = vechinv(S_i, J);
                vec gam_k = gamR.col(k);
                mat W(J, nw * J, fill::zeros);
                for (int j = 0; j < J; ++j)
                {
                    W(j, span(j * nw, (j + 1) * nw - 1)) = wstar_i;
                }
                vec resid_i = arma::trans(resid.row(i)) - W * gam_k;
                mat qq = ntk * resid_i * resid_i.t() + (ntk - 1.0) * Stk;
                rowvec delta_i = delta.row(i);
                rowvec r_i = Rhoinv.row(i);
                mat Rhoinv_tk = vecrinv(r_i, J);
                Rhoinv_tk.diag().fill(1.0);

                for (int j = 0; j < J; ++j)
                {
                    auto fx_delta = [&](double delta_input[]) -> double
                    {
                        return -loglik_delta_mcorr(delta_input[0], delta_i, j, Rhoinv_tk, qq, a0, b0, ntk);
                    };

                    double dstar = std::log(delta_i(j));
                    double start[] = {dstar};
                    double xmin[] = {0.0};
                    double ynewlo = 0.0;
                    double reqmin = arma::datum::eps;
                    int konvge = 5;
                    int kcount = 1000;
                    double step[] = {0.2};
                    int icount = 0;
                    int numres = 0;
                    int ifault = 0;
                    nelmin(fx_delta, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);
                    double xmax = xmin[0];
                    double minll = ynewlo;

                    mat cl(5, 3, fill::zeros);
                    vec dl(5, fill::zeros);
                    double step_size = delta_stepsize;

                    bool cont_flag = true;
                    while (cont_flag)
                    {
                        for (int iii = 0; iii < 5; ++iii)
                        {
                            double e1 = static_cast<double>(iii - 2);
                            cl(iii, 0) = std::pow(xmax + e1 * step_size, 2.0);
                            cl(iii, 1) = xmax + e1 * step_size;
                            cl(iii, 2) = 1.0;
                            dl(iii) = -loglik_delta_mcorr(xmax + e1 * step_size, delta_i, j, Rhoinv_tk, qq, a0, b0, ntk);
                        }
                        if (any(dl < minll))
                        {
                            step_size *= 1.2;
                        }
                        else
                        {
                            cont_flag = false;
                        }
                    }

                    vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
                    double sigmaa = std::sqrt(0.5 / fl(0));

                    // log-likelihood difference
                    double dstar_prop = ::norm_rand() * sigmaa + xmax;
                    double ll_diff = loglik_delta_mcorr(dstar_prop, delta_i, j, Rhoinv_tk, qq, a0, b0, ntk) - loglik_delta_mcorr(dstar, delta_i, j, Rhoinv_tk, qq, a0, b0, ntk) - 0.5 * (std::pow(dstar - xmax, 2.0) - std::pow(dstar_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
                    if (std::log(::unif_rand()) < ll_diff)
                    {
                        delta_i(j) = std::exp(dstar_prop);
                        delta(i, j) = delta_i(j);
                        mat siginvm = arma::diagmat(1.0 / delta_i);
                        mat Siginv_new = siginvm * Rhoinv_tk * siginvm;
                        Siginv_lt.row(i) = arma::trans(vech(Siginv_new));

                        ++delta_rates(i, j);
                    }
                }
            }

            // Update Rho
            Rcpp::Rcout << "Rho" << std::endl;
            for (int i = 0; i < N; ++i)
            {
                vec U(neta, fill::randn);
                vec etap = eta + SS * U;
                // log-likelihood difference
                double ll_diff = loglik_eta_mcorr(etap, WCovariate, ZCovariate, Trial, Second, gamR, delta, resid, Npt, S, phi_all, transform_type, I_JJm12, TOL) -
                                 loglik_eta_mcorr(eta, WCovariate, ZCovariate, Trial, Second, gamR, delta, resid, Npt, S, phi_all, transform_type, I_JJm12, TOL);
                bool eta_changed = false;
                if (std::log(::unif_rand()) < ll_diff)
                {
                    eta = etap;
                    ++eta_rates;
                    eta_changed = true;
                }
                else
                {
                    // delayed rejection
                    U.randn();
                    vec zzz = eta + std::sqrt(0.5) * (SS * U);

                    vec ystar = zzz - (etap - eta);
                    double log1pxy = std::log1p(-std::min(1.0, std::exp(ll_diff)));
                    double ll_diff_zystar = loglik_eta_mcorr(zzz, WCovariate, ZCovariate, Trial, Second, gamR, delta, resid, Npt, S, phi_all, transform_type, I_JJm12, TOL) - 
                                            loglik_eta_mcorr(ystar, WCovariate, ZCovariate, Trial, Second, gamR, delta, resid, Npt, S, phi_all, transform_type, I_JJm12, TOL);
                    double log1pzystar = std::log1p(-std::min(1.0, std::exp(ll_diff_zystar)));
                    double ll_diff_zx = loglik_eta_mcorr(zzz, WCovariate, ZCovariate, Trial, Second, gamR, delta, resid, Npt, S, phi_all, transform_type, I_JJm12, TOL) - 
                                        loglik_eta_mcorr(eta, WCovariate, ZCovariate, Trial, Second, gamR, delta, resid, Npt, S, phi_all, transform_type, I_JJm12, TOL);
                    ll_diff = ll_diff_zx + log1pzystar - log1pxy;
                    if (std::log(::unif_rand()) < ll_diff)
                    {
                        eta = zzz;
                        ++eta_rates;
                        eta_changed = true;
                    }
                }
                double alpha_n = std::min(1.0, std::exp(ll_diff));
                adapt_S(SS, U, alpha_n, alpha_star, icount_mh, gam_exp);
                
                // Update Sigmainvs
                if (eta_changed) {
                    for (int i = 0; i < N; ++i)
                    {
                        rowvec z_tk = ZCovariate.row(i);
                        mat Zmat = arma::kron(I_JJm12, z_tk);
                        vec xi_tk = Zmat * eta;
                        xi.row(i) = xi_tk.t();

                        // Calculate rho_tk from xi_tk
                        // Store the off-diagonal elements as a vector in vRho
                        mat rho_tk(J, J, fill::zeros);
                        if (transform_type == 1)
                        {
                            // rho is modeled as a matrix exponential of vecrinv(xi_tk)
                            mat varrho_tk = find_diag(xi_tk, TOL);
                            rho_tk = arma::expmat_sym(varrho_tk);
                            rho_tk.diag().fill(1.0);
                        }
                        else if (transform_type == 2)
                        {
                            // rho is modeled as recovered from the partial correlation matrix
                            // partial correlation matrix is further modeled as tanh(vecrinv(xi_tk))
                            vec tanh_xi_tk = arma::tanh(xi_tk);
                            mat pRho_tk = vecrinv(tanh_xi_tk, J);
                            rho_tk = pRho_to_Rho(pRho_tk);
                        }
                        else if (transform_type == 3)
                        {
                            // rho is modeled by its Cholesky factor
                            // The Cholesky factor is further modeled as hyperspherical coordinates (or hyperangles)
                            mat tmp = M_PI / (1.0 + arma::exp(-xi_tk));
                            mat varrho_tk = vecrinv(tmp, J);
                            mat B_tk = constructB(varrho_tk);
                            rho_tk = arma::trimatl(B_tk) * arma::trimatu(B_tk.t());
                        }

                        Rho.row(i) = arma::trans(vecr(rho_tk));
                        Rhoinv.row(i) = arma::trans(vecr(rho_tk.i()));

                        mat Sigma = arma::diagmat(delta.row(i)) * rho_tk * arma::diagmat(delta.row(i));
                        mat Sigmainv = arma::diagmat(1.0 / delta.row(i)) * rho_tk.i() * arma::diagmat(1.0 / delta.row(i));
                        Sig_lt.row(i) = arma::trans(vech(Sigma));
                        Siginv_lt.row(i) = arma::trans(vech(Sigmainv));
                    }
                }
            }

            // Update Rtk
            Rcpp::Rcout << "Rtk" << std::endl;
            for (int i = 0; i < N; ++i)
            {
                rowvec y_i = Outcome.row(i);
                rowvec vrtk = vRtk.row(i);
                mat V = arma::diagmat(SD.row(i));
                double ntk = Npt(i);
                mat Siginv = vechinv(arma::trans(Siginv_lt.row(i)), J);
                mat VSV = V * Siginv * V;
                for (int kk = 0; kk < J * (J - 1) / 2; ++kk)
                {
                    int iR = J - 2 - static_cast<int>(std::sqrt(-8.0 * static_cast<double>(kk) + 4.0 * static_cast<double>(J * (J - 1)) - 7.0) / 2.0 - 0.5); // row index
                    int iC = kk + iR + 1 - (J * (J - 1)) / 2 + ((J - iR) * ((J - iR) - 1)) / 2;                                                              // column index
                    auto fx_rstar = [&](double rstar_input[]) -> double
                    {
                        return -loglik_rik(rstar_input[0], vrtk, kk, J, iR, iC, ntk, VSV);
                    };

                    double zstar = vrtk(kk);
                    double start[] = {zstar};
                    double xmin[] = {0.0};
                    double ynewlo = 0.0;
                    double reqmin = arma::datum::eps;
                    int konvge = 5;
                    int kcount = 1000;
                    double step[] = {0.02};
                    int icount = 0;
                    int numres = 0;
                    int ifault = 0;
                    nelmin(fx_rstar, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);
                    double xmax = xmin[0];
                    double minll = ynewlo;

                    mat cl(5, 3, fill::zeros);
                    vec dl(5, fill::zeros);
                    double step_size = R_stepsize;

                    bool cont_flag = true;
                    while (cont_flag)
                    {
                        for (int iii = 0; iii < 5; ++iii)
                        {
                            double e1 = static_cast<double>(iii - 2);
                            cl(iii, 0) = std::pow(xmax + e1 * step_size, 2.0);
                            cl(iii, 1) = xmax + e1 * step_size;
                            cl(iii, 2) = 1.0;
                            dl(iii) = -loglik_rik(xmax + e1 * step_size, vrtk, kk, J, iR, iC, ntk, VSV);
                        }
                        if (any(dl < minll))
                        {
                            step_size *= 1.2;
                        }
                        else
                        {
                            cont_flag = false;
                        }
                    }

                    vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
                    double sigmaa = std::sqrt(0.5 / fl(0));

                    double zstar_prop = ::norm_rand() * sigmaa + xmax;
                    // log-likelihood difference
                    double ll_diff = loglik_rik(zstar_prop, vrtk, kk, J, iR, iC, ntk, VSV) - loglik_rik(zstar, vrtk, kk, J, iR, iC, ntk, VSV) - 0.5 * (std::pow(zstar - xmax, 2.0) - std::pow(zstar_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);

                    if (std::log(::unif_rand()) < ll_diff)
                    {
                        vrtk(kk) = zstar_prop;
                        vRtk(i, kk) = zstar_prop;
                        ++vR_rates(i, kk);
                    }
                }
                rowvec s_tk = SD.row(i);
                mat Vtk = arma::diagmat(s_tk);
                mat pRtk = vecrinv(arma::trans(arma::tanh(vRtk.row(i))), J);
                pRtk.diag().fill(1.0);
                mat Rtk = pRho_to_Rho(pRtk);
                S.row(i) = arma::trans(vech(Vtk * Rtk * Vtk));
            }

            // Update phis
            Rcpp::Rcout << "phis" << std::endl;
            phi_all /= phi_global; // remove phi_global from phi_all
            double b1_global = b1 + 0.5 * arma::accu(arma::square(eta) / phi_all);
            phi_global = b1_global / ::Rf_rgamma(a1 + 0.5 * static_cast<double>(neta), 1.0);
            phi_all *= phi_global; // update phi_all with new phi_global


            // Update phi_group
            for (int l = 0; l < JJm12; ++l) {
                
                double b2_group = b2 + 0.5 * arma::accu(arma::square(eta(arma::span(l*zcols, (l+1)*zcols-1))) / phi_local(arma::span(l*zcols, (l+1)*zcols-1))) / phi_global;

                phi_group(l) = b2_group / ::Rf_rgamma(a2_group, 1.0);
            }

            for (int l = 0; l < JJm12; ++l) {
                double divisor = phi_group(l) * phi_global;
                for (int i = l*zcols; i < (l+1)*zcols; ++i) {
                    phi_local(i) = (b3 + 0.5 * eta(i) * eta(i) / divisor) / ::Rf_rgamma(a3_local, 1.0);
                }
            }

            // Update phi_all
            phi_all = phi_local * phi_global;
            for (int i = 0; i < JJm12; ++i)
            {
                phi_all(arma::span(i * zcols, (i + 1) * zcols - 1)) *= phi_group(i);
            }

            prog.increment();
        }
    }

    if (verbose)
    {
        Rcout << "Sampling" << endl;
    }
    {
        Progress prog(nkeep, verbose);
        for (int ikeep = 0; ikeep < nkeep; ++ikeep) {
            if (Progress::check_abort()) {
				return Rcpp::List::create(Rcpp::Named("error") = "user interrupt aborted");
			}
            mat resid_ikeep(N, J, fill::zeros);
            for (int iskip = 0; iskip < nskip; ++iskip) {
                ++icount_mh;

                // Update theta
                mat Sig_theta(nt, nt, fill::zeros);
                Sig_theta.diag().fill(1.0 / c0);
                vec mu_theta(nt, fill::zeros);
                for (int k = 0; k < K; ++k)
                {
                    uvec idxk = idxks(k);
                    int n_k = idxk.n_elem;
                    mat XSX(nt, nt, fill::zeros);
                    mat WSX(nw * J, nt, fill::zeros);
                    vec WSy(nw * J, fill::zeros);
                    vec XSy(nt, fill::zeros);
                    mat Sig_gamk = Omegainv;
                    for (int i = 0; i < n_k; ++i)
                    {
                        int i_k = idxk(i);
                        rowvec x_i = XCovariate.row(i_k);
                        rowvec w_i = WCovariate.row(i_k);
                        rowvec wstar_i(nw, fill::zeros);
                        if (Second(i_k) == 0)
                        {
                            wstar_i.head(nn) = w_i;
                        }
                        else
                        {
                            wstar_i.tail(nn) = w_i;
                        }
                        rowvec y_i = Outcome.row(i_k);
                        double ntk = Npt(i_k);
                        mat Siginv = vechinv(arma::trans(Siginv_lt.row(i_k)), J);
                        mat X(J, xcols * J, fill::zeros);
                        mat W(J, nw * J, fill::zeros);
                        for (int j = 0; j < J; ++j)
                        {
                            X(j, arma::span(j * xcols, (j + 1) * xcols - 1)) = x_i;
                            W(j, arma::span(j * nw, (j + 1) * nw - 1)) = wstar_i;
                        }
                        mat Xstar = arma::join_horiz(X, W);
                        XSX += ntk * (Xstar.t() * Siginv * Xstar);
                        XSy += ntk * (Xstar.t() * Siginv * y_i.t());
                        mat WS = ntk * (W.t() * Siginv);
                        Sig_gamk += WS * W;
                        WSX += WS * Xstar;
                        WSy += WS * y_i.t();
                    }
                    mat Sig_gamk_inv = arma::inv(Sig_gamk);
                    Sig_theta += XSX - WSX.t() * Sig_gamk_inv * WSX;
                    mu_theta += XSy - WSX.t() * Sig_gamk_inv * WSy;
                }
                Sig_theta = 0.5 * (Sig_theta + Sig_theta.t());
                mat Sig_theta_chol = arma::chol(Sig_theta);
                vec ttmp(nt);
                std::generate(ttmp.begin(), ttmp.end(), ::norm_rand);
                theta = arma::solve(arma::trimatu(Sig_theta_chol), arma::solve(arma::trimatl(Sig_theta_chol.t()), mu_theta) + ttmp);

                for (int i = 0; i < N; ++i)
                {
                    rowvec x_i = XCovariate.row(i);
                    rowvec w_i = WCovariate.row(i);
                    rowvec y_i = Outcome.row(i);
                    rowvec wstar_i(nw, fill::zeros);
                    if (Second(i) == 0)
                    {
                        wstar_i.head(nn) = w_i;
                    }
                    else
                    {
                        wstar_i.tail(nn) = w_i;
                    }
                    mat X(J, xcols * J, fill::zeros);
                    mat W(J, nw * J, fill::zeros);
                    for (int j = 0; j < J; ++j)
                    {
                        X(j, arma::span(j * xcols, (j + 1) * xcols - 1)) = x_i;
                        W(j, arma::span(j * nw, (j + 1) * nw - 1)) = wstar_i;
                    }
                    mat Xstar = arma::join_horiz(X, W);
                    resid.row(i) = arma::trans(y_i.t() - Xstar * theta);
                }

                // Update gamR
                for (int k = 0; k < K; ++k)
                {
                    mat Siggam = Omegainv;
                    vec mugam(nw * J, fill::zeros);
                    uvec idxk = idxks(k);
                    int n_k = idxk.n_elem;
                    for (int i = 0; i < n_k; ++i)
                    {
                        int i_k = idxk(i);
                        rowvec w_i = WCovariate.row(i_k);
                        rowvec wstar_i(nw, fill::zeros);
                        if (Second(i_k) == 0)
                        {
                            wstar_i.head(nn) = w_i;
                        }
                        else
                        {
                            wstar_i.tail(nn) = w_i;
                        }
                        rowvec resid_i = resid.row(i_k);
                        double ntk = Npt(i_k);
                        mat Siginv = vechinv(arma::trans(Siginv_lt.row(i_k)), J);
                        mat W(J, nw * J, fill::zeros);
                        for (int j = 0; j < J; ++j)
                        {
                            W(j, arma::span(j * nw, (j + 1) * nw - 1)) = wstar_i;
                        }
                        mat WS = W.t() * Siginv;
                        Siggam += ntk * (WS * W);
                        mugam += ntk * (WS * resid_i.t());
                    }
                    Siggam = 0.5 * (Siggam + Siggam.t());
                    mat SiggamChol = arma::chol(Siggam);
                    vec gtmp(nw * J);
                    std::generate(gtmp.begin(), gtmp.end(), ::norm_rand);
                    gamR.col(k) = arma::solve(arma::trimatu(SiggamChol), arma::solve(arma::trimatl(SiggamChol.t()), mugam) + gtmp);
                    for (int j = 0; j < J; ++j)
                    {
                        for (int j2 = 0; j2 < nn; ++j2)
                        {
                            gamR(nw * j + j2, k) = (1.0 - onstat(k)) * gamR(nw * j + j2, k);
                            gamR(nw * j + nn + j2, k) = onstat(k) * gamR(nw * j + nn + j2, k);
                        }
                    }
                }

                // Update Omega
                for (int jj = 0; jj < J; ++jj)
                {
                    mat gamstar = gamR.rows(nw * jj, nw * jj + nn - 1);
                    mat qq1 = Omega0inv + (gamstar * gamstar.t());
                    gamstar = gamR.rows(nw * jj + nn, nw * (jj + 1) - 1);
                    mat qq2 = Omega0inv + (gamstar * gamstar.t());
                    mat ominv1 = arma::wishrnd(arma::inv(qq1), shape_omega1);
                    mat ominv2 = arma::wishrnd(arma::inv(qq2), shape_omega2);
                    mat om1 = arma::inv(ominv1);
                    mat om2 = arma::inv(ominv2);
                    Omegainv(arma::span(nw * jj, nw * jj + nn - 1), arma::span(nw * jj, nw * jj + nn - 1)) = ominv1;
                    Omegainv(arma::span(nw * jj + nn, nw * (jj + 1) - 1), arma::span(nw * jj + nn, nw * (jj + 1) - 1)) = ominv2;
                    Omega(arma::span(nw * jj, nw * jj + nn - 1), arma::span(nw * jj, nw * jj + nn - 1)) = om1;
                    Omega(arma::span(nw * jj + nn, nw * (jj + 1) - 1), arma::span(nw * jj + nn, nw * (jj + 1) - 1)) = om2;
                }

                // Update Sigma
                // Update sig2
                for (int i = 0; i < N; ++i)
                {
                    int k = Trial(i);
                    double ntk = Npt(i);
                    rowvec w_i = WCovariate.row(i);
                    rowvec wstar_i(nw, fill::zeros);
                    if (Second(i) == 0)
                    {
                        wstar_i.head(nn) = w_i;
                    }
                    else
                    {
                        wstar_i.tail(nn) = w_i;
                    }
                    rowvec S_i = S.row(i);
                    mat Stk = vechinv(S_i, J);
                    vec gam_k = gamR.col(k);
                    mat W(J, nw * J, fill::zeros);
                    for (int j = 0; j < J; ++j)
                    {
                        W(j, span(j * nw, (j + 1) * nw - 1)) = wstar_i;
                    }
                    vec resid_i = arma::trans(resid.row(i)) - W * gam_k;
                    resid_ikeep.row(i) = resid_i.t();
                    mat qq = ntk * resid_i * resid_i.t() + (ntk - 1.0) * Stk;
                    rowvec delta_i = delta.row(i);
                    rowvec r_i = Rhoinv.row(i);
                    mat Rhoinv_tk = vecrinv(r_i, J);
                    Rhoinv_tk.diag().fill(1.0);

                    for (int j = 0; j < J; ++j)
                    {
                        auto fx_delta = [&](double delta_input[]) -> double
                        {
                            return -loglik_delta_mcorr(delta_input[0], delta_i, j, Rhoinv_tk, qq, a0, b0, ntk);
                        };

                        double dstar = std::log(delta_i(j));
                        double start[] = {dstar};
                        double xmin[] = {0.0};
                        double ynewlo = 0.0;
                        double reqmin = arma::datum::eps;
                        int konvge = 5;
                        int kcount = 1000;
                        double step[] = {0.2};
                        int icount = 0;
                        int numres = 0;
                        int ifault = 0;
                        nelmin(fx_delta, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);
                        double xmax = xmin[0];
                        double minll = ynewlo;

                        mat cl(5, 3, fill::zeros);
                        vec dl(5, fill::zeros);
                        double step_size = delta_stepsize;

                        bool cont_flag = true;
                        while (cont_flag)
                        {
                            for (int iii = 0; iii < 5; ++iii)
                            {
                                double e1 = static_cast<double>(iii - 2);
                                cl(iii, 0) = std::pow(xmax + e1 * step_size, 2.0);
                                cl(iii, 1) = xmax + e1 * step_size;
                                cl(iii, 2) = 1.0;
                                dl(iii) = -loglik_delta_mcorr(xmax + e1 * step_size, delta_i, j, Rhoinv_tk, qq, a0, b0, ntk);
                            }
                            if (any(dl < minll))
                            {
                                step_size *= 1.2;
                            }
                            else
                            {
                                cont_flag = false;
                            }
                        }

                        vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
                        double sigmaa = std::sqrt(0.5 / fl(0));

                        // log-likelihood difference
                        double dstar_prop = ::norm_rand() * sigmaa + xmax;
                        double ll_diff = loglik_delta_mcorr(dstar_prop, delta_i, j, Rhoinv_tk, qq, a0, b0, ntk) - loglik_delta_mcorr(dstar, delta_i, j, Rhoinv_tk, qq, a0, b0, ntk) - 0.5 * (std::pow(dstar - xmax, 2.0) - std::pow(dstar_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
                        if (std::log(::unif_rand()) < ll_diff)
                        {
                            delta_i(j) = std::exp(dstar_prop);
                            delta(i, j) = delta_i(j);
                            mat siginvm = arma::diagmat(1.0 / delta_i);
                            mat Siginv_new = siginvm * Rhoinv_tk * siginvm;
                            Siginv_lt.row(i) = arma::trans(vech(Siginv_new));

                            ++delta_rates(i, j);
                        }
                    }
                }

                // Update Rho
                for (int i = 0; i < N; ++i)
                {
                    vec U(neta, fill::randn);
                    vec etap = eta + SS * U;
                    // log-likelihood difference
                    double ll_diff = loglik_eta_mcorr(etap, WCovariate, ZCovariate, Trial, Second, gamR, delta, resid, Npt, S, phi_all, transform_type, I_JJm12, TOL) -
                                     loglik_eta_mcorr(eta, WCovariate, ZCovariate, Trial, Second, gamR, delta, resid, Npt, S, phi_all, transform_type, I_JJm12, TOL);
                    bool eta_changed = false;
                    if (std::log(::unif_rand()) < ll_diff)
                    {
                        eta = etap;
                        ++eta_rates;
                        eta_changed = true;
                    }
                    else
                    {
                        // delayed rejection
                        U.randn();
                        vec zzz = eta + std::sqrt(0.5) * (SS * U);

                        vec ystar = zzz - (etap - eta);
                        double log1pxy = std::log1p(-std::min(1.0, std::exp(ll_diff)));
                        double ll_diff_zystar = loglik_eta_mcorr(zzz, WCovariate, ZCovariate, Trial, Second, gamR, delta, resid, Npt, S, phi_all, transform_type, I_JJm12, TOL) -
                                                loglik_eta_mcorr(ystar, WCovariate, ZCovariate, Trial, Second, gamR, delta, resid, Npt, S, phi_all, transform_type, I_JJm12, TOL);
                        double log1pzystar = std::log1p(-std::min(1.0, std::exp(ll_diff_zystar)));
                        double ll_diff_zx = loglik_eta_mcorr(zzz, WCovariate, ZCovariate, Trial, Second, gamR, delta, resid, Npt, S, phi_all, transform_type, I_JJm12, TOL) -
                                            loglik_eta_mcorr(eta, WCovariate, ZCovariate, Trial, Second, gamR, delta, resid, Npt, S, phi_all, transform_type, I_JJm12, TOL);
                        ll_diff = ll_diff_zx + log1pzystar - log1pxy;
                        if (std::log(::unif_rand()) < ll_diff)
                        {
                            eta = zzz;
                            ++eta_rates;
                            eta_changed = true;
                        }
                    }
                    double alpha_n = std::min(1.0, std::exp(ll_diff));
                    adapt_S(SS, U, alpha_n, alpha_star, icount_mh, gam_exp);

                    // Update Sigmainvs
                    if (eta_changed)
                    {
                        for (int i = 0; i < N; ++i)
                        {
                            rowvec z_tk = ZCovariate.row(i);
                            mat Zmat = arma::kron(I_JJm12, z_tk);
                            vec xi_tk = Zmat * eta;
                            xi.row(i) = xi_tk.t();

                            // Calculate rho_tk from xi_tk
                            // Store the off-diagonal elements as a vector in vRho
                            mat rho_tk(J, J, fill::zeros);
                            if (transform_type == 1)
                            {
                                // rho is modeled as a matrix exponential of vecrinv(xi_tk)
                                mat varrho_tk = find_diag(xi_tk, TOL);
                                rho_tk = arma::expmat_sym(varrho_tk);
                                rho_tk.diag().fill(1.0);
                            }
                            else if (transform_type == 2)
                            {
                                // rho is modeled as recovered from the partial correlation matrix
                                // partial correlation matrix is further modeled as tanh(vecrinv(xi_tk))
                                vec tanh_xi_tk = arma::tanh(xi_tk);
                                mat pRho_tk = vecrinv(tanh_xi_tk, J);
                                rho_tk = pRho_to_Rho(pRho_tk);
                            }
                            else if (transform_type == 3)
                            {
                                // rho is modeled by its Cholesky factor
                                // The Cholesky factor is further modeled as hyperspherical coordinates (or hyperangles)
                                mat tmp = M_PI / (1.0 + arma::exp(-xi_tk));
                                mat varrho_tk = vecrinv(tmp, J);
                                mat B_tk = constructB(varrho_tk);
                                rho_tk = arma::trimatl(B_tk) * arma::trimatu(B_tk.t());
                            }

                            Rho.row(i) = arma::trans(vecr(rho_tk));
                            Rhoinv.row(i) = arma::trans(vecr(rho_tk.i()));

                            mat Sigma = arma::diagmat(delta.row(i)) * rho_tk * arma::diagmat(delta.row(i));
                            mat Sigmainv = arma::diagmat(1.0 / delta.row(i)) * rho_tk.i() * arma::diagmat(1.0 / delta.row(i));
                            Sig_lt.row(i) = arma::trans(vech(Sigma));
                            Siginv_lt.row(i) = arma::trans(vech(Sigmainv));
                        }
                    }
                }

                // Update Rtk
                for (int i = 0; i < N; ++i)
                {
                    rowvec y_i = Outcome.row(i);
                    rowvec vrtk = vRtk.row(i);
                    mat V = arma::diagmat(SD.row(i));
                    double ntk = Npt(i);
                    mat Siginv = vechinv(arma::trans(Siginv_lt.row(i)), J);
                    mat VSV = V * Siginv * V;
                    for (int kk = 0; kk < J * (J - 1) / 2; ++kk)
                    {
                        int iR = J - 2 - static_cast<int>(std::sqrt(-8.0 * static_cast<double>(kk) + 4.0 * static_cast<double>(J * (J - 1)) - 7.0) / 2.0 - 0.5); // row index
                        int iC = kk + iR + 1 - (J * (J - 1)) / 2 + ((J - iR) * ((J - iR) - 1)) / 2;                                                              // column index
                        auto fx_rstar = [&](double rstar_input[]) -> double
                        {
                            return -loglik_rik(rstar_input[0], vrtk, kk, J, iR, iC, ntk, VSV);
                        };

                        double zstar = vrtk(kk);
                        double start[] = {zstar};
                        double xmin[] = {0.0};
                        double ynewlo = 0.0;
                        double reqmin = arma::datum::eps;
                        int konvge = 5;
                        int kcount = 1000;
                        double step[] = {0.02};
                        int icount = 0;
                        int numres = 0;
                        int ifault = 0;
                        nelmin(fx_rstar, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);
                        double xmax = xmin[0];
                        double minll = ynewlo;

                        mat cl(5, 3, fill::zeros);
                        vec dl(5, fill::zeros);
                        double step_size = R_stepsize;

                        bool cont_flag = true;
                        while (cont_flag)
                        {
                            for (int iii = 0; iii < 5; ++iii)
                            {
                                double e1 = static_cast<double>(iii - 2);
                                cl(iii, 0) = std::pow(xmax + e1 * step_size, 2.0);
                                cl(iii, 1) = xmax + e1 * step_size;
                                cl(iii, 2) = 1.0;
                                dl(iii) = -loglik_rik(xmax + e1 * step_size, vrtk, kk, J, iR, iC, ntk, VSV);
                            }
                            if (any(dl < minll))
                            {
                                step_size *= 1.2;
                            }
                            else
                            {
                                cont_flag = false;
                            }
                        }

                        vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
                        double sigmaa = std::sqrt(0.5 / fl(0));

                        double zstar_prop = ::norm_rand() * sigmaa + xmax;
                        // log-likelihood difference
                        double ll_diff = loglik_rik(zstar_prop, vrtk, kk, J, iR, iC, ntk, VSV) - loglik_rik(zstar, vrtk, kk, J, iR, iC, ntk, VSV) - 0.5 * (std::pow(zstar - xmax, 2.0) - std::pow(zstar_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);

                        if (std::log(::unif_rand()) < ll_diff)
                        {
                            vrtk(kk) = zstar_prop;
                            vRtk(i, kk) = zstar_prop;
                            ++vR_rates(i, kk);
                        }
                    }
                    rowvec s_tk = SD.row(i);
                    mat Vtk = arma::diagmat(s_tk);
                    mat pRtk = vecrinv(arma::trans(arma::tanh(vRtk.row(i))), J);
                    pRtk.diag().fill(1.0);
                    mat Rtk = pRho_to_Rho(pRtk);
                    S.row(i) = arma::trans(vech(Vtk * Rtk * Vtk));
                }

                // Update phis
                phi_all /= phi_global; // remove phi_global from phi_all
                double b1_global = b1 + 0.5 * arma::accu(arma::square(eta) / phi_all);
                phi_global = b1_global / ::Rf_rgamma(a1 + 0.5 * static_cast<double>(neta), 1.0);
                phi_all *= phi_global; // update phi_all with new phi_global

                // Update phi_group
                for (int l = 0; l < JJm12; ++l)
                {

                    double b2_group = b2 + 0.5 * arma::accu(arma::square(eta(arma::span(l * zcols, (l + 1) * zcols - 1))) / phi_local(arma::span(l * zcols, (l + 1) * zcols - 1))) / phi_global;

                    phi_group(l) = b2_group / ::Rf_rgamma(a2_group, 1.0);
                }

                for (int l = 0; l < JJm12; ++l)
                {
                    double divisor = phi_group(l) * phi_global;
                    for (int i = l * zcols; i < (l + 1) * zcols; ++i)
                    {
                        phi_local(i) = (b3 + 0.5 * eta(i) * eta(i) / divisor) / ::Rf_rgamma(a3_local, 1.0);
                    }
                }

                // Update phi_all
                phi_all = phi_local * phi_global;
                for (int i = 0; i < JJm12; ++i)
                {
                    phi_all(arma::span(i * zcols, (i + 1) * zcols - 1)) *= phi_group(i);
                }
            }

            theta_save.col(ikeep) = theta;
            mat Rtk(arma::size(vRtk), fill::zeros);
			for (int i = 0; i < N; ++i) {
				mat RR = vecrinv(trans(arma::tanh(vRtk.row(i))), J);
				RR.diag().fill(1.0);
				mat R = pRho_to_Rho(RR);
				Rtk.row(i) = arma::trans(vecr(R));
			}
			Rtk_save.slice(ikeep) = Rtk;
            resid_save.slice(ikeep) = resid_ikeep;
            delta_save.slice(ikeep) = delta;
            eta_save.col(ikeep) = eta;

            prog.increment();
        }
    }
    return Rcpp::List::create(
			Rcpp::Named("resid") = resid_save,
			Rcpp::Named("theta") = theta_save,
			Rcpp::Named("delta") = delta_save,
			Rcpp::Named("eta") = eta_save,
			Rcpp::Named("R") = Rtk_save,
			Rcpp::Named("delta_acceptance") = delta_rates / static_cast<double>(ndiscard + nskip*nkeep),
			Rcpp::Named("eta_acceptance") = eta_rates / static_cast<double>(ndiscard + nskip*nkeep),
		    Rcpp::Named("vR_acceptance") = vR_rates / static_cast<double>(ndiscard + nskip*nkeep)
		);
}
