#' Print results
#' 
#' @param object the output model from fitting a meta analysis/regression model
#' @return does not return anything; print a summary of the output
#' @importFrom coda mcmc HPDinterval
#' @method print bayesnmr
#' @export
"print.bayesnmr" <- function(object, conf.level=0.95, HPD=TRUE, ...) {
	cat("Bayesian Network Meta-Regression Hierarchical Models\nUsing Heavy-Tailed Multivariate Random Effects\nwith Covariate-Dependent Variances\n")
	cat("\n")
	cat("Model:\n")
	cat("  (Aggregate mean)\n    y_kt = x_kt'beta + tau_kt * gamma_kt + N(0, sigma_kt^2 / n_kt)\n")
	cat("  (Sample Variance)\n    (n_kt - 1) S^2 / sigma_kt^2 ~ chi^2(n_kt - 1)\n")
	cat("  (Random effects)\n    ")
	cat("[gam | Rho,nu] ~ MVT(0, E_k' Rho E_k, nu)\n")
	cat("Priors:\n")
	cat("  beta       ~ MVN(0, c01 * I_p), c01=", object$prior$c01, "\n")
	cat("  phi        ~ MVN(0, c02 * I_q), c02=", object$prior$c02, "\n")
	cat("  p(sigma^2) ~ 1/sigma^2 * I(sigma^2 > 0)\n")
	cat("  p(Rho)     ~ 1\n")

	cat("-------------------------------------------------\n")

	cat("Number of studies:    ", object$K, "\n")
	cat("Number of arms:       ", length(object$y), "\n")
	cat("Number of treatments: ", object$nT, "\n")
	beta <- list()
	phi <- list()
	gam <- list()
	sig2 <- list()
	Rho <- list()
	beta$mean <- rowMeans(object$mcmc.draws$beta)
	beta$sd <- apply(object$mcmc.draws$beta, 1, sd)
	phi$mean <- rowMeans(object$mcmc.draws$phi)
	phi$sd <- apply(object$mcmc.draws$phi, 1, sd)
	gam$mean <- rowMeans(object$mcmc.draws$gam)
	gam$sd <- apply(object$mcmc.draws$gam, 1, sd)

	sig.level <- 1 - conf.level

	if (HPD) {
		beta.hpd <- coda::HPDinterval(mcmc(t(object$mcmc.draws$beta), end=object$mcmc$nkeep), prob=conf.level)
		beta$lower <- beta.hpd[,1]
		beta$upper <- beta.hpd[,2]

		phi.hpd <- coda::HPDinterval(mcmc(t(object$mcmc.draws$phi), end=object$mcmc$nkeep), prob=conf.level)
		phi$lower <- phi.hpd[,1]
		phi$upper <- phi.hpd[,2]

		gam.hpd <- coda::HPDinterval(mcmc(t(object$mcmc.draws$gam), end=object$mcmc$nkeep), prob=conf.level)
		gam$lower <- gam.hpd[,1]
		gam$upper <- gam.hpd[,2]
	} else {
		beta$lower <- apply(object$mcmc.draws$beta, 1, function(xx) quantile(xx, prob = sig.level/2))
		beta$upper <- apply(object$mcmc.draws$beta, 1, function(xx) quantile(xx, prob = 1-sig.level/2))

		phi$lower <- apply(object$mcmc.draws$phi, 1, function(xx) quantile(xx, prob = sig.level/2))
		phi$upper <- apply(object$mcmc.draws$phi, 1, function(xx) quantile(xx, prob = 1-sig.level/2))

		gam$lower <- apply(object$mcmc.draws$gam, 1, function(xx) quantile(xx, prob = sig.level/2))
		gam$upper <- apply(object$mcmc.draws$gam, 1, function(xx) quantile(xx, prob = 1-sig.level/2))
	}
	beta_print <- cbind(beta$mean, beta$sd, beta$lower, beta$upper)
	phi_print <- cbind(phi$mean, phi$sd, phi$lower, phi$upper)
	gam_print <- cbind(gam$mean, gam$sd, gam$lower, gam$upper)
	p_print <- rbind(beta_print, phi_print, gam_print)
	if (HPD) {
		colnames(p_print) <- c("Est.", "Std.Dev", "HPD(Lower)", "HPD(Upper)")
	} else {
		colnames(p_print) <- c("Est.", "Std.Dev", "CI(Lower)", "CI(Upper)")
	}
	rownames(p_print) <- c(paste0("beta", 1:length(beta$mean)),
						   paste0("phi", 1:length(phi$mean)),
						   paste0("gam", 1:length(gam$mean)))
	print(p_print, justify="left", digits=2)
	cat("-----------------------------\n")
	cat("*Credible level: ", conf.level, "\n")
}
