#' Print results
#' 
#' @param object the output model from fitting a meta analysis/regression model
#' @return does not return anything; print a summary of the output
#' @importFrom coda mcmc HPDinterval
#' @method print bayes.parobs
#' @export
"print.bayes.parobs" <- function(object, conf.level=0.95, HPD=TRUE, ...) {
	cat("Bayesian Inference for Multivariate Meta-Regression\nWith a Partially Observed Within-Study Sample Covariance Matrix\n")
	cat("\n")
	cat("Model:\n")
	cat("  (Aggregate mean)\n    y_tk = X_tk * beta + W_tk * gamma + N(0, Sigma_kt / n_kt)\n")
	cat("  (Sample Variance)\n    (n_kt - 1) S_tk ~ Wishart(n_kt - 1, Sigma_tk)\n")
	cat("  (Random effects)\n    ")
	cat("[gamma_k | Omega] ~ N(0, Omega)\n")
	cat("Priors:\n")
	cat("   beta ~ MVN(0, c0 * I_p), c0=", object$prior$c0, "\n")
	cat("   Omega_j^{-1} ~ Wishart(dj0, Omega0), dj0=", object$prior$dj0, "\n")
	if (object$fmodel == 1) {
		cat("   Sigma_tk = diag(sig_{tk,11}^2, ..., sig_{tk,JJ}^2)\n")
		cat("   where sig_{tk,jj}^2 ~ IG(s0, d0), s0=", object$prior$s0, ", d0=", object$prior$d0, "\n")
	} else if (object$fmodel == 2) {
		cat("   Sigma_tk = Sigma, where Sigma ~ Wishart(s0, Sigma0)\n")
		cat("   s0=", object$prior$s0, "\n")
	} else if (object$fmodel == 3) {
		cat("   Sigma_{tk} = sig_{tk} * Rho * sig_{tk},\n")
		cat("   where p(Rho) = 1, and sig_{tk,jj} ~ IG(s0, d0)\n")
	} else if (object$fmodel == 4) {
		cat("   Sigma_{tk}^{-1} ~ Wishart(nu0, (nu0-J-1)*Sigma),\n")
		cat("   Sigma ~ Wishart(d0, Sigma0),\n")
		cat("   where nu0=", object$prior$nu0, ", d0=", object$prior$d0, "\n")
	}

	cat("-------------------------------------------------\n")

	mm <- matrix("", 2, 2)
	mm[1,1] <- "Number of trials"
	mm[1,2] <- as.character(object$K)
	mm[2,1] <- "Number of treatments"
	mm[2,2] <- as.character(object$T)
	print(mm, justify="left")

	theta <- list()
	theta$mean <- rowMeans(object$mcmc.draws$theta)
	theta$sd <- apply(object$mcmc.draws$theta, 1, sd)

	sig.level <- 1 - conf.level

	if (HPD) {
		theta.hpd <- coda::HPDinterval(mcmc(t(object$mcmc.draws$theta), end=object$mcmc$nkeep), prob=conf.level)
		theta$lower <- theta.hpd[,1]
		theta$upper <- theta.hpd[,2]
	} else {
		theta$lower <- apply(object$mcmc.draws$theta, 1, function(xx) quantile(xx, prob = sig.level/2))
		theta$upper <- apply(object$mcmc.draws$theta, 1, function(xx) quantile(xx, prob = 1-sig.level/2))
	}
	theta_print <- cbind(theta$mean, theta$sd, theta$lower, theta$upper)
	# p_print <- rbind(beta_print, phi_print, gam_print)
	if (HPD) {
		colnames(theta_print) <- c("Est.", "Std.Dev", "HPD(Lower)", "HPD(Upper)")
	} else {
		colnames(theta_print) <- c("Est.", "Std.Dev", "CI(Lower)", "CI(Upper)")
	}
	rownames(theta_print) <- c(paste0("theta", 1:length(theta$mean)))
	# rownames(p_print) <- c(paste0("beta", 1:length(beta$mean)),
	# 					   paste0("phi", 1:length(phi$mean)),
	# 					   paste0("gam", 1:length(gam$mean)))
	print(theta_print, justify="left", digits=2)
	cat("-----------------------------\n")
	cat("*Credible level: ", conf.level, "\n")
}
