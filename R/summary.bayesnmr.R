#' Summarize results
#' 
#' @param object the output model from fitting a network meta analysis/regression model
#' @param ... additional arguments for print
#' @return does not return anything; print a summary of the output
#' @importFrom coda mcmc HPDinterval
#' @importFrom stats quantile sd
#' @method summary bayesnmr
#' @export
"summary.bayesnmr" <- function(object, ...) {
	digits <- max(3, getOption("digits") - 3)
	theta <- list()
	phi <- list()
	gam <- list()
	sig2 <- list()
	Rho <- list()
	param <- object$mcmc.draws$theta
	if (object$scale_x) {
		xcols <- ncol(object$Covariate)
		tlength <- nrow(param)
		trlength <- tlength - xcols
		tscale <- c(apply(object$Covariate, 2, sd), rep(1, trlength))
	} else {
		tlength <- nrow(param)
		tscale <- rep(1, tlength)
	}
	theta.post <- vapply(1:object$mcmc$nkeep, function(ikeep) {
		param[,ikeep] / tscale
	}, FUN.VALUE = numeric(tlength))
	theta$mean <- rowMeans(param)
	theta$sd <- apply(param, 1, sd)
	phi$mean <- rowMeans(object$mcmc.draws$phi)
	phi$sd <- apply(object$mcmc.draws$phi, 1, sd)
	gam$mean <- rowMeans(object$mcmc.draws$gam)
	gam$sd <- apply(object$mcmc.draws$gam, 1, sd)

	conf.level <- 0.95
	sig.level <- 1 - conf.level

	if (HPD) {
		theta.hpd <- coda::HPDinterval(mcmc(t(param), end=object$mcmc$nkeep), prob=conf.level)
		theta$lower <- theta.hpd[,1]
		theta$upper <- theta.hpd[,2]
	} else {
		theta$lower <- apply(param, 1, function(xx) quantile(xx, prob = sig.level/2))
		theta$upper <- apply(param, 1, function(xx) quantile(xx, prob = 1-sig.level/2))
	}
	r <- cbind(theta$mean, theta$sd, theta$lower, theta$upper)
	if (HPD) {
		colnames(r) <- c("Post.Mean", "Std.Dev", "HPD(Lower)", "HPD(Upper)")
	} else {
		colnames(r) <- c("Post.Mean", "Std.Dev", "CI(Lower)", "CI(Upper)")
	}
	cat("\nPosterior inference in network meta-regression models\n")
	cat("Fixed-effects:\n")
	r <- round(r, digits=digits)
	print.default(r, print.gap = 2)
	invisible(r)
}
