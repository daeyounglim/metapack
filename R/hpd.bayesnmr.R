#' get the highest posterior density (HPD) interval
#' @param object the output model from fitting a (network) meta analysis/regression model
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the probability which the HPD interval will cover
#' @param HPD a logical value indicating whether HPD or equal-tailed credible interval should be computed; by default, TRUE
#' @importFrom coda mcmc HPDinterval
#' @return dataframe containing HPD intervals for the parameters
#' @method hpd bayesnmr
#' @export

"hpd.bayesnmr" <- function(object, parm, level = 0.95, HPD = TRUE) {
	if (object$scale_x) {
		xcols <- ncol(object$XCovariate)
		tlength <- nrow(object$mcmc.draws$theta)
		trlength <- tlength - xcols
		tscale <- c(unname(attr(object$XCovariate, "scaled:scale")), rep(1, trlength))
	} else {
		tlength <- nrow(object$mcmc.draws$theta)
		tscale <- rep(1, tlength)
	}
	theta.post <- vapply(1:object$mcmc$nkeep, function(ikeep) {
		object$mcmc.draws$theta[,ikeep] / tscale
	}, FUN.VALUE = numeric(tlength))
	if (missing(parm)) {
		out <- list()
		if (HPD) {
			out$theta <- coda::HPDinterval(coda::mcmc(t(theta.post), end=object$mcmc$nkeep), prob=level)
			out$phi <- coda::HPDinterval(coda::mcmc(t(object$mcmc.draws$phi), end=object$mcmc$nkeep), prob=level)
			out$gam <- coda::HPDinterval(coda::mcmc(t(object$mcmc.draws$gam), end=object$mcmc$nkeep), prob=level)
			out$sig2 <- coda::HPDinterval(coda::mcmc(t(object$mcmc.draws$sig2), end=object$mcmc$nkeep), prob=level)
			out$Rho <- hpdarray(object$mcmc.draws$Rho, level = level)
			if (object$control$sample_df) {
				out$df <- coda::HPDinterval(coda::mcmc(object$mcmc.draws$df, end=object$mcmc$nkeep), prob=level)
			}
			attr(out, "type") <- "HPD"
		} else {
			sig.level <- 1 - level
			out$theta <- apply(theta.post, 1, function(xx) quantile(xx, prob = c(sig.level/2, 1 - sig.level/2)))
			out$phi <- apply(object$mcmc.draws$phi, 1, function(xx) quantile(xx, prob = c(sig.level/2, 1 - sig.level/2)))
			out$gam <- apply(object$mcmc.draws$gam, 1, function(xx) quantile(xx, prob = c(sig.level/2, 1 - sig.level/2)))
			out$sig2 <- apply(object$mcmc.draws$sig2, 1, function(xx) quantile(xx, prob = c(sig.level/2, 1 - sig.level/2)))
			out$Rho <- ciarray(object$mcmc.draws$Rho, level = level)
			if (object$control$sample_df) {
				out$df <- quantile(object$mcmc.draws$df, c(sig.level/2, 1 - sig.level/2))
			}
			attr(out, "type") <- "equal-tailed CI"
		}
		class(out) <- "bayesnmr.hpd"
		return(out)
	} else {
		cl <- list()
		if (HPD) {
			cl$theta <-quote(coda::HPDinterval(coda::mcmc(t(object$mcmc.draws$theta), end=object$mcmc$nkeep), prob=level))
			cl$phi <- quote(coda::HPDinterval(coda::mcmc(t(object$mcmc.draws$phi), end=object$mcmc$nkeep), prob=level))
			cl$gam <- quote(coda::HPDinterval(coda::mcmc(t(object$mcmc.draws$gam), end=object$mcmc$nkeep), prob=level))
			cl$sig2 <- quote(coda::HPDinterval(coda::mcmc(t(object$mcmc.draws$sig2), end=object$mcmc$nkeep), prob=level))
			cl$Rho <- quote(hpdarray(object$mcmc.draws$Rho, level = level))
			if (object$control$sample_df) {
				cl$df <- coda::HPDinterval(coda::mcmc(object$mcmc.draws$df, end=object$mcmc$nkeep), prob=level)
			}
		} else {
			sig.level <- 1 - level
			cl$theta <- quote(apply(theta.post, 1, function(xx) quantile(xx, prob = c(sig.level/2, 1 - sig.level/2))))
			cl$phi <- quote(apply(object$mcmc.draws$phi, 1, function(xx) quantile(xx, prob = c(sig.level/2, 1 - sig.level/2))))
			cl$gam <- quote(apply(object$mcmc.draws$gam, 1, function(xx) quantile(xx, prob = c(sig.level/2, 1 - sig.level/2))))
			cl$sig2 <- quote(apply(object$mcmc.draws$sig2, 1, function(xx) quantile(xx, prob = c(sig.level/2, 1 - sig.level/2))))
			cl$Rho <- quote(ciarray(object$mcmc.draws$Rho, level = level))
			if (object$control$sample_df) {
				cl$df <- quote(quantile(object$mcmc.draws$df, c(sig.level/2, 1 - sig.level/2)))
			}
		}
		return(eval(cl[[parm]]))
	}
}