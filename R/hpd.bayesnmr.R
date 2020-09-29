#' get the highest posterior density (HPD) interval
#' @param object the output model from fitting a meta analysis/regression model
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the probability which the HPD interval will cover
#' @param ... additional arguments for confint
#' @importFrom coda mcmc HPDinterval
#' @return dataframe containing HPD intervals for the parameters
#' @method confint bayesnmr
#' @export

"confint.bayesnmr" <- function(object, parm, level = 0.95, ...) {
	if (missing(parm)) {
		out <- list()
		out$theta <- coda::HPDinterval(mcmc(t(object$mcmc.draws$theta), end=object$mcmc$nkeep), prob=level)
		out$phi <- coda::HPDinterval(mcmc(t(object$mcmc.draws$phi), end=object$mcmc$nkeep), prob=level)
		out$gam <- coda::HPDinterval(mcmc(t(object$mcmc.draws$gam), end=object$mcmc$nkeep), prob=level)
		out$sig2 <- coda::HPDinterval(mcmc(t(object$mcmc.draws$sig2), end=object$mcmc$nkeep), prob=level)
		out$Rho <- coda::HPDinterval(mcmc(t(apply(object$mcmc.draws$Rho, 3, function(xx) xx[lower.tri(xx)])), end=object$mcmc$nkeep), prob=level)
		class(out) <- "bayesnmr.hpd"
		return(out)
	} else {
		cl <- list()
		cl$theta <-quote(coda::HPDinterval(mcmc(t(object$mcmc.draws$theta), end=object$mcmc$nkeep), prob=level))
		cl$phi <- quote(coda::HPDinterval(mcmc(t(object$mcmc.draws$phi), end=object$mcmc$nkeep), prob=level))
		cl$gam <- quote(coda::HPDinterval(mcmc(t(object$mcmc.draws$gam), end=object$mcmc$nkeep), prob=level))
		cl$sig2 <- quote(coda::HPDinterval(mcmc(t(object$mcmc.draws$sig2), end=object$mcmc$nkeep), prob=level))
		cl$Rho <- quote(coda::HPDinterval(mcmc(t(apply(object$mcmc.draws$Rho, 3, function(xx) xx[lower.tri(xx)])), end=object$mcmc$nkeep), prob=level))
		return(eval(cl[[parm]]))
	}
}