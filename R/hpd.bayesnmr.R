#' get the highest posterior density (HPD) interval
#' @param object the output model from fitting a meta analysis/regression model
#' @param prob the probability which the HPD interval will cover
#' @importFrom coda mcmc HPDinterval
#' @return dataframe containing HPD intervals for the parameters
#' @method confint bayesnmr
#' @export

"confint.bayesnmr" <- function(object, prob = 0.95) {
	out <- list()
	out$theta.hpd <- coda::HPDinterval(mcmc(t(object$mcmc.draws$theta), end=object$mcmc$nkeep), prob=prob)
	out$phi.hpd <- coda::HPDinterval(mcmc(t(object$mcmc.draws$phi), end=object$mcmc$nkeep), prob=prob)
	out$gam.hpd <- coda::HPDinterval(mcmc(t(object$mcmc.draws$gam), end=object$mcmc$nkeep), prob=prob)
	out$sig2.hpd <- coda::HPDinterval(mcmc(t(object$mcmc.draws$sig2), end=object$mcmc$nkeep), prob=prob)
	out$Rho.hpd <- coda::HPDinterval(mcmc(t(apply(object$mcmc.draws$Rho, 3, function(xx) xx[lower.tri(xx)])), end=object$mcmc$nkeep), prob=prob)
	class(out) <- "bayesnmr.hpd"
	out
}