#' get the highest posterior density (HPD) interval
#' @param object the output model from fitting a meta analysis/regression model
#' @param prob the probability which the HPD interval will cover
#' @importFrom coda mcmc HPDinterval
#' @return dataframe containing HPD intervals for the parameters
#' @method hpd bayes.parobs
#' @export

"hpd.bayes.parobs" <- function(object, prob = 0.95) {
	out <- list()

	out$theta.hpd <- coda::HPDinterval(mcmc(t(object$mcmc.draws$theta), end=object$mcmc$nkeep), prob=prob)
	out$Omega.hpd <- hpdarray(object$mcmc.draws$Omega, conf.level = prob)
	out$Sigma.hpd <- hpdarray(object$mcmc.draws$Sigma, conf.level = prob)
	if (object$fmodel >= 2) {
		out$R.hpd <- hpdarray(object$mcmc.draws$R, conf.level = prob)
		if (object$fmodel == 3) {
			out$delta.hpd <- hpdarray(object$mcmc.draws$delta, conf.level = prob)
			out$Rho.hpd <- hpdarray(object$mcmc.draws$Rho, conf.level = prob)
		} else if (object$fmodel == 4) {
			out$Delta.hpd <- coda::HPDinterval(mcmc(t(object$mcmc.draws$Delta), end=object$mcmc$nkeep), prob=prob)
			out$Rho.hpd <- hpdarray(object$mcmc.draws$Rho, conf.level = prob)
			out$Sigma0.hpd <- hpdarray(object$mcmc.draws$Sigma0, conf.level = prob)
		}
	}
	class(out) <- "bayes.parobs.hpd"
	out
}