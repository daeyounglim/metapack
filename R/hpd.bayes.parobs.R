#' get the highest posterior density (HPD) interval
#' @param object the output model from fitting a meta analysis/regression model
#' @param prob the probability which the HPD interval will cover
#' @importFrom coda mcmc HPDinterval
#' @return dataframe containing HPD intervals for the parameters
#' @method confint bayes.parobs
#' @export

"confint.bayes.parobs" <- function(object, level = 0.95) {
	out <- list()

	out$ypred.hpd <- hpdarray(object$mcmc.draws$ypred, conf.level=level)
	out$theta.hpd <- coda::HPDinterval(mcmc(t(object$mcmc.draws$theta), end=object$mcmc$nkeep), prob=level)
	out$Omega.hpd <- hpdarray(object$mcmc.draws$Omega, conf.level = level)
	out$Sigma.hpd <- hpdarray(object$mcmc.draws$Sigma, conf.level = level)
	if (object$fmodel >= 2) {
		out$R.hpd <- hpdarray(object$mcmc.draws$R, conf.level = level)
		if (object$fmodel == 3) {
			out$delta.hpd <- hpdarray(object$mcmc.draws$delta, conf.level = level)
			out$Rho.hpd <- hpdarray(object$mcmc.draws$Rho, conf.level = level)
		} else if (object$fmodel == 4) {
			out$delta.hpd <- coda::HPDinterval(mcmc(t(object$mcmc.draws$delta), end=object$mcmc$nkeep), prob=level)
			out$Rho.hpd <- hpdarray(object$mcmc.draws$Rho, conf.level = level)
			out$Sigma0.hpd <- hpdarray(object$mcmc.draws$Sigma0, conf.level = level)
		}
	}
	class(out) <- "bayes.parobs.hpd"
	out
}