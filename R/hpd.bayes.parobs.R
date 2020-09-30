#' get the highest posterior density (HPD) interval
#' @param object the output model from fitting a meta analysis/regression model
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the probability which the HPD interval will cover
#' @param ... additional arguments for confint
#' @importFrom coda mcmc HPDinterval
#' @return dataframe containing HPD intervals for the parameters
#' @method confint bayes.parobs
#' @export

"confint.bayes.parobs" <- function(object, parm, level = 0.95, ...) {
	if (missing(parm)) {
		out <- list()
		out$ypred <- hpdarray(object$mcmc.draws$ypred, conf.level=level)
		out$theta <- coda::HPDinterval(mcmc(t(object$mcmc.draws$theta), end=object$mcmc$nkeep), prob=level)
		out$Omega <- hpdarray(object$mcmc.draws$Omega, conf.level = level)
		out$Sigma <- hpdarray(object$mcmc.draws$Sigma, conf.level = level)
		if (object$fmodel >= 2) {
			out$R <- hpdarray(object$mcmc.draws$R, conf.level = level)
			if (object$fmodel == 4) {
				out$delta <- hpdarray(object$mcmc.draws$delta, conf.level = level)
				out$Rho <- hpdarray(object$mcmc.draws$Rho, conf.level = level)
			} else if (object$fmodel == 5) {
				out$delta <- coda::HPDinterval(mcmc(t(object$mcmc.draws$delta), end=object$mcmc$nkeep), prob=level)
				out$Rho <- hpdarray(object$mcmc.draws$Rho, conf.level = level)
				out$Sigma0 <- hpdarray(object$mcmc.draws$Sigma0, conf.level = level)
			}
		}
		class(out) <- "bayes.parobs.hpd"
		return(out)
	} else {
		cl <- list()
		cl$ypred <- quote(hpdarray(object$mcmc.draws$ypred, conf.level=level))
		cl$theta <- quote(coda::HPDinterval(mcmc(t(object$mcmc.draws$theta), end=object$mcmc$nkeep), prob=level))
		cl$Omega <- quote(hpdarray(object$mcmc.draws$Omega, conf.level = level))
		cl$Sigma <- quote(hpdarray(object$mcmc.draws$Sigma, conf.level = level))
		if (object$fmodel >= 2) {
			cl$R <- quote(hpdarray(object$mcmc.draws$R, conf.level = level))
			if (object$fmodel == 4) {
				cl$delta <- quote(hpdarray(object$mcmc.draws$delta, conf.level = level))
				cl$Rho <- quote(hpdarray(object$mcmc.draws$Rho, conf.level = level))
			} else if (object$fmodel == 5) {
				cl$delta <- quote(coda::HPDinterval(mcmc(t(object$mcmc.draws$delta), end=object$mcmc$nkeep), prob=level))
				cl$Rho <- quote(hpdarray(object$mcmc.draws$Rho, conf.level = level))
				cl$Sigma0 <- quote(hpdarray(object$mcmc.draws$Sigma0, conf.level = level))
			}
		}
		return(eval(cl[[parm]]))
	}
}