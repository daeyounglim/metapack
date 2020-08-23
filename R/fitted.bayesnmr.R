#' get fitted values
#' 
#' @param object the output model from fitting a meta analysis/regression model
#' @return a list of fitted values
#' @importFrom coda mcmc HPDinterval
#' @method fitted bayesnmr
#' @export
"fitted.bayesnmr" <- function(object, conf.level = 0.95, HPD = TRUE, ...) {
	
	out <- list()

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
	sig2$mean <- rowMeans(object$mcmc.draws$sig2)
	sig2$sd <- apply(object$mcmc.draws$sig2, 1, sd)
	Rho$mean <- apply(object$mcmc.draws$Rho, c(1,2), mean)
	Rho$sd <- apply(object$mcmc.draws$Rho, c(1,2), sd)

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

		sig2.hpd <- coda::HPDinterval(mcmc(t(object$mcmc.draws$sig2), end=object$mcmc$nkeep), prob=conf.level)
		sig2$lower <- sig2.hpd[,1]
		sig2$upper <- sig2.hpd[,2]

		Rho.hpd <- coda::HPDinterval(mcmc(t(apply(object$mcmc.draws$Rho, 3, function(xx) xx[lower.tri(xx)])), end=object$mcmc$nkeep), prob=conf.level)
		Rho$lower <- Rho.hpd[,1]
		Rho$upper <- Rho.hpd[,2]
	} else {
		beta$lower <- apply(object$mcmc.draws$beta, 1, function(xx) quantile(xx, prob = sig.level/2))
		beta$upper <- apply(object$mcmc.draws$beta, 1, function(xx) quantile(xx, prob = 1-sig.level/2))

		phi$lower <- apply(object$mcmc.draws$phi, 1, function(xx) quantile(xx, prob = sig.level/2))
		phi$upper <- apply(object$mcmc.draws$phi, 1, function(xx) quantile(xx, prob = 1-sig.level/2))

		gam$lower <- apply(object$mcmc.draws$gam, 1, function(xx) quantile(xx, prob = sig.level/2))
		gam$upper <- apply(object$mcmc.draws$gam, 1, function(xx) quantile(xx, prob = 1-sig.level/2))

		sig2$lower <- apply(object$mcmc.draws$sig2, 1, function(xx) quantile(xx, prob = sig.level/2))
		sig2$upper <- apply(object$mcmc.draws$sig2, 1, function(xx) quantile(xx, prob = 1-sig.level/2))

		Rho$lower <- apply(object$mcmc.draws$Rho, 3, function(xx) quantile(xx, prob = sig.level/2))[lower.tri(object$mcmc.draws$Rho[,,1])]
		Rho$upper <- apply(object$mcmc.draws$Rho, 3, function(xx) quantile(xx, prob = 1-sig.level/2))[lower.tri(object$mcmc.draws$Rho[,,1])]
	}

	out <- object
	out$conf.level <- conf.level
	out$hpd <- hpd
	out$beta <- beta
	out$phi <- phi
	out$gam <- gam
	out$sig2 <- sig2
	out$Rho <- Rho
	class(out) <- "fitted.bayesnmr"
	out
}