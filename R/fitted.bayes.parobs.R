#' get fitted values
#' 
#' @param object the output model from fitting a meta analysis/regression model
#' @param conf.level credible level for interval estimation; set to 0.95 by default
#' @param HPD a logical argument indicating whether HPD intervals should be computed; if FALSE, equal-tail credible intervals are computed
#' @param ... additional arguments for fitted
#' @return a list of fitted values
#' @importFrom coda mcmc HPDinterval
#' @method fitted bayes.parobs
#' @export
"fitted.bayes.parobs" <- function(object, conf.level = 0.95, HPD = TRUE, ...) {
	out <- list()
	fmodel <- object$fmodel
	ypred <- list()
	ypred$mean <- apply(object$mcmc.draws$ypred, c(1,2), mean)
	ypred$sd <- apply(object$mcmc.draws$ypred, c(1,2), sd)
	if (object$fmodel == 1) {
		theta <- list()
		theta$mean <- rowMeans(object$mcmc.draws$theta)
		theta$sd <- apply(object$mcmc.draws$theta, 1, sd)
		Omega <- list()
		Omega$mean <- apply(object$mcmc.draws$Omega, c(1,2), mean)
		Omega$sd <- apply(object$mcmc.draws$Omega, c(1,2), sd)
		Sigma <- list()
		Sigma$mean <- apply(object$mcmc.draws$Sigma, c(1,2), mean)
		Sigma$sd <- apply(object$mcmc.draws$Sigma, c(1,2), sd)
	} else if (object$fmodel == 2) {
		theta <- list()
		theta$mean <- rowMeans(object$mcmc.draws$theta)
		theta$sd <- apply(object$mcmc.draws$theta, 1, sd)
		Omega <- list()
		Omega$mean <- apply(object$mcmc.draws$Omega, c(1,2), mean)
		Omega$sd <- apply(object$mcmc.draws$Omega, c(1,2), sd)
		Sigma <- list()
		Sigma$mean <- apply(object$mcmc.draws$Sigma, c(1,2), mean)
		Sigma$sd <- apply(object$mcmc.draws$Sigma, c(1,2), sd)
		R <- list()
		R$mean <- apply(object$mcmc.draws$R, c(1,2), mean)
		R$sd <- apply(object$mcmc.draws$R, c(1,2), sd)
	} else if (object$fmodel == 3) {
		theta <- list()
		theta$mean <- rowMeans(object$mcmc.draws$theta)
		theta$sd <- apply(object$mcmc.draws$theta, 1, sd)
		Omega <- list()
		Omega$mean <- apply(object$mcmc.draws$Omega, c(1,2), mean)
		Omega$sd <- apply(object$mcmc.draws$Omega, c(1,2), sd)
		Sigma <- list()
		Sigma$mean <- apply(object$mcmc.draws$Sigma, c(1,2), mean)
		Sigma$sd <- apply(object$mcmc.draws$Sigma, c(1,2), sd)
		delta <- list()
		delta$mean <- apply(object$mcmc.draws$delta, c(1,2), mean)
		delta$sd <- apply(object$mcmc.draws$delta, c(1,2), sd)
		Rho <- list()
		Rho$mean <- apply(object$mcmc.draws$Rho, c(1,2), mean)
		Rho$sd <- apply(object$mcmc.draws$Rho, c(1,2), sd)
		R <- list()
		R$mean <- apply(object$mcmc.draws$R, c(1,2), mean)
		R$sd <- apply(object$mcmc.draws$R, c(1,2), sd)
	} else if (object$fmodel == 4) {
		theta <- list()
		theta$mean <- rowMeans(object$mcmc.draws$theta)
		theta$sd <- apply(object$mcmc.draws$theta, 1, sd)
		Omega <- list()
		Omega$mean <- apply(object$mcmc.draws$Omega, c(1,2), mean)
		Omega$sd <- apply(object$mcmc.draws$Omega, c(1,2), sd)
		Sigma <- list()
		Sigma$mean <- apply(object$mcmc.draws$Sigma, c(1,2), mean)
		Sigma$sd <- apply(object$mcmc.draws$Sigma, c(1,2), sd)
		Delta <- list()
		Delta$mean <- rowMeans(object$mcmc.draws$delta)
		Delta$sd <- apply(object$mcmc.draws$delta, 1, sd)
		Rho <- list()
		Rho$mean <- apply(object$mcmc.draws$Rho, c(1,2), mean)
		Rho$sd <- apply(object$mcmc.draws$Rho, c(1,2), sd)
		Sigma0 <- list()
		Sigma0$mean <- apply(object$mcmc.draws$Sigma0, c(1,2), mean)
		Sigma0$sd <- apply(object$mcmc.draws$Sigma0, c(1,2), sd)
		R <- list()
		R$mean <- apply(object$mcmc.draws$R, c(1,2), mean)
		R$sd <- apply(object$mcmc.draws$R, c(1,2), sd)
	}
	sig.level <- 1 - conf.level

	if (HPD) {
		theta.hpd <- coda::HPDinterval(mcmc(t(object$mcmc.draws$theta), end=object$mcmc$nkeep), prob=conf.level)
		theta$lower <- theta.hpd[,1]
		theta$upper <- theta.hpd[,2]

		Omega.hpd <- hpdarray(object$mcmc.draws$Omega, conf.level = conf.level)
		Omega$lower <- Omega.hpd[,,1]
		Omega$upper <- Omega.hpd[,,2]
		Sigma.hpd <- hpdarray(object$mcmc.draws$Sigma, conf.level = conf.level)
		Sigma$lower <- Sigma.hpd[,,1]
		Sigma$upper <- Sigma.hpd[,,2]

		if (fmodel >= 2) {
			R.hpd <- hpdarray(object$mcmc.draws$R, conf.level = conf.level)
			R$lower <- R.hpd[,,1]
			R$upper <- R.hpd[,,2]

			if (fmodel == 3) {
				delta.hpd <- hpdarray(object$mcmc.draws$delta, conf.level = conf.level)
				delta$lower <- delta.hpd[,,1]
				delta$upper <- delta.hpd[,,2]

				Rho.hpd <- hpdarray(object$mcmc.draws$Rho, conf.level = conf.level)
				Rho$lower <- Rho.hpd[,,1]
				Rho$upper <- Rho.hpd[,,2]
			} else if (fmodel == 4) {
				Delta.hpd <- coda::HPDinterval(mcmc(t(object$mcmc.draws$delta), end=object$mcmc$nkeep), prob=conf.level)
				Delta$lower <- Delta.hpd[,1]
				Delta$upper <- Delta.hpd[,2]

				Rho.hpd <- hpdarray(object$mcmc.draws$Rho, conf.level = conf.level)
				Rho$lower <- Rho.hpd[,,1]
				Rho$upper <- Rho.hpd[,,2]
				
				Sigma0.hpd <- hpdarray(object$mcmc.draws$Sigma0, conf.level = conf.level)
				Sigma0$lower <- Sigma0.hpd[,,1]
				Sigma0$upper <- Sigma0.hpd[,,2]
			}
		}
	} else {
		theta$lower <- apply(object$mcmc.draws$theta, 1, function(xx) quantile(xx, prob = sig.level/2))
		theta$upper <- apply(object$mcmc.draws$theta, 1, function(xx) quantile(xx, prob = 1-sig.level/2))
		Omega$lower <- apply(object$mcmc.draws$Omega, 3, function(xx) quantile(xx, prob = sig.level/2))[lower.tri(object$mcmc.draws$Omega[,,1])]
		Omega$upper <- apply(object$mcmc.draws$Omega, 3, function(xx) quantile(xx, prob = 1-sig.level/2))[lower.tri(object$mcmc.draws$Omega[,,1])]
		Sigma$lower <- apply(object$mcmc.draws$Sigma, 3, function(xx) quantile(xx, prob = sig.level/2))[lower.tri(object$mcmc.draws$Sigma[,,1])]
		Sigma$upper <- apply(object$mcmc.draws$Sigma, 3, function(xx) quantile(xx, prob = 1-sig.level/2))[lower.tri(object$mcmc.draws$Sigma[,,1])]

		if (fmodel >= 2) {
			R$lower <- apply(object$mcmc.draws$R, 3, function(xx) quantile(xx, prob = sig.level/2))
			R$upper <- apply(object$mcmc.draws$R, 3, function(xx) quantile(xx, prob = 1-sig.level/2))

			if (fmodel == 3) {
				delta$lower <- apply(object$mcmc.draws$delta, 3, function(xx) quantile(xx, prob = sig.level/2))
				delta$upper <- apply(object$mcmc.draws$delta, 3, function(xx) quantile(xx, prob = 1-sig.level/2))

				Rho$lower <- apply(object$mcmc.draws$Rho, 3, function(xx) quantile(xx, prob = sig.level/2))
				Rho$upper <- apply(object$mcmc.draws$Rho, 3, function(xx) quantile(xx, prob = 1-sig.level/2))
			} else if (fmodel == 4) {
				Delta$lower <- apply(object$mcmc.draws$delta, 1, function(xx) quantile(xx, prob = sig.level/2))
				Delta$upper <- apply(object$mcmc.draws$delta, 1, function(xx) quantile(xx, prob = 1-sig.level/2))

				Rho$lower <- apply(object$mcmc.draws$Rho, 3, function(xx) quantile(xx, prob = sig.level/2))
				Rho$upper <- apply(object$mcmc.draws$Rho, 3, function(xx) quantile(xx, prob = 1-sig.level/2))
				
				Sigma0$lower <- apply(object$mcmc.draws$Sigma0, 3, function(xx) quantile(xx, prob = sig.level/2))
				Sigma0$upper <- apply(object$mcmc.draws$Sigma0, 3, function(xx) quantile(xx, prob = 1-sig.level/2))
			}
		}
	}

	out <- object
	out$conf.level <- conf.level
	out$hpd <- HPD
	out$ypred <- ypred
	out$theta <- theta
	out$Sigma <- Sigma
	out$Omega <- Omega
	if (fmodel >= 2) {
		out$R <- R
		if (fmodel == 3) {
			out$delta <- delta
			out$Rho <- Rho
		} else if (fmodel == 4) {
			out$Delta <- Delta
			out$Rho <- Rho
			out$Sigma0 <- Sigma0
		}
	}

	class(out) <- "fitted.bayes.parobs"
	out
}