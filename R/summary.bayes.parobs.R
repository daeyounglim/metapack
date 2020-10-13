#' `summary` method for class "`bayes.parobs`"
#' @param object the output model from fitting a meta analysis/regression model
#' @param ... additional arguments for summary
#' @importFrom coda mcmc HPDinterval
#' @return print summary for the model fit
#' @method summary bayes.parobs
#' @md
#' @export
"summary.bayes.parobs" <- function(object, ...) {
	digits <- max(3, getOption("digits") - 3)
	if (class(object) != "bayes.parobs") {
		stop("'summary.bayes.parobs' designed for 'bayes.parobs' objects")
	}
	theta <- list()
	theta$mean <- rowMeans(object$mcmc.draws$theta)
	theta$sd <- apply(object$mcmc.draws$theta, 1, sd)

	sig.level <- 1 - 0.95

	theta.hpd <- coda::HPDinterval(mcmc(t(object$mcmc.draws$theta), end=object$mcmc$nkeep), prob=0.95)
	theta$lower <- theta.hpd[,1]
	theta$upper <- theta.hpd[,2]
	r <- cbind(theta$mean, theta$sd, theta$lower, theta$upper)
	colnames(r) <- c("Post.Mean", "Std.Dev", "HPD(Lower)", "HPD(Upper)")
	xcc <- if (!is.null(colnames(object$XCovariate))) colnames(object$XCovariate) else paste0("beta_", 1:ncol(object$XCovariate))
	wcc <- if (!is.null(colnames(object$WCovariate))) colnames(object$WCovariate) else paste0("gam_", 1:ncol(object$WCovariate))

	J <- ncol(object$Outcome)
	if (is.null(object$group)) {
		rownames(r) <- c(paste0(rep(xcc, J), "_", rep(1:J, each=length(xcc))), paste0(rep(wcc, J), "_", rep(1:J, each=length(wcc))))
	} else {
		# rownames(r) <- c(paste0(rep(xcc, J), "_", rep(1:J, each=length(xcc))), paste0(rep(wcc, J),"*(1-2nd)", "_", rep(1:J, each = length(wcc))), paste0(rep(wcc, J),"*2nd", "_", rep(1:J, each = length(wcc))))
		rownames(r) <- c(paste0(rep(xcc, J), "_", rep(1:J, each=length(xcc))),
						 paste0(rep(wcc, 2*J), rep(rep(c("*(1-2nd)", "*2nd"), each = length(wcc)), J), "_", rep(1:J, each = 2*length(wcc))))
	}
	cat("\nPosterior inference in multivariate meta-regression models\n")
	cat("Fixed-effects:\n")
	r <- round(r, digits=digits)
	print.default(r, print.gap = 2)
	invisible(r)
}