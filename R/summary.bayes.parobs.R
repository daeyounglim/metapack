#' `summary` method for class "`bayes.parobs`"
#' @param object the output model from fitting a meta analysis/regression model
#' @param ... additional arguments for summary
#' @return print summary for the model fit
#' @md
#' @export
"summary.bayes.parobs" <- function(object, ...) {
	digits <- max(3, getOption("digits") - 3)
	if (!inherits(object, "bayes.parobs")) {
		stop(paste(sQuote('summary.bayes.parobs'), "designed for", sQuote('bayes.parobs'), "objects"))
	}

	if (inherits(object, "bsynthesis")) {
		cat("\nCall:\n", paste(deparse(object$call), sep = "\n", 
	        collapse = "\n"), "\n", sep = "")
	}
	if (object$scale_x) {
		J <- ncol(object$Outcome)
		xcols <- ncol(object$XCovariate)
		tlength <- nrow(object$mcmc.draws$theta)
		trlength <- tlength - xcols * J
		tscale <- c(rep(unname(attributes(object$XCovariate)$`scaled:scale`), J), rep(1, trlength))
	} else {
		tlength <- nrow(object$mcmc.draws$theta)
		tscale <- rep(1, tlength)
	}
	theta <- list()
	theta.post <- vapply(1:object$mcmc$nkeep, function(ikeep) {
		object$mcmc.draws$theta[,ikeep] / tscale
	}, FUN.VALUE = numeric(tlength))
	theta <- list()
	theta$mean <- rowMeans(theta.post)
	theta$sd <- apply(theta.post, 1, sd)

	sig.level <- 1 - 0.95

	theta.hpd <- mhpd(theta.post, 0.95)
	theta$lower <- theta.hpd[,1]
	theta$upper <- theta.hpd[,2]
	r <- cbind(theta$mean, theta$sd, theta$lower, theta$upper)
	colnames(r) <- c("Post.Mean", "Std.Dev", "HPD(Lower)", "HPD(Upper)")
	xcc <- if (!is.null(colnames(object$XCovariate))) colnames(object$XCovariate) else paste0("beta", 1:ncol(object$XCovariate))
	wcc <- if (!is.null(colnames(object$WCovariate))) colnames(object$WCovariate) else paste0("gam", 1:ncol(object$WCovariate))

	J <- ncol(object$Outcome)
	if (is.null(object$group)) {
		rownames(r) <- c(paste0(rep(xcc, J), "_", rep(1:J, each=length(xcc))), paste0(rep(wcc, J), "_", rep(1:J, each=length(wcc))))
	} else {
		rownames(r) <- c(paste0(rep(xcc, J), "_", rep(1:J, each=length(xcc))),
						 paste0(rep(wcc, 2*J), rep(rep(c("*(1-2nd)", "*2nd"), each = length(wcc)), J), "_", rep(1:J, each = 2*length(wcc))))
	}
	cat("Fixed-effects:\n")
	r <- round(r, digits=digits)
	print.default(r, print.gap = 2)
	invisible(r)
}
