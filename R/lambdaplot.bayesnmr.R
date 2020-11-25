#' plot the Gamma Quantile-Quantile for gamma(nu/2, nu/2) with NMA object
#' @param x the output model from fitting a network meta analysis/regression model
#' @param xlab plot labels
#' @param ylab plot labels
#' @param type type of plot; either boxplot or qqplot
#' @param ... additional arguments for plot
#' @importFrom stats qqplot ppoints
#' @importFrom graphics par boxplot
#' @method lambdaplot bayesnmr
#' @export
"lambdaplot.bayesnmr" <- function(x, xlab = NULL, ylab = NULL, type = "boxplot", ...) {
	if (type == 'qqplot') {
		old_pars <- par(mfcol = c(2, 2))
		on.exit(par(old_pars))
		if (is.null(ylab)) {
			ylab_s <- "Theoretical Gamma Quantiles"
		} else {
			ylab_s <- ylab
		}
		for (k in 1:(x$K)) {
			if (is.null(xlab)) {
				xlab_s <- bquote("Sample Quantiles " ~ lambda[.(k)])
			} else {
				xlab_s <- xlab
			}
			lamk <- x$mcmc.draws$lam[k,]
			if (x$control$sample_df) {
				df <- mean(x$mcmc.draws$df)
			} else {
				df <- x$prior$df
			}
			stats::qqplot(y = qgamma(ppoints(lamk), shape = 0.5 * df, rate = 0.5 * df),
				x = lamk, xlab = xlab_s, ylab = ylab_s, ...)
			abline(a = 0, b = 1, col = "dodgerblue", lwd = 2)
		}
	} else if (type == "boxplot") {
		if (x$control$sample_df) {
			lam <- sapply(1:(x$mcmc$nkeep), function(ikeep) {
				(x$mcmc.draws$lam[,ikeep] - 1) * sqrt(0.5 * x$mcmc.draws$df[ikeep])
			})
		} else {
			lam <- (x$mcmc.draws$lam - 1) * sqrt(0.5 * x$prior$df)
		}
		if (is.null(xlab)) {
			xlab_s <- "Trial"
		} else {
			xlab_s <- xlab
		}
		if (is.null(ylab)){
			ylab_s <- expression(lambda[k]^"\u2217")
		} else {
			ylab_s <- ylab
		}
		boxplot(t(lam), xlab = xlab_s, ylab = ylab_s)
	} else {
		stop("Invalid type. Choose between `boxplot` and `qqplot`.")
	}
}