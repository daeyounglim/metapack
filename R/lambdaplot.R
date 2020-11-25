#' plot the Gamma Quantile-Quantile for gamma(nu/2, nu/2) with NMA object
#' @param x the output model from fitting a network meta analysis/regression model
#' @param xlab plot labels
#' @param ylab plot labels
#' @param ... additional arguments for plot
#' @export
"lambdaplot" <- function(x, xlab = NULL, ylab = NULL, type = "boxplot", ...) {
	UseMethod("lambdaplot", x)
} 