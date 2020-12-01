#' plot the side-by-side box plot or gamma Quantile-Quantile for gamma(nu/2, nu/2) with NMA object
#' @param x the output model from fitting a network meta analysis/regression model
#' @param xlab plot x label
#' @param ylab plot y label
#' @param type type of plot; either boxplot or qqplot
#' @param ... additional arguments for plot
#' @export
"lambdaplot" <- function(x, xlab = NULL, ylab = NULL, type = "boxplot", ...) {
	UseMethod("lambdaplot", x)
} 