#' get goodness of fit 
#' @param x the output model from fitting a meta analysis/regression model
#' @param ... additional parameters for plot
#' @importFrom graphics par
#' @importFrom stats density 
#' @method plot bayesnmr
#' @export
#' 
#' 
'plot.bayesnmr' <- function(x, ...) {
	nkeep <- x$mcmc$nkeep
	param <- x$mcmc.draws$theta
	p <- nrow(param)
	if (x$scale_x) {
		xcols <- ncol(x$Covariate)
		tlength <- nrow(param)
		trlength <- tlength - xcols
		tscale <- c(apply(x$Covariate, 2, sd), rep(1, trlength))
	} else {
		tlength <- nrow(param)
		tscale <- rep(1, tlength)
	}
	theta.post <- vapply(1:x$mcmc$nkeep, function(ikeep) {
		param[,ikeep] / tscale
	}, FUN.VALUE = numeric(tlength))
	wname <- rownames(param)	
	if (p == 2) {
		par(mfcol = c(2, 2))
		for (i in 1:p) {
			plot(1:nkeep, param[i, ], xlab = "Iteration", ylab = "", main = wname[i], type = "l")
			plot(density(param[i, ]), main = "")
		}
	} else if (p == 3) {
		par(mfcol = c(2, 3))
		for (i in 1:p) {
			plot(1:nkeep, param[i, ], xlab = "Iteration", ylab = "", main = wname[i], type = "l")
			plot(density(param[i, ]), main = "")
		}
	} else if (p == 4) {
		par(mfcol = c(2, 2))
		for (i in 1:p) {
			plot(1:nkeep, param[i, ], xlab = "Iteration", ylab = "", main = wname[i], type = "l")
			plot(density(param[i, ]), main = "")
			par(...)
		}
	} else {
		par(mfcol = c(2, 3))
		for (i in 1:p) {
			plot(1:nkeep, param[i, ], xlab = "Iteration", ylab = "", main = wname[i], type = "l")
			plot(density(param[i, ]), main = "")
			par(...)
		}
	}
}