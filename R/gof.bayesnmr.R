#' get goodness of fit 
#' @param object the output model from fitting a meta analysis/regression model
#' @param type the type of goodness of fit to compute; DIC or LPML
#' @param verbose FALSE by default; If TRUE, then progress bar will appear
#' @param ncores the number of CPU cores to use for parallel processing; it must not exceed the number of existing cores
#' @param h the interval width for trapezoidal rule
#' @importFrom parallel detectCores
#' @method gof bayesnmr
#' @export

"gof.bayesnmr" <- function(object, type="lpml", verbose=FALSE, ncores=NULL, h = 0.5) {
	y <- object$Outcome
	npt <- object$Npt
	x <- object$Covariate
	z <- object$z
	ids <- object$Trial
	iarm <- object$Treat
	K <- object$K
	nT <- object$nT
	nkeep <- object$mcmc$nkeep
	nu <- object$prior$df

	if (!is.null(ncores)) {
		ncores_ <- parallel::detectCores()
		if (ncores > ncores_) {
			stop(paste0("The number of cores must not exceed ", ncores_))
		}
	} else {
		ncores <- parallel::detectCores()
	}


	if (type == "dic") {
		gof <- .Call(`_metapack_calc_modelfit_dic`,
					 as.double(y),
					 as.matrix(x),
					 as.matrix(z),
					 as.integer(ids),
					 as.integer(iarm),
					 as.double(npt),
					 as.double(object$mcmc.draws$df),
					 as.double(nu),
					 as.matrix(object$mcmc.draws$theta),
					 as.matrix(object$mcmc.draws$sig2),
					 as.matrix(object$mcmc.draws$phi),
					 as.matrix(object$mcmc.draws$lam),
					 as.array(object$mcmc.draws$Rho),
					 as.integer(K),
					 as.integer(nT),
					 as.integer(nkeep),
					 as.logical(object$control$sample_df),
					 as.logical(verbose),
					 as.integer(ncores))
	} else if (type == "lpml") {
		gof <- .Call(`_metapack_calc_modelfit_lpml`,
					 as.double(y),
					 as.matrix(x),
					 as.matrix(z),
					 as.integer(ids),
					 as.integer(iarm),
					 as.double(npt),
					 as.double(object$mcmc.draws$df),
					 as.double(nu),
					 as.matrix(object$mcmc.draws$theta),
					 as.matrix(object$mcmc.draws$sig2),
					 as.matrix(object$mcmc.draws$phi),
					 as.matrix(object$mcmc.draws$lam),
					 as.array(object$mcmc.draws$Rho),
					 as.integer(K),
					 as.integer(nT),
					 as.integer(nkeep),
					 as.logical(object$control$sample_df),
					 as.logical(verbose),
					 as.integer(ncores))
	}

	class(gof) <- "gofnmr"
	gof
}