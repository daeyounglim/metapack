#' get compute the model comparison measures
#' @param object the output model from fitting a meta analysis/regression model
#' @param type the type of model comparison measures; DIC or LPML
#' @param verbose FALSE by default; If TRUE, then progress bar will appear
#' @param ncores the number of CPU cores to use for parallel processing; it must not exceed the number of existing cores
#' @importFrom parallel detectCores
#' @method model.comp bayesnmr
#' @export

"model.comp.bayesnmr" <- function(object, type="lpml", verbose=FALSE, ncores=NULL) {
	y <- object$Outcome
	npt <- object$Npt
	x <- object$XCovariate
	z <- object$ZCovariate
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
		tt <- table(iarm)
		idx <- which(tt == 1)
		if (length(idx) > 0) {
			nm <- names(tt)[idx]
			l <- c()
			for (i in 1:length(ids)) {
				if (iarm[i] %in% nm) {
					l <- c(l, ids[i])
				}
			}
			l <- unique(l)
			l <- l + 1
		}
		gof$lpml <- sum(gof$logcpo[-l])
	} else if (type == "pearson") {
		gof <- .Call('_metapack_calc_modelfit_pearson',
					 as.double(y),
					 as.matrix(x),
					 as.matrix(z),
					 as.integer(ids),
					 as.integer(iarm),
					 as.double(npt),
					 as.matrix(object$mcmc.draws$theta),
					 as.matrix(object$mcmc.draws$sig2),
					 as.matrix(object$mcmc.draws$phi),
					 as.matrix(object$mcmc.draws$lam),
					 as.array(object$mcmc.draws$Rho),
					 as.integer(K),
					 as.integer(nT),
					 as.integer(nkeep),
					 as.logical(verbose),
					 as.integer(ncores))
	}
	class(gof) <- "gofnmr"
	gof
}