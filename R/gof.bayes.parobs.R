#' get goodness of fit 
#' @param object the output model from fitting a meta analysis/regression model
#' @param type the type of goodness of fit to compute; DIC or LPML
#' @param verbose FALSE by default; If TRUE, then progress bar will appear
#' @method gof bayes.parobs
#' @export
"gof.bayes.parobs" <- function(object, type="lpml", ncores=NULL, verbose=FALSE) {
	nkeep <- object$mcmc$nkeep
	
	Sigmahat <- apply(object$mcmc.draws$Sigma, c(1,2), mean)
	Omegahat <- apply(object$mcmc.draws$Omega, c(1,2), mean)
	thetahat <- rowMeans(object$mcmc.draws$theta)

	if (!is.null(ncores)) {
		ncores_ <- parallel::detectCores()
		if (ncores > ncores_) {
			stop(paste0("The number of cores must not exceed ", ncores_))
		}
	} else {
		ncores <- parallel::detectCores()
	}

	if (type == "dic") {
		gof <- .Call(`_metapack_dic_parcov`,
					  as.matrix(object$Outcome),
			  		  as.matrix(object$XCovariate),
			  		  as.matrix(object$WCovariate),
			  		  as.vector(object$Npt),
			  		  object$mcmc.draws$Sigma,
			  		  object$mcmc.draws$Omega,
			  		  as.matrix(object$mcmc.draws$theta),
			  		  as.double(thetahat),
			  		  as.matrix(Sigmahat),
			  		  as.matrix(Omegahat),
			  		  as.integer(object$fmodel),
			  		  as.integer(nkeep),
			  		  as.logical(verbose),
			  		  as.integer(ncores))
	} else if (type == "lpml") {
		gof <- .Call(`_metapack_lpml_parcov`,
					  as.matrix(object$Outcome),
			  		  as.matrix(object$XCovariate),
			  		  as.matrix(object$WCovariate),
			  		  as.vector(object$Npt),
			  		  object$mcmc.draws$Sigma,
			  		  object$mcmc.draws$Omega,
			  		  as.matrix(object$mcmc.draws$theta),
			  		  as.double(thetahat),
			  		  as.matrix(Sigmahat),
			  		  as.matrix(Omegahat),
			  		  as.integer(object$fmodel),
			  		  as.integer(nkeep),
			  		  as.logical(verbose),
			  		  as.integer(ncores))
	}

	class(gof) <- "gofparobs"
	gof
}