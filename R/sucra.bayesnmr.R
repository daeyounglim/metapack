#' get surface under the cumulative ranking curve (SUCRA)
#' @param object the output model from fitting a network meta analysis/regression model
#' @return a matrix containing SUCRA
#' @method sucra bayesnmr
#' @export
"sucra.bayesnmr" <- function(object) {
	nT <- object$nT
	nx <- ncol(object$Covariate)
	gam.posterior <- object$mcmc.draws$gam[(nx+1):(nx+nT),] # treatment effects
	nkeep <- object$mcmc$nkeep

	TRT.rank=matrix(0, nkeep, nT)
	for (i in 1:nkeep){
	  	TRT.rank[i,] <- sort(gam.posterior[,i],index.return=TRUE)$ix
	}

	R = matrix(0, nT, nT)
	for (j in 1:nT) {
	  	for (k in 1:nT) {
	    	R[j,k] <- sum(as.numeric(TRT.rank[,j] == k))
	  	}
	}
	row.sum <- rowSums(R)
	Rank.prob <- R / row.sum

	cumeffectiveness <- t(apply(Rank.prob, 2, cumsum))
	SUCRA <- numeric(nT)
	for (i in 1:nT) {
		SUCRA[i] <- sum(cumeffectiveness[i, 1:(nT-1)]) / (nT-1)
	}
	out <- list(SUCRA=SUCRA, rankProb=Rank.prob, names=object$TrtLabels)
	class(out) <- "sucra"
	out
} 