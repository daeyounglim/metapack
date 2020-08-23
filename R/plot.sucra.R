#' plot the surface under the cumulative ranking curve (SUCRA)
#' @param object the output model from fitting a network meta analysis/regression model
#' @method plot sucra
#' @export
"plot.sucra" <- function(object, ...) {
	nT <- length(object$SUCRA)
	cumeffectiveness <- apply(object$rankProb, 2, cumsum)
	names <- object$names
	oask <- devAskNewPage(TRUE)
	for (TRT in 1:nT) {
		Area=round(object$SUCRA[TRT],3)
		plot(1:nT,type="n",ylim=c(0,1), xaxt="n",
		   xlab=paste("Rank of",as.character(names[TRT])), 
		   ylab="Probability", main = paste0("SUCRA = ", Area), ...)
		axis(1,at=1:nT,label=names)
		lines(lwd=2,lty=1,c(1,c(1:c(nT)),nT), 
		    cumeffectiveness[c(1,1:c(nT),c(nT)),TRT], col=rgb(0, 157, 114, maxColorValue = 255))
		lines(lwd=2,lty=3,col="#eab159",c(1,c(1:c(nT)),nT), 
		    object$rankProb[c(1,1:c(nT),c(nT)),TRT])
	}
	on.exit(devAskNewPage(oask))
}