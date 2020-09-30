#' plot the surface under the cumulative ranking curve (SUCRA)
#' @param x the output model from fitting a network meta analysis/regression model
#' @param legend.position the position of the legend that will be passed onto ggplot
#' @param ... additional arguments for plot
#' @importFrom grDevices devAskNewPage rgb
#' @importFrom graphics axis lines plot
#' @method plot sucra
#' @importFrom ggplot2 ggplot geom_line aes labs
#' @export
"plot.sucra" <- function(x, legend.position = "none", ...) {
	nT <- length(x$SUCRA)
	cumeffectiveness <- apply(x$rankProb, 2, cumsum)
	names <- x$names
	oask <- devAskNewPage(TRUE)
	for (TRT in 1:nT) {
		Area=round(x$SUCRA[TRT], 3)
		ddd <- data.frame(trt = names, cdf = cumeffectiveness[,TRT], pdf = x$rankProb[,TRT])
		bb <- ggplot(ddd, aes(x = trt)) +
			 geom_line(aes(y = cdf), color = "#eab159", size = 1) +
			 geom_line(aes(y = pdf), color = rgb(0, 157, 114, maxColorValue = 255), linetype = "twodash", size = 1) +
			 theme_gray() + theme(legend.position = legend.position) +
			 ylab("SUCRA") + labs(title = paste0("Trt (", names[TRT], "): ", Area))
		print(bb)
		# plot(1:nT,type="n",ylim=c(0,1), xaxt="n",
		#    xlab=paste("Rank of",as.character(names[TRT])), 
		#    ylab="Probability", main = paste0("SUCRA = ", Area), ...)
		# axis(1,at=1:nT,labels=names)
		# lines(lwd=2,lty=1,c(1,c(1:c(nT)),nT), 
		#     cumeffectiveness[c(1,1:c(nT),c(nT)),TRT], col=rgb(0, 157, 114, maxColorValue = 255))
		# lines(lwd=2,lty=3,col="#eab159",c(1,c(1:c(nT)),nT), 
		#     x$rankProb[c(1,1:c(nT),c(nT)),TRT])

	}
	on.exit(devAskNewPage(oask))
}