#' plot the surface under the cumulative ranking curve (SUCRA)
#' @param x the output model from fitting a network meta analysis/regression model
#' @param legend.position the position of the legend that will be passed onto ggplot
#' @param ... additional arguments for plot
#' @importFrom grDevices devAskNewPage rgb
#' @importFrom graphics axis lines plot
#' @method plot sucra
#' @importFrom ggplot2 ggplot geom_line aes labs theme_gray theme ylab
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
	}
	on.exit(devAskNewPage(oask))
}