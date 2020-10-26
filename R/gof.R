#' get goodness of fit: DIC or LPML
#' @param object the output model from fitting a meta analysis/regression model
#' @param type the type of goodness of fit to compute; DIC or LPML
#' @param verbose FALSE by default; If TRUE, then progress bar will appear
#' @param ncores the number of CPU cores to use for parallel processing; it must not exceed the number of existing cores
#' @param h the interval width for trapezoidal rule
#' @return dataframe containing the goodness of fit measure
#' @export
"gof" <- function(object, type = "lpml", verbose=FALSE, ncores=NULL, h = 0.5) {
	UseMethod("gof", object)
} 