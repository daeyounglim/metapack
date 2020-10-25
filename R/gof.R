#' get goodness of fit: DIC or LPML
#' @param object the output model from fitting a meta analysis/regression model
#' @param type the type of goodness of fit to compute; DIC or LPML
#' @param ncores the number of CPU cores to use for parallel processing; it must not exceed the number of existing cores
#' @param verbose FALSE by default; If TRUE, then progress bar will appear
#' @return dataframe containing the goodness of fit measure
#' @export
"gof" <- function(object, type = "lpml", ncores=NULL, verbose=FALSE) {
	UseMethod("gof", object)
} 