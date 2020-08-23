#' get goodness of fit: DIC or LPML
#' @param object the output model from fitting a meta analysis/regression model
#' @param type the type of goodness of fit to compute; DIC or LPML
#' @param verbose FALSE by default; If TRUE, then progress bar will appear
#' @return dataframe containing the goodness of fit measure
#' @export
"gof" <- function(object, type = "lpml", verbose=FALSE) {
	UseMethod("gof", object)
} 