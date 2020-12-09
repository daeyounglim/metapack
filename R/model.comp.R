#' compute the model comparison measures: DIC or LPML
#' @param object the output model from fitting a meta analysis/regression model
#' @param type the type of model comparison measure to compute; DIC or LPML
#' @param verbose FALSE by default; If TRUE, then progress bar will appear
#' @param ncores the number of CPU cores to use for parallel processing; it must not exceed the number of existing cores
#' @return dataframe containing the compute the model comparison measures
#' @export
"model.comp" <- function(object, type = "lpml", verbose=FALSE, ncores=NULL) {
	UseMethod("model.comp", object)
} 