#' get surface under the cumulative ranking curve (SUCRA)
#' @param object the output model from fitting a network meta analysis/regression model
#' @return a matrix containing SUCRA
#' @export
"sucra" <- function(object) {
	UseMethod("sucra", object)
} 