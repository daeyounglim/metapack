#' get the highest posterior density (HPD) interval
#' @param object the output model from fitting a meta analysis/regression model
#' @param prob the probability which the HPD interval will cover
#' @importFrom coda mcmc HPDinterval
#' @return dataframe containing HPD intervals for the parameters
#' @export
"hpd" <- function(object, prob) {
	UseMethod("hpd", object)
} 