#' get the highest posterior density (HPD) interval
#' @param object the output model from fitting a (network) meta analysis/regression model
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the probability which the HPD interval will cover
#' @param HPD a logical value indicating whether HPD or equal-tailed credible interval should be computed; by default, TRUE
#' @return dataframe containing HPD intervals for the parameters
#' @export

"hpd" <- function(object, parm, level = 0.95, HPD = TRUE) {
    UseMethod("hpd", object)
} 