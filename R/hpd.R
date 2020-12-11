#' get the highest posterior density (HPD) interval
#' @param object the output model from fitting a (network) meta analysis/regression model
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the probability which the HPD interval will cover
#' @param HPD a logical value indicating whether HPD or equal-tailed credible interval should be computed; by default, TRUE
#' @references
#' Chen, M. H., & Shao, Q. M. (1999). Monte Carlo estimation of Bayesian credible and HPD intervals. *Journal of Computational and Graphical Statistics*, **8(1)**, 69-92.
#' @md
#' @return dataframe containing HPD intervals for the parameters
#' @export

"hpd" <- function(object, parm, level = 0.95, HPD = TRUE) {
    UseMethod("hpd", object)
} 
