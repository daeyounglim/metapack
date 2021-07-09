#' Wrapper function for worker functions: bayes.parobs, bayes.nmr
#' 
#' This is the one function to rule them all. All other worker functions will be subsumed by this function, so that users can forget about the implementation details.
#' @author Daeyoung Lim, \email{daeyoung.lim@uconn.edu}
#' @param formula an object of class of \link[Formula]{Formula}: a symbolic description of the meta-analytic model to be fitted. The list of models includes multivariate meta-regression and univariate network meta-regression. More will be added moving forward.
#' @param data a data frame, list, or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables in the model. If not found in `data`, the variables are taken from `environment(formula)`, typically the environment from which `synthesize` is called.
#' @param prior an optional object that contains the hyperparameter values for the model
#' @param control an optional object that contains the control tuning parameters for the Metropolis-Hastings algorithm
#' @import Formula
#' @export
'synthesize' <- function(formula, data, prior = list(), control = list(), ...) {
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m)]

    f <- Formula(formula)
    mf[[1]] <- as.name("model.frame")
    mf$formula <- f
    mf <- eval(mf, parent.frame())


    if (is.null(y <- model.response(mf))) {
        y <- model.part(f, data = mf, lhs = 1)
    }
    x <- model.matrix(f, data = mf, rhs = 1)
    z <- model.matrix(f, data = mf, rhs = 2)
    browser()
}