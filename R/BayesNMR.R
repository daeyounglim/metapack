#' Fit Bayesian Network Meta-Regression Hierarchical Models Using Heavy-Tailed Multivariate Random Effects with Covariate-Dependent Variances
#'
#' This is a function the fits the model introduced in *Bayesian Network Meta-Regression Models Using Heavy-Tailed Multivariate Random Effects with Covariate-Dependent Variances (under revision)*. The first seven arguments are required except `ZCovariate`. If not provided, `ZCovariate` will be assigned a vector of ones, `rep(1, length(Outcome))`. `ZCovariate` is the centerpiece of the modeling of variances and the heavy-tailed random effects distribution. 
#' @author Daeyoung Lim, \email{daeyoung.lim@uconn.edu}
#' @param Outcome the aggregate mean of the responses for each arm of every study.
#' @param SD the standard deviation of the responses for each arm of every study.
#' @param XCovariate the aggregate covariates for the fixed effects.
#' @param ZCovariate the aggregate covariates associated with the variance of the random effects.
#' @param Trial the study/trial identifiers. The elements within will be coerced to consecutive integers.
#' @param Treat the treatment identifiers for trial arm. This is equivalent to the arm labels in each study. The elements within will be coerced to consecutive integers
#' @param Npt the number of observations/participants for a unique `(t,k)`, or each arm of every trial.
#' @param prior (Optional) a list of hyperparameters. The hyperparameters include `df`, `c01`, `c02`, `a4`, `b4`, `a5`, and `b5`. `df` indicates the degrees of freedom whose value is 20. The hyperparameters `a*` and `b*` will take effect only if `sample_df=TRUE`. See `control`.
#' @param mcmc (Optional) a list of MCMC specification. `ndiscard` is the number of burn-in iterations. `nskip` configures the thinning of the MCMC. For instance, if `nskip=5`, `bayes.nmr` will save the posterior sample every 5 iterations. `nkeep` is the size of the posterior sample. The total number of iterations will be `ndiscard + nskip * nkeep`.
#' @param control (Optional) a list of parameters for [the Metropolis-Hastings algorithm](https://en.wikipedia.org/wiki/Metropolis-Hastings_algorithm). `lambda`, `phi`, and `Rho` are sampled through the localized Metropolis algorithm. `*_stepsize` with the asterisk replaced with one of the names above specifies the stepsize for determining the sample evaluation points in the localized Metropolis algorithm. `sample_Rho` can be set to `FALSE` to suppress the sampling of `Rho`. When `sample_Rho` is `FALSE`, `Rho` will be fixed using the value given by the `init` argument, which defaults to an equicorrelation matrix of \eqn{0.5\boldsymbol{I}+0.5\boldsymbol{1}\boldsymbol{1}^\prime}{0.5*I + 0.5*11'} where \eqn{\boldsymbol{1}}{1} is the vector of ones. When `sample_df` is `TRUE`, `df` will be sampled.
#' @param scale_x (Optional) a logical variable indicating whether `XCovariate` should be scaled/standardized. The effect of setting this to `TRUE` is not limited to merely standardizing `XCovariate`. The following generic functions will scale the posterior sample of `theta` back to its original unit: `plot`, `fitted`, `summary`, and `print`. That is `theta[j] <- theta[j] / sd(XCovariate[,j])`. 
#' @param init (Optional) a list of initial values for the parameters to be sampled: `theta`, `phi`, `sig2`, and `Rho`.
#' @param Treat_order (Optional) a vector of unique treatments to be used for renumbering the `Treat` vector. The first element will be assigned treatment zero, potentially indicating placebo. If not provided, the numbering will default to an alphabetical/numerical order.
#' @param Trial_order (Optional) a vector unique trials. The first element will be assigned trial zero. If not provided, the numbering will default to an alphabetical/numerical order.
#' @param verbose (Optional) a logical value indicating whether to print the progress bar during the MCMC sampling.
#' @return `bayes.nmr` returns an object of class `"bayesnmr"`. The functions `summary` or `print` are used to obtain and print a summary of the results. The generic accessor function `fitted` extracts the posterior mean, posterior standard deviation, and the interval estimates of the value returned by `bayes.nmr`.
#' 
#' An object of class `bayes.nmr` is a list containing the following components:
#' 
#' + `Outcome` - the aggregate response used in the function call.
#' + `SD` - the standard deviation used in the function call.
#' + `Npt` - the number of participants for `(t,k)` used in the function call.
#' + `XCovariate` - the aggregate design matrix for fixed effects used in the function call. Depending on `scale_x`, this may differ from the matrix provided at function call.
#' + `ZCovariate` - the aggregate design matrix for random effects. `bayes.nmr` will assign `rep(1, length(Outcome))` if it was not provided at function call.
#' + `Trial` - the *renumbered* trial indicators. Depending on `Trial_order`, it may differ from the vector provided at function call.
#' + `Treat` - the *renumbered* treatment indicators. Depending on `Treat_order`, it may differ from the vector provided at function call.
#' + `TrtLabels` - the vector of treatment labels corresponding to the renumbered `Treat`. This is equivalent to `Treat_order` if it was given at function call.
#' + `TrialLabels` - the vector of trial labels corresponding to the renumbered `Trial`. This is equivalent to `Trial_order` if it was given at function call.
#' + `K` - the total number of trials.
#' + `nT` - the total number of treatments.
#' + `scale_x` - a Boolean indicating whether `XCovariate` has been scaled/standardized.
#' + `prior` - the list of hyperparameters used in the function call.
#' + `control` - the list of tuning parameters used for MCMC in the function call.
#' + `mcmctime` - the elapsed time for the MCMC algorithm in the function call. This does not include all the other preprocessing and post-processing outside of MCMC.
#' + `mcmc` - the list of MCMC specification used in the function call.
#' + `mcmc.draws` - the list containing the MCMC draws. The posterior sample will be accessible here.
#' 
#' @references 
#' Li, H., Chen, M. H., Ibrahim, J. G., Kim, S., Shah, A. K., Lin, J., & Tershakovec, A. M. (2019). Bayesian inference for network meta-regression using multivariate random effects with applications to cholesterol lowering drugs. *Biostatistics*, **20(3)**, 499-516.
#' 
#' Li, H., Lim, D., Chen, M. H., Ibrahim, J. G., Kim, S., Shah, A. K., Lin, J. (2021). Bayesian network meta-regression hierarchical models using heavy-tailed multivariate random effects with covariate-dependent variances. Submitted.
#' 
#' @examples
#' library(metapack)
#' data(TNM)
#' groupInfo <- list(c("PBO"), c("R"))
#' nz <- length(groupInfo)
#' ns <- nrow(TNM)
#' XCovariate <- model.matrix(~ 0 + bldlc + bhdlc + btg + age +
#'  white + male + bmi + potencymed + potencyhigh + durat, data = TNM)
#' XCovariate <- scale(XCovariate, center = TRUE, scale = FALSE)
#' ZCovariate <- matrix(0, ns, nz)
#' for (j in 1:length(groupInfo)) {
#'     for (i in 1:ns) {
#'         if (TNM$Treat[i] %in% groupInfo[[j]]) {
#'             ZCovariate[i, j] <- 1
#'         }
#'     }
#' }
#' addz <- scale(cbind(TNM$bldlc, TNM$btg), center=TRUE, scale=TRUE)
#' ZCovariate <- cbind(1, ZCovariate, addz)
#' theta_init <- c(0.05113, -1.38866, 1.09817, -0.85855, -1.12056, -1.14133,
#'              -0.22435, 3.63453, -2.09322, 1.07858, 0.80566, -40.76753,
#'              -45.07127, -28.27232, -44.14054, -28.13203, -19.19989,
#'              -47.21824, -51.31234, -48.46266, -47.71443)
#' set.seed(2797542)
#' fit <- bayes.nmr(TNM$ptg, TNM$sdtg, XCovariate, ZCovariate, TNM$Trial,
#'     TNM$Treat, TNM$Npt, prior = list(c01 = 1.0e05, c02 = 4, df = 3),
#'     mcmc = list(ndiscard = 1, nskip = 1, nkeep = 1),
#'     init = list(theta = theta_init),
#'     Treat_order = c("PBO", "S", "A", "L", "R", "P", "E", "SE",
#'          "AE", "LE", "PE"),
#'     scale_x = TRUE, verbose = FALSE)
#' @importFrom stats model.matrix
#' @importFrom methods is
#' @md
#' @export
bayes.nmr <- function(Outcome, SD, XCovariate, ZCovariate, Trial, Treat, Npt, prior = list(), mcmc = list(), scale_x = FALSE, control = list(), init = list(), Treat_order = NULL, Trial_order = NULL, verbose = FALSE) {
  if (!is(Outcome, "vector")) {
    tmp <- try(Outcome <- as.vector(Outcome))
    if (is(tmp, "try-error")) {
      stop("Outcome must be a vector or able to be coerced to a vector")
    }
  }
  if (!is(SD, "vector")) {
    tmp <- try(SD <- as.vector(SD))
    if (is(tmp, "try-error")) {
      stop("SD must be a vector or able to be coerced to a vector")
    }
  }
  if (!is(XCovariate, "matrix")) {
    tmp <- try(XCovariate <- model.matrix(~ 0 + ., data = XCovariate), silent = TRUE)
    if (is(tmp, "try-error")) {
      stop("XCovariate must be a matrix or able to be coerced to a matrix")
    }
  }
  if (missing(ZCovariate)) {
    ZCovariate <- rep(1, nrow(XCovariate))
  } else {
    if (!is(ZCovariate, "matrix")) {
      tmp <- try(ZCovariate <- model.matrix(~ 0 + ., data = ZCovariate), silent = TRUE)
      if (is(tmp, "try-error")) {
        stop("ZCovariate must be a matrix or able to be coerced to a matrix")
      }
    }
  }
  if (!is(Treat, "vector")) {
    tmp <- try(Treat <- as.vector(Treat))
    if (is(tmp, "try-error")) {
      stop("Treat must be a vector or able to be coerced to a vector")
    }
  }
  if (!is(Trial, "vector")) {
    tmp <- try(Trial <- as.vector(Trial))
    if (is(tmp, "try-error")) {
      stop("Trial must be a vector or able to be coerced to a vector")
    }
  }
  if (!is(Npt, "numeric")) {
    tmp <- try(Npt <- as.numeric(Npt))
    if (is(tmp, "try-error")) {
      stop("Npt must be numeric or able to be coerced to numeric")
    }
  }
  if (any(is.na(Outcome)) | any(is.na(XCovariate)) | any(is.na(ZCovariate)) | any(is.na(Treat)) | any(is.na(Trial)) | any(is.na(Npt))) {
    stop("Missing data (NA) detected. Handle missing data (e.g., delete missings, delete variables, imputation) before passing it as an argument")
  }

  mcvals <- list(ndiscard = 5000L, nskip = 1L, nkeep = 20000L)
  mcvals[names(mcmc)] <- mcmc
  ndiscard <- mcvals$ndiscard
  nskip <- mcvals$nskip
  if (nskip < 1) {
    stop("The thinning cannot be smaller than 1.")
  }
  nkeep <- mcvals$nkeep

  priorvals <- list(df = 20, c01 = 1.0e05, c02 = 4, a4 = 1, b4 = 0.1, a5 = 0.1, b5 = 0.1)
  priorvals[names(prior)] <- prior
  df <- priorvals$df
  c01 <- priorvals$c01
  c02 <- priorvals$c02
  a4 <- priorvals$a4
  b4 <- priorvals$b4
  a5 <- priorvals$a5
  b5 <- priorvals$b5

  if (is.null(Treat_order)) {
    Treat.order <- sort(unique(Treat))
  } else {
    Treat.order <- Treat_order
  }
  Treat.n <- relabel.vec(Treat, Treat.order) - 1 # relabel the treatment numbers

  if (is.null(Trial_order)) {
    Trial.order <- sort(unique(Trial))
  } else {
    Trial.order <- Trial_order
  }
  Trial.n <- relabel.vec(Trial, Trial.order) - 1 # relabel the trial numbers

  nx <- ncol(XCovariate)
  nz <- ncol(ZCovariate)
  ns <- length(Outcome)
  K <- length(unique(Trial))
  nT <- length(unique(Treat))
  xn <- if (!is.null(colnames(XCovariate))) colnames(XCovariate) else paste0("beta", 1:ncol(XCovariate))

  if (scale_x) {
    XCovariate_ <- scale(XCovariate, center = FALSE, scale = TRUE)
  } else {
    XCovariate_ <- XCovariate
  }

  nt <- nx + nT
  XSX <- matrix(0, nt, nt)
  XSy <- numeric(nt)
  for (k in 1:K) {
    idx <- which(Trial.n == k-1)
    Xk <- XCovariate_[idx,]
    idx_l = length(idx)
    Ek <- matrix(0, nT, idx_l)
    iarm_k = Treat.n[idx]
    for (j in 1:idx_l) {
      Ek[iarm_k[j]+1,j] = 1
    }
    Xstar <- cbind(Xk, t(Ek))
    Skinv <- diag(Npt[idx] / SD[idx]^2, nrow=idx_l)
    XSX <- XSX + crossprod(Xstar, Skinv %*% Xstar)
    XSy <- XSy + crossprod(Xstar, Skinv %*% Outcome[idx])
  }
  thetahat <- solve(XSX + 0.1 * diag(1, nrow=nt), XSy)

  zn <- if (!is.null(Treat_order)) Treat_order else paste0("gam", 1:nT)
  init_final <- list(theta = thetahat, phi = numeric(nz), sig2 = SD^2, Rho = diag(1, nrow = nT))
  init_final$Rho[upper.tri(init_final$Rho)] <- 0.5
  init_final$Rho[lower.tri(init_final$Rho)] <- 0.5
  init_final[names(init)] <- init
  Rho_init <- init_final$Rho
  if (any(eigen(Rho_init, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
    stop("The initial value for Omega is not positive definite")
  }

  ctrl <- list(
    lambda_stepsize = 0.5,
    phi_stepsize = 0.5,
    Rho_stepsize = 0.2,
    sample_Rho = TRUE,
    sample_df = FALSE
  )
  ctrl[names(control)] <- control
  lambda_stepsize <- ctrl$lambda_stepsize
  phi_stepsize <- ctrl$phi_stepsize
  Rho_stepsize <- ctrl$Rho_stepsize
  sample_Rho <- ctrl$sample_Rho
  sample_df <- ctrl$sample_df

  if (is.infinite(df) && sample_df) {
    stop("The degrees of freedom cannot be sampled for a normal random-effects model")
  }


  mcmctime <- system.time({
    fout <- .Call(
      `_metapack_BayesNMR`,
      as.double(Outcome),
      as.double(SD),
      as.matrix(XCovariate_),
      as.matrix(ZCovariate),
      as.integer(Trial.n),
      as.integer(Treat.n),
      as.double(Npt),
      as.double(df),
      as.double(1 / c01),
      as.double(1 / c02),
      as.double(a4),
      as.double(b4),
      as.double(a5),
      as.double(b5),
      as.integer(K),
      as.integer(nT),
      as.integer(ndiscard),
      as.integer(nskip),
      as.integer(nkeep),
      as.logical(verbose),
      as.double(init_final$theta),
      as.double(init_final$phi),
      as.double(init_final$sig2),
      as.matrix(Rho_init),
      as.double(lambda_stepsize),
      as.double(phi_stepsize),
      as.double(Rho_stepsize),
      as.logical(sample_Rho),
      as.logical(sample_df)
    )
  })
  rownames(fout$theta) <- c(xn, zn)
  out <- list(
    Outcome = Outcome,
    SD = SD,
    Npt = Npt,
    XCovariate = XCovariate_,
    ZCovariate = ZCovariate,
    Treat = Treat.n,
    Trial = Trial.n,
    TrtLabels = Treat.order,
    TrialLabels = Trial.order,
    K = K,
    nT = nT,
    scale_x = scale_x,
    prior = priorvals,
    control = ctrl,
    mcmctime = mcmctime,
    mcmc = mcvals,
    mcmc.draws = fout
  )
  class(out) <- "bayesnmr"
  return(out)
}
