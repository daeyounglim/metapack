#' Fit Bayesian Network Meta-Regression Hierarchical Models Using Heavy-Tailed Multivariate Random Effects with Covariate-Dependent Variances
#'
#' This is a function for running the Markov chain Monte Carlo algorithm for the BNMHHtMRe Model. The first seven arguments are required.
#' @author Daeyoung Lim, \email{daeyoung.lim@uconn.edu}
#' @param Outcome aggregate mean of the responses for each arm of each study
#' @param SD standard deviation of the responses for each arm of each study
#' @param Covariate aggregate covariates for the mean component
#' @param Trial study/trial identifiers; will be coerced to consecutive integers
#' @param Treat treatment identifiers for the corresponding trial arm; equivalent to the arm number of each study; will be coerced to consecutive integers
#' @param Npt number of observations/participants per trial
#' @param groupInfo list of grouping information; the control(baseline) group must start from 0; the aggregate covariates 'z' explaining the variance of the random effect of the t-th treatment will be construct based on this grouping information
#' @param prior list of hyperparameters when not given, algorithm will run in default setting.
#' @param mcmc list of MCMC-related parameters: number of burn-ins (ndiscard), number of thinning(nskip), and posterior sample size (nkeep)
#' @param add.z additional covariates other than the grouping vectors that should be column-concatenated to 'Z'. This should have the same number of rows as 'Outcome', and 'Covariate'
#' @param verbose logical variable for printing progress bar. Default to FALSE.
#' @param control list of parameters for localized Metropolis algorithm: the step sizes for lambda, phi, and Rho (lambda_stepsize, phi_stepsize, Rho_stepsize); All default to 0.5 except Rho_stepsize if not set; Rho_stepsize defaults to 0.2; sample_Rho is a logical value, by default TRUE; if sample_Rho=FALSE, MCMC sampling of Rho is suppressed; if sample_df is set to TRUE, the degrees of freedom for the t random effects will be treated as unknown and sampled in the MCMC algorithm
#' @param Treat_order a vector of unique treatments to be used for renumbering the 'Treat' vector; the first element will be assigned treatment zero, potentially indicating placebo; if not provided, the numbering will default to an alphabetical/numerical order
#' @param Trial_order a vector of unique trials; the first element will be assigned trial zero; if not provided, the numbering will default to an alphabetical/numerical order
#' @param init initial values for theta (ns + nT dimensional) and phi. Dimensions must be conformant.
#' @return a dataframe with input arguments, posterior samples, Metropolis algorithm acceptance rates, etc
#' @examples
#' \dontrun{
#' data(TNM)
#' groupInfo <- list(c(0, 1), c(2, 3), c(4)) # define the variance structure
#' x <- TNM[, 6:15]
#' x <- scale(x, center = TRUE, scale = TRUE)
#' beta_init <- c(0.05113, -1.38866, 1.09817, -0.85855, -1.12056, -1.14133, -0.22435, 3.63453,
#'              -2.09322, 1.07858, 0.80566, -40.76753, -45.07127, -28.27232, -44.14054,
#'              -28.13203, -19.19989, -47.21824, -51.31234, -48.46266, -47.71443)
#' fit <- bayes.nmr(TNM$ptg, TNM$sdtg, x, TNM$Trial, TNM$Treat, TNM$Npt, groupInfo,
#'   add.z = cbind(TNM$bldlc, TNM$btg),
#'   prior = list(c01 = 1.0e05, c02 = 4, df = 3),
#'   mcmc = list(ndiscard = 2500, nskip = 1, nkeep = 10000),
#'   init = list(beta = beta_init),
#'   Treat_order = c("PBO", "S", "A", "L", "R", "P", "E", "SE", "AE", "LE", "PE"),
#'   verbose = TRUE
#' )
#' }
#' @importFrom stats model.matrix
#' @importFrom methods is
#' @export
bayes.nmr <- function(Outcome, SD, Covariate, Trial, Treat, Npt, groupInfo, prior = list(), mcmc = list(), add.z = NULL, control = list(), init = list(), Treat_order = NULL, Trial_order = NULL, verbose = FALSE) {
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
  if (!is(Covariate, "matrix")) {
    tmp <- try(Covariate <- model.matrix(~ 0 + ., data = Covariate), silent = TRUE)
    if (is(tmp, "try-error")) {
      stop("Covariate must be a matrix or able to be coerced to a matrix")
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
  if (any(is.na(Outcome)) | any(is.na(Covariate)) | any(is.na(Treat)) | any(is.na(Trial)) | any(is.na(Npt))) {
    stop("Missing data (NA) detected. Handle missing data (e.g., delete missings, delete variables, imputation) before passing it as an argument")
  }

  mcvals <- list(ndiscard = 5000L, nskip = 1L, nkeep = 20000L)
  mcvals[names(mcmc)] <- mcmc
  ndiscard <- mcvals$ndiscard
  nskip <- mcvals$nskip
  nkeep <- mcvals$nkeep
  if (length(groupInfo) == 0) {
    stop("groupInfo is missing.")
  }
  if (min(c(unlist(groupInfo))) != 0) {
    warning("Baseline treatment for groupInfo should start from 0.\nAssuming baseline is not ")
  }

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

  nx <- ncol(Covariate)
  nz <- length(groupInfo)
  ns <- length(Outcome)
  K <- length(unique(Trial))
  nT <- length(unique(Treat))
  z <- matrix(0, ns, nz)
  if (nz > 0) {
    for (j in 1:nz) {
      for (i in 1:ns) {
        if (Treat.n[i] %in% groupInfo[[j]]) {
          z[i, j] <- 1
        }
      }
    }
  }
  if (!is.null(add.z)) {
    z <- cbind(z, add.z)
  }
  init_final <- list(theta = numeric(nx + nT), phi = numeric(ncol(z)), sig2 = rep(1, ns), Rho = diag(1, nrow = nT))
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
    stop("Cannot sample degrees of freedom for a normal random effects model")
  }

  mcmctime <- system.time({
    fout <- .Call(
      `_metapack_BayesNMR`,
      as.double(Outcome),
      as.double(SD),
      as.matrix(Covariate),
      as.matrix(z),
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

  out <- list(
    Outcome = Outcome,
    SD = SD,
    Npt = Npt,
    Covariate = Covariate,
    z = z,
    Trial = Trial.n,
    Treat = Treat.n,
    TrtLabels = Treat.order,
    K = K,
    nT = nT,
    groupInfo = groupInfo,
    prior = priorvals,
    control = ctrl,
    mcmctime = mcmctime,
    mcmc = mcvals,
    mcmc.draws = fout
  )
  class(out) <- "bayesnmr"
  out
}
