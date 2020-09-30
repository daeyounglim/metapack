#' Fit Bayesian Inference for Multivariate Meta-Regression With a Partially Observed Within-Study Sample Covariance Matrix
#'
#' This is a function for running the Markov chain Monte Carlo algorithm for the BMVMR_POCov Model. The first six arguments are required.
#' fmodel can be one of 5 numbers: 1, 2, 3, 4, and 5. The first model, fmodel = 1 denoted by M1, indicates that the \eqn{\Sigma_{tk}}{Sig.tk}
#' are diagonal matrices with zero covariances. M2 indicates that \eqn{\Sigma_{tk}}{Sig.tk} are all equivalent but allowed to be full symmetric
#' positive definite. M3 is where \eqn{\Sigma_{tk}}{Sig.tk} are allowed to differ across treatments, i.e., \eqn{\Sigma_{tk}=\Sigma_t}{Sig.tk=Sig.t}.
#' M4 assumes thata the correlation matrix, \eqn{\rho}{Rho}, is identical for all trials/treatments, but the variances are allowed to vary.
#' Finally, M5 assumes a hierarchical model where \eqn{(\Sigma_{tk}\mid \Sigma)}{(Sig.tk|Sig)} follows an inverse-Wishart distribution with fixed
#' degrees of freedom and scale matrix \eqn{\Sigma}{Sig}. \eqn{\Sigma}{Sig} then follows another inverse-Wishart distribution with fixed parameters.
#' @author Daeyoung Lim, \email{daeyoung.lim@uconn.edu}
#' @param Outcome aggregate mean of the responses for each arm of each study
#' @param SD standard deviation of the responses for each arm of each study
#' @param XCovariate aggregate covariates for the fixed effects
#' @param WCovariate aggregate covariates for the random effects
#' @param Treat treatment applied; equivalent to the arm number of each study; the number of unique treatments must be equal across trials
#' @param Trial trial identifier
#' @param Npt number of observations per trial
#' @param fmodel the model number; defaults to M1
#' @param prior list of hyperparameters; when not given, algorithm will run in default setting
#' @param mcmc list of MCMC-related parameters: number of burn-ins (ndiscard), number of thinning(nskip), and posterior sample size (nkeep)
#' @param control list of parameters for localized Metropolis algorithm: the step sizes for R, Rho, delta, and Delta (R_stepsize, Rho_stepsize, delta_stepsize); If not provided, default to 0.02, 0.02, and 0.2, respectively; sample_Rho is a logical value, by default TRUE; if sample_Rho=FALSE, MCMC sampling of Rho is suppressed in fmodel=3
#' @param init initial values for the parameters. Dimensions must be conformant.
#' @param Treat_order a vector of unique treatments to be used for renumbering the 'Treat' vector; the first element will be assigned treatment zero, potentially indicating placebo; if not provided, the numbering will default to an alphabetical/numerical order
#' @param Trial_order a vector of unique trials; the first element will be assigned trial zero; if not provided, the numbering will default to an alphabetical/numerical order
#' @param verbose logical variable for printing progress bar. Default to FALSE.
#' @return a dataframe with input arguments, posterior samples, Metropolis algorithm acceptance rates, etc
#' @examples
#' \dontrun{
#' data("cholesterol")
#' Outcome <- cbind(cholesterol$pldlc, cholesterol$phdlc, cholesterol$ptg)
#' SD <- cbind(cholesterol$sdldl, cholesterol$sdhdl, cholesterol$sdtg)
#' Trial <- cholesterol$Trial
#' Treat <- cholesterol$trt
#' Npt <- cholesterol$Npt
#' XCovariate <- cbind(cholesterol$bldlc, cholesterol$bhdlc,
#' 		cholesterol$btg, cholesterol$age, cholesterol$durat, cholesterol$white,
#' 		cholesterol$male, cholesterol$dm)
#' WCovariate <- cbind(1 - cholesterol$onstat, cholesterol$trt * (1 - cholesterol$onstat),
#' 		cholesterol$onstat, cholesterol$trt * cholesterol$onstat)
#' 
#' fmodel <- 3
#' fit <- bayes.parobs(Outcome, SD, scale(XCovariate, scale = TRUE, center = TRUE),
#' 			WCovariate, Treat, Trial, Npt, fmodel,
#'   		mcmc = list(ndiscard = 100000, nskip = 1, nkeep = 20000),
#' 			control = list(delta_stepsize = 0.1,
#' 			rho_stepsize = 0.05, R_stepsize = 0.05), verbose = TRUE
#' )
#' }
#' @importFrom stats model.matrix
#' @importFrom methods is
#' @export
bayes.parobs <- function(Outcome, SD, XCovariate, WCovariate, Treat, Trial, Npt, fmodel = 1, prior = list(), mcmc = list(), control = list(), init = list(), Treat_order = NULL, Trial_order = NULL, verbose = FALSE) {
  if (!is(Outcome, "matrix")) {
    tmp <- try(Outcome <- model.matrix(~ 0 + ., data = Outcome), silent = TRUE)
    if (is(tmp, "try-error")) {
      stop("Outcome must be a matrix or able to be coerced to a matrix")
    }
  }
  if (!is(SD, "matrix")) {
    tmp <- try(SD <- model.matrix(~ 0 + ., data = SD), silent = TRUE)
    if (is(tmp, "try-error")) {
      stop("SD must be a matrix or able to be coerced to a matrix")
    }
  }
  if (!is(XCovariate, "matrix")) {
    tmp <- try(XCovariate <- model.matrix(~ 0 + ., data = XCovariate), silent = TRUE)
    if (is(tmp, "try-error")) {
      stop("XCovariate must be a matrix or able to be coerced to a matrix")
    }
  }
  if (!is(WCovariate, "matrix")) {
    tmp <- try(WCovariate <- model.matrix(~ 0 + ., data = WCovariate), silent = TRUE)
    if (is(tmp, "try-error")) {
      stop("WCovariate must be a matrix or able to be coerced to a matrix")
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
  if (any(is.na(Outcome)) | any(is.na(SD)) | any(is.na(XCovariate)) | any(is.na(WCovariate)) | any(is.na(Treat)) | any(is.na(Trial))) {
    stop("Missing data (NA) detected. Handle missing data (e.g., delete missings, delete variables, imputation) before passing it as an argument")
  }

  mcvals <- list(ndiscard = 5000L, nskip = 1L, nkeep = 20000L)
  mcvals[names(mcmc)] <- mcmc
  ndiscard <- mcvals$ndiscard
  nskip <- mcvals$nskip
  nkeep <- mcvals$nkeep

  ctrl <- list(R_stepsize = 0.02, Rho_stepsize = 0.02, delta_stepsize = 0.2, sample_Rho = TRUE)
  ctrl[names(control)] <- control
  R_stepsize <- ctrl$R_stepsize
  Rho_stepsize <- ctrl$Rho_stepsize
  delta_stepsize <- ctrl$delta_stepsize
  sample_Rho <- ctrl$sample_Rho # for fmodel = 3

  J <- ncol(Outcome)
  nw <- ncol(WCovariate)
  priorvals <- list(c0 = 1.0e05, dj0 = 0.1 + nw, d0 = 0.1 + J, s0 = 0.1, a0 = 0.1, b0 = 0.1, Omega0 = diag(10, nw), Sigma0 = diag(10, J), nu0 = 10)
  priorvals[names(prior)] <- prior
  a0 <- priorvals$a0
  b0 <- priorvals$b0
  c0 <- priorvals$c0
  dj0 <- priorvals$dj0
  d0 <- priorvals$d0
  s0 <- priorvals$s0
  Omega0 <- priorvals$Omega0
  Sigma0 <- priorvals$Sigma0
  nu0 <- priorvals$nu0

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

  K <- length(unique(Trial))
  T <- length(unique(Treat))

  xcols <- ncol(XCovariate)

  init_final <- list(
    theta = numeric((xcols + nw) * J),
    gamR = matrix(0, nw * J, K), Omega = diag(1, nrow = nw * J),
    Rho = diag(1, nrow = J)
  )
  init_final[names(init)] <- init
  theta_init <- init_final$theta
  gamR_init <- init_final$gamR
  Omega_init <- init_final$Omega
  Rho_init <- init_final$Rho

  if (any(eigen(Omega_init, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
    stop("The initial value for Omega is not positive definite")
  }
  if (any(eigen(Rho_init, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
    stop("The initial value for Omega is not positive definite")
  }

  mcmctime <- system.time({
    if (fmodel == 1) {
      fout <- .Call(
        `_metapack_fmodel1`,
        as.matrix(Outcome),
        as.matrix(SD),
        as.matrix(XCovariate),
        as.matrix(WCovariate),
        as.integer(Treat.n),
        as.integer(Trial.n),
        as.double(Npt),
        as.double(c0),
        as.double(dj0),
        as.double(a0),
        as.double(b0),
        as.matrix(Omega0),
        as.integer(K),
        as.integer(T),
        as.integer(ndiscard),
        as.integer(nskip),
        as.integer(nkeep),
        as.double(theta_init),
        as.matrix(gamR_init),
        as.matrix(Omega_init),
        as.logical(verbose)
      )
    } else if (fmodel == 2) {
      fout <- .Call(
        `_metapack_fmodel2`,
        as.matrix(Outcome),
        as.matrix(SD),
        as.matrix(XCovariate),
        as.matrix(WCovariate),
        as.integer(Treat.n),
        as.integer(Trial.n),
        as.double(Npt),
        as.double(c0),
        as.double(dj0),
        as.double(s0),
        as.matrix(Omega0),
        as.matrix(Sigma0),
        as.integer(K),
        as.integer(T),
        as.integer(ndiscard),
        as.integer(nskip),
        as.integer(nkeep),
        as.double(theta_init),
        as.matrix(gamR_init),
        as.matrix(Omega_init),
        as.double(R_stepsize),
        as.logical(verbose)
      )
	} else if (fmodel == 3) {
		fout <- .Call(
        `_metapack_fmodel2p5`,
        as.matrix(Outcome),
        as.matrix(SD),
        as.matrix(XCovariate),
        as.matrix(WCovariate),
        as.integer(Treat.n),
        as.integer(Trial.n),
        as.double(Npt),
        as.double(c0),
        as.double(dj0),
        as.double(s0),
        as.matrix(Omega0),
        as.matrix(Sigma0),
        as.integer(K),
        as.integer(T),
        as.integer(ndiscard),
        as.integer(nskip),
        as.integer(nkeep),
        as.double(theta_init),
        as.matrix(gamR_init),
        as.matrix(Omega_init),
        as.double(R_stepsize),
        as.logical(verbose)
      )
    } else if (fmodel == 4) {
      fout <- .Call(
        `_metapack_fmodel3`,
        as.matrix(Outcome),
        as.matrix(SD),
        as.matrix(XCovariate),
        as.matrix(WCovariate),
        as.integer(Treat.n),
        as.integer(Trial.n),
        as.double(Npt),
        as.double(c0),
        as.double(dj0),
        as.double(a0),
        as.double(b0),
        as.matrix(Omega0),
        as.integer(K),
        as.integer(T),
        as.integer(ndiscard),
        as.integer(nskip),
        as.integer(nkeep),
        as.double(delta_stepsize),
        as.double(Rho_stepsize),
        as.double(R_stepsize),
        as.double(theta_init),
        as.matrix(gamR_init),
        as.matrix(Omega_init),
        as.matrix(Rho_init),
        as.logical(sample_Rho),
        as.logical(verbose)
      )
    } else if (fmodel == 5) {
      fout <- .Call(
        `_metapack_fmodel4`,
        as.matrix(Outcome),
        as.matrix(SD),
        as.matrix(XCovariate),
        as.matrix(WCovariate),
        as.integer(Treat.n),
        as.integer(Trial.n),
        as.double(Npt),
        as.double(c0),
        as.double(dj0),
        as.double(d0),
        as.double(nu0),
        as.matrix(Sigma0),
        as.matrix(Omega0),
        as.integer(K),
        as.integer(T),
        as.integer(ndiscard),
        as.integer(nskip),
        as.integer(nkeep),
        as.double(delta_stepsize),
        as.double(Rho_stepsize),
        as.double(R_stepsize),
        as.double(theta_init),
        as.matrix(gamR_init),
        as.matrix(Omega_init),
        as.logical(verbose)
      )
    } else {
      stop("`fmodel` is invalid. Please pick from {1, 2, 3, 4, 5}.")
    }
  })
  if (!is.null(colnames(XCovariate)) && !is.null(colnames(WCovariate))) {
    rownames(fout$theta) <- c(rep(colnames(XCovariate), J), rep(colnames(WCovariate), J))
  }
  out <- list(
    Outcome = Outcome,
    SD = SD,
    Npt = Npt,
    XCovariate = XCovariate,
    WCovariate = WCovariate,
    Treat = Treat.n,
    Trial = Trial.n,
    TrtLabels = Treat.order,
    TrialLabels = Trial.order,
    K = K,
    T = T,
    fmodel = fmodel,
    prior = priorvals,
    mcmctime = mcmctime,
    mcmc = mcvals,
    mcmc.draws = fout
  )
  class(out) <- "bayes.parobs"
  out
}
