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
#' @param prior (Optional) list of hyperparameters; when not given, algorithm will run in default setting
#' @param mcmc (Optional) list of MCMC-related parameters: number of burn-ins (ndiscard), number of thinning(nskip), and posterior sample size (nkeep)
#' @param control (Optional) list of parameters for localized Metropolis algorithm: the step sizes for R, Rho, delta, and Delta (R_stepsize, Rho_stepsize, delta_stepsize); If not provided, default to 0.02, 0.02, and 0.2, respectively; sample_Rho is a logical value, by default TRUE; if sample_Rho=FALSE, MCMC sampling of Rho is suppressed in fmodel=3
#' @param init (Optional) initial values for the parameters. Dimensions must be conformant.
#' @param Treat_order (Optional) a vector of unique treatments to be used for renumbering the 'Treat' vector; the first element will be assigned treatment zero, potentially indicating placebo; if not provided, the numbering will default to an alphabetical/numerical order
#' @param Trial_order (Optional) a vector of unique trials; the first element will be assigned trial zero; if not provided, the numbering will default to an alphabetical/numerical order
#' @param group (Optional) a vector of binary group indicators; it must be binary indicating groups that have different random effects
#' @param group_order (Optional) a vector of unique group labels; the first element will be assigned zero; if not provided, the numbering will default to an alphabetical/numerical order
#' @param scale_x (Optional) a logical variable whether `XCovariate` should be scaled; if `TRUE`, `theta` will be scaled back to its original scale after posterior sampling
#' @param verbose (Optional) a logical variable for printing progress bar. Default to FALSE.
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
#' WCovariate <- cbind(1, cholesterol$trt)
#' 
#' fmodel <- 3
#' fit <- bayes.parobs(Outcome, SD, scale(XCovariate, scale = TRUE, center = TRUE),
#' 			WCovariate, Treat, Trial, Npt, fmodel,
#'   		mcmc = list(ndiscard = 100000, nskip = 1, nkeep = 20000),
#' 			control = list(delta_stepsize = 0.1,
#' 			rho_stepsize = 0.05, R_stepsize = 0.05),
#'      group = cholesterol$onstat,
#'      verbose = TRUE)
#' }
#' @importFrom stats model.matrix optim
#' @importFrom methods is
#' @importFrom Matrix nearPD
#' @md
#' @export
bayes.parobs <- function(Outcome, SD, XCovariate, WCovariate, Treat, Trial, Npt, fmodel = 1, prior = list(), mcmc = list(), control = list(), init = list(), Treat_order = NULL, Trial_order = NULL, group = NULL, group_order = NULL, scale_x = FALSE, verbose = FALSE) {
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

  nr <- 1
  second.exist <- FALSE
  if (!is.null(group)) {
    second.exist <- TRUE
    sl <- unique(group)
    if (length(sl) != 2) {
      stop("The second-line trial indicators must be binary.")
    }
    nr <- 2
    if (is.null(group_order)) {
      group.order <- sort(sl)
    } else {
      group.order <- group_order
    }
    group.n <- relabel.vec(group, group.order) - 1
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
  sample_Rho <- ctrl$sample_Rho # for fmodel = 4

  J <- ncol(Outcome)
  nw <-  ncol(WCovariate)
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
    theta = numeric((xcols + nw*nr) * J),
    gamR = matrix(0, nw * nr * J, K), Omega = diag(1, nrow = nw*nr*J),
    Rho = diag(1, nrow = J)
  )
  init_final[names(init)] <- init
  theta_init <- init_final$theta
  gamR_init <- init_final$gamR
  Omega_init <- init_final$Omega
  Rho_init <- init_final$Rho

  # if (fmodel == 4) {
  #   N <- nrow(Outcome)
  #   pars <- vecr(Rho_init)
  #   sumNpt <- sum(Npt)
  #   qq <- matrix(0, J, J)
  #   for (i in 1:N) {
  #     k <- Trial.n[i] + 1
  #     ntk <- Npt[i]
  #     x_i <- XCovariate[i,]
  #     w_i <- WCovariate[i,]
  #     pRR <- vecrinv(tanh(rep(0.5, (J*(J-1))/2)), as.integer(J))
  #     diag(pRR) <- 1
  #     RR <- pRho_to_Rho(pRR)
  #     gam_k <- gamR_init[,k]
  #     V <- diag(SD[i,], nrow=J)
  #     X <- matrix(0, J, xcols * J)
  #     W <- matrix(0, J, nw * J)
  #     for (j in 1:J) {
  #       X[j, ((j-1)*xcols+1):(j*xcols)] <- x_i
  #       W[j, ((j-1)*nw+1):(j*nw)] <- w_i
  #     }
  #     Xstar <- cbind(X,W)
  #     ypred_i <- drop(Xstar %*% theta_init)
  #     resid_i <- Outcome[i,] - ypred_i - drop(W %*% gam_k)
  #     dAd <- ntk * tcrossprod(resid_i) + (ntk - 1) * V %*% RR %*% V
  #     siginvm <- diag(1 / SD[i,], nrow=J)
  #     qq <- qq + (siginvm %*% dAd %*% siginvm)
  #   }
  #   fx_vrho <- function(vRho) {
  #     z <- tanh(vRho)
  #     pRRho <- vecrinv(z, as.integer(J))
  #     diag(pRRho) <- 1
  #     Rhop <- pRho_to_Rho(pRRho)
  #     Rhop <- as.matrix(Matrix::nearPD(Rhop)$mat)
  #     Rhopinv <- chol2inv(chol(Rhop))
  #     loglik <- -0.5 * sum(qq * Rhopinv) - 0.5 * sumNpt * determinant(Rhop)$modulus[1]
  #     for (i in 1:J) {
  #       ii <- i - 1
  #       iR <- J - 2 - floor(sqrt(-8.0*ii + 4.0*(J*(J-1))-7.0)/2.0 - 0.5);
  #       iC <- ii + iR + 1 - (J*(J-1))/2 + ((J-iR)*((J-iR)-1))/2;
  #       loglik <- loglik + 0.5 * (J + 1 - abs(iC-iR)) * log1p(-z[i]^2)
  #     }
  #     loglik
  #   }
  #   ff <- optim(pars, fx_vrho, control=list(fnscale=-1), hessian=TRUE)
  #   fisher_info <- solve(-ff$hessian)
  #   fisher_info <- Matrix::nearPD(fisher_info)$mat
  #   fisher_chol <- t(chol(fisher_info))
  # }

  if (any(eigen(Omega_init, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
    stop("The initial value for Omega is not positive definite")
  }
  if (any(eigen(Rho_init, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
    stop("The initial value for Omega is not positive definite")
  }
  
  if (scale_x) {
    XCovariate_ <- scale(XCovariate, center = FALSE, scale = TRUE)
  } else {
    XCovariate_ <- XCovariate
  }

  mcmctime <- system.time({
    if (second.exist) {
      if (fmodel == 1) {
        fout <- .Call(
          `_metapack_fmodel1p`,
          as.matrix(Outcome),
          as.matrix(SD),
          as.matrix(XCovariate_),
          as.matrix(WCovariate),
          as.integer(Treat.n),
          as.integer(Trial.n),
          as.integer(group.n),
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
          `_metapack_fmodel2p`,
          as.matrix(Outcome),
          as.matrix(SD),
          as.matrix(XCovariate_),
          as.matrix(WCovariate),
          as.integer(Treat.n),
          as.integer(Trial.n),
          as.integer(group.n),
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
          as.double(R_stepsize),
          as.double(theta_init),
          as.matrix(gamR_init),
          as.matrix(Omega_init),
          as.logical(verbose)
        )
    } else if (fmodel == 3) {
      fout <- .Call(
          `_metapack_fmodel2p5p`,
          as.matrix(Outcome),
          as.matrix(SD),
          as.matrix(XCovariate_),
          as.matrix(WCovariate),
          as.integer(Treat.n),
          as.integer(Trial.n),
          as.integer(group.n),
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
          as.double(R_stepsize),
          as.double(theta_init),
          as.matrix(gamR_init),
          as.matrix(Omega_init),
          as.logical(verbose)
        )
      } else if (fmodel == 4) {
        fout <- .Call(
          `_metapack_fmodel3pp`,
          as.matrix(Outcome),
          as.matrix(SD),
          as.matrix(XCovariate_),
          as.matrix(WCovariate),
          as.integer(Treat.n),
          as.integer(Trial.n),
          as.integer(group.n),
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
          `_metapack_fmodel4p`,
          as.matrix(Outcome),
          as.matrix(SD),
          as.matrix(XCovariate_),
          as.matrix(WCovariate),
          as.integer(Treat.n),
          as.integer(Trial.n),
          as.integer(group.n),
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
    } else {
      if (fmodel == 1) {
        fout <- .Call(
          `_metapack_fmodel1`,
          as.matrix(Outcome),
          as.matrix(SD),
          as.matrix(XCovariate_),
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
          as.matrix(XCovariate_),
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
          as.double(R_stepsize),
          as.double(theta_init),
          as.matrix(gamR_init),
          as.matrix(Omega_init),
          as.logical(verbose)
        )
  	} else if (fmodel == 3) {
  		fout <- .Call(
          `_metapack_fmodel2p5`,
          as.matrix(Outcome),
          as.matrix(SD),
          as.matrix(XCovariate_),
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
          as.double(R_stepsize),
          as.double(theta_init),
          as.matrix(gamR_init),
          as.matrix(Omega_init),
          as.logical(verbose)
        )
      } else if (fmodel == 4) {
        fout <- .Call(
          `_metapack_fmodel3`,
          as.matrix(Outcome),
          as.matrix(SD),
          as.matrix(XCovariate_),
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
          as.matrix(XCovariate_),
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
    }
  })
  xcc <- if (!is.null(colnames(XCovariate))) colnames(XCovariate) else paste0("beta", 1:ncol(XCovariate))
  wcc <- if (!is.null(colnames(WCovariate))) colnames(WCovariate) else paste0("gam", 1:ncol(WCovariate))
  if (!is.null(colnames(XCovariate)) && !is.null(colnames(WCovariate))) {
    if (!is.null(group)) {
      rownames(fout$theta) <- c(rep(colnames(XCovariate), J), paste0(rep(colnames(WCovariate), J),"*(1-2nd)"), paste0(rep(colnames(WCovariate), J),"*2nd"))
    } else {
      rownames(fout$theta) <- c(paste0(rep(xcc, J), "_", rep(1:J, each=length(xcc))),
             paste0(rep(wcc, 2*J), rep(rep(c("*(1-2nd)", "*2nd"), each = length(wcc)), J), "_", rep(1:J, each = 2*length(wcc))))
    }
  }
  out <- list(
    Outcome = Outcome,
    SD = SD,
    Npt = Npt,
    XCovariate = XCovariate,
    WCovariate = WCovariate,
    Treat = Treat.n,
    Trial = Trial.n,
    group = group,
    TrtLabels = Treat.order,
    TrialLabels = Trial.order,
    GroupLabels = group.order,
    K = K,
    T = T,
    fmodel = fmodel,
    scale_x = scale_x,
    prior = priorvals,
    mcmctime = mcmctime,
    mcmc = mcvals,
    mcmc.draws = fout
  )
  class(out) <- "bayes.parobs"
  out
}