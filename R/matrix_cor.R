#' `meta_matrix_cor` fits the meta-regression model where the correlation matrices are covariate-modeled
#' 
#' @description
#' `meta_matrix_cor` extends the meta-regression model so that the correlation matrix can be modeled using covariates. There are three options for modeling the correlation matrix: (1) matrix exponential, (2) partial correlation, and (3) Cholesky factor further transformed into hyperspherical coordinates. Note that the latter two are not *permutation invariant*. This means that \eqn{(y_1, y_2, y_3)} and \eqn{(y_2, y_1, y_3)} (or any permutation of the random variables) will have a nonsignificant effect on the estimation and the statistical inference.
#' 
#' @author Daeyoung Lim, \email{daeyoung.lim@uconn.edu}
#' @param Outcome the aggregate mean of the responses for each arm of every study.
#' @param SD the standard deviation of the responses for each arm of every study.
#' @param XCovariate the aggregate covariates for the fixed effects.
#' @param WCovariate the aggregate covariates for the random effects.
#' @param ZCovariate the aggregate covariates for modeling the correlation matrices.
#' @param Treat the treatment identifiers. This is equivalent to the arm number of each study. The number of unique treatments must be equal across trials. The elements within will be coerced to consecutive integers.
#' @param Trial the trial identifiers. This is equivalent to the arm labels in each study. The elements within will be coerced to consecutive integers
#' @param Npt the number of observations/participants for a unique `(t,k)`, or each arm of every trial.
#' @param transform_type the model number. The possible values for `transform_type` are 1 to 3, each indicating a different prior specification for \eqn{\rho_{tk}}. It will default to M1, `transform_type=1` if not specified at function call. See the following model descriptions.
#' 
#' + `transform_type=1` - \eqn{\rho_{tk}} is modeled as a matrix exponential of a real symmetric matrix whose upper-diagonal elements are unconstrained and modeled as a determistic function of the covariates.
#' + `transform_type=2` - \eqn{\rho_{tk}} is modeled as a correlation matrix whose partial correlation matrix is Fisher's z-transformed so the upper-diagonal elements are unconstrained. The free variables are then modeled as a determistic function of the covariates.
#' + `transform_type=3` - \eqn{\rho_{tk}} is modeled as a correlation matrix whose Cholesky factor is modeled as hyperspherical coordinates. The hyperspherical coordinates are a deterministic function of the covariates.

#' @param prior (Optional) a list of hyperparameters. Despite `theta` in every model, each `fmodel`, along with the `group` argument, requires a different set of hyperparameters.
#' @param mcmc (Optional) a list for MCMC specification. `ndiscard` is the number of burn-in iterations. `nskip` configures the thinning of the MCMC. For instance, if `nskip=5`, `bayes.parobs` will save the posterior sample every 5 iterations. `nkeep` is the size of the posterior sample. The total number of iterations will be `ndiscard + nskip * nkeep`.
#' @param control (Optional) a list of tuning parameters for [the Metropolis-Hastings algorithm](https://en.wikipedia.org/wiki/Metropolis-Hastings_algorithm). `R`, and `delta` are sampled through either localized Metropolis algorithm or delayed rejection robust adaptive Metropolis algorithm. `*_stepsize` with the asterisk replaced with one of the names above specifies the stepsize for determining the sample evaluation points in the localized Metropolis algorithm.
#' @param init (Optional) a list of initial values for the parameters to be sampled: `theta`, `gamR`, and `Omega`.
#' @param Treat_order (Optional) a vector of unique treatments to be used for renumbering the `Treat` vector. The first element will be assigned treatment zero, potentially indicating placebo. If not provided, the numbering will default to an alphabetical/numerical order.
#' @param Trial_order (Optional) a vector of unique trials. The first element will be assigned zero. If not provided, the numbering will default to an alphabetical/numerical order.
#' @param group (Optional) a vector containing binary variables for \eqn{u_{tk}}. If not provided, `bayes.parobs` will assume that there is no grouping and set \eqn{u_{tk}=0} for all `(t,k)`.
#' @param group_order (Optional) a vector of unique group labels. The first element will be assigned zero. If not provided, the numbering will default to an alphabetical/numerical order. `group_order` will take effect only if `group` is provided by the user.
#' @param scale_x (Optional) a logical variable indicating whether `XCovariate` should be scaled/standardized. The effect of setting this to `TRUE` is not limited to merely standardizing `XCovariate`. The following generic functions will scale the posterior sample of `theta` back to its original unit: `plot`, `fitted`, `summary`, and `print`.
#' @param verbose (Optional) a logical variable indicating whether to print the progress bar during the MCMC sampling.

#' @examples
#' set.seed(2797542)
#' data("cholesterol")
#' f_1 <- 'pldlc + phdlc + ptg | sdldl + sdhdl + sdtg ~ 0 + bldlc + bhdlc + btg +
#'   age + durat + white + male + dm + ns(n) | treat | treat + trial + onstat'
#' out_1 <- bmeta_analyze(as.formula(f_1), data = cholesterol,
#'   prior = list(model="NoRecovery"),
#'   mcmc = list(ndiscard = 3, nskip = 1, nkeep = 1),
#'   control=list(scale_x = TRUE, verbose=FALSE))
#'
#' set.seed(2797542)
#' data("TNM")
#' TNM$group <- factor(match(TNM$treat, c("PBO", "R"), nomatch = 0))
#' f_2 <- 'ptg | sdtg ~
#'   0 + bldlc + bhdlc + btg + age + white + male + bmi +
#'   potencymed + potencyhigh + durat + ns(n) |
#'   scale(bldlc) + scale(btg) + group | treat  + trial'
#' out_2 <- bmeta_analyze(as.formula(f_2), data = TNM,
#'   mcmc = list(ndiscard = 1, nskip = 1, nkeep = 1),
#'   control=list(scale_x = TRUE, verbose=FALSE))
#' @references 
#' Yao, H., Kim, S., Chen, M. H., Ibrahim, J. G., Shah, A. K., & Lin, J. (2015). Bayesian inference for multivariate meta-regression with a partially observed within-study sample covariance matrix. *Journal of the American Statistical Association*, **110(510)**, 528-544.
#' 
#' Li, H., Chen, M. H., Ibrahim, J. G., Kim, S., Shah, A. K., Lin, J., & Tershakovec, A. M. (2019). Bayesian inference for network meta-regression using multivariate random effects with applications to cholesterol lowering drugs. *Biostatistics*, **20(3)**, 499-516.
#' 
#' Li, H., Lim, D., Chen, M. H., Ibrahim, J. G., Kim, S., Shah, A. K., & Lin, J. (2021). Bayesian network meta-regression hierarchical models using heavy-tailed multivariate random effects with covariate-dependent variances. *Statistics in Medicine*.
#' @md
#' @export
'meta_matrix_cor' <- function(Outcome, SD, XCovariate, WCovariate, ZCovariate, Treat, Trial, Npt, transform_type = 1, prior = list(), mcmc = list(), control = list(), init = list(), Treat_order = NULL, Trial_order = NULL, group = NULL, group_order = NULL, scale_x = FALSE, verbose = FALSE) {
    # initial version: DY 2021/09/21
    # - DY Sep 21 2021: first version drafted
    if (!is(Outcome, "matrix")) {
        tmp <- try(Outcome <- model.matrix(~0 + ., data = Outcome), silent = TRUE)
        if (is(tmp, "try-error")) {
            stop("Outcome must be a matrix or able to be coerced to a matrix")
        }
    }
    if (!is(SD, "matrix")) {
        tmp <- try(SD <- model.matrix(~0 + ., data = SD), silent = TRUE)
        if (is(tmp, "try-error")) {
            stop("SD must be a matrix or able to be coerced to a matrix")
        }
    }
    if (!is(XCovariate, "matrix")) {
        tmp <- try(XCovariate <- model.matrix(~0 + ., data = XCovariate), silent = TRUE)
        if (is(tmp, "try-error")) {
            stop("XCovariate must be a matrix or able to be coerced to a matrix")
        }
    }
    if (!is(WCovariate, "matrix")) {
        tmp <- try(WCovariate <- model.matrix(~0 + ., data = WCovariate), silent = TRUE)
        if (is(tmp, "try-error")) {
            stop("WCovariate must be a matrix or able to be coerced to a matrix")
        }
    }
    if (!is(ZCovariate, "matrix")) {
        tmp <- try(ZCovariate <- model.matrix(~0 + ., data = ZCovariate), silent = TRUE)
        if (is(tmp, "try-error")) {
            stop("ZCovariate must be a matrix or able to be coerced to a matrix")
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
            stop(paste(sQuote("group"), "must be binary"))
        }
        nr <- 2
        if (is.null(group_order)) {
            group.order <- sort(sl)
        } else {
            group.order <- group_order
        }
        group.n <- relabel.vec(group, group.order) - 1
    } else {
        group.n <- NULL
        group.order <- NULL
    }


    mcvals <- list(ndiscard = 5000L, nskip = 1L, nkeep = 20000L)
    mcvals[names(mcmc)] <- mcmc
    ndiscard <- mcvals$ndiscard
    nskip <- mcvals$nskip
    if (nskip < 1) {
    stop(paste0(sQuote("nskip"), "can't be smaller than 1"))
    }
    nkeep <- mcvals$nkeep

    ctrl <- list(R_stepsize = 0.02, delta_stepsize = 0.2, TOL = 1.0e-6)
    ctrl[names(control)] <- control
    R_stepsize <- ctrl$R_stepsize
    delta_stepsize <- ctrl$delta_stepsize
    TOL <- ctrl$TOL

    J <- ncol(Outcome)
    nw <- ncol(WCovariate)

  priorvals <- list(c0 = 1.0e05, dj0 = 0.1 + nw, a0 = 0.1, b0 = 0.1, a1 = 0.1, b1 = 0.1,
                      a2 = 0.1, b2 = 0.1, a3 = 0.1, b3 = 0.1,
                      Omega0 = diag(10, nw))
    priorvals[names(prior)] <- prior
    a0 <- priorvals$a0
    b0 <- priorvals$b0
    a1 <- priorvals$a1
    b1 <- priorvals$b1
    a2 <- priorvals$a2
    b2 <- priorvals$b2
    a3 <- priorvals$a3
    b3 <- priorvals$b3
    c0 <- priorvals$c0
    dj0 <- priorvals$dj0
    Omega0 <- priorvals$Omega0

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
        theta = numeric((xcols + nw * nr) * J),
        gamR = matrix(0, nw * nr * J, K), Omega = diag(1, nrow = nw * nr * J)
    )
    init_final[names(init)] <- init
    theta_init <- init_final$theta
    gamR_init <- init_final$gamR
    Omega_init <- init_final$Omega
    if (length(theta_init) != (xcols + nw * nr) * J) {
        stop(paste("theta initialized with length", sQuote(length(theta_init)), "but", sQuote((xcols + nw * nr) * J), "wanted"))
    }
    if (dim(Omega_init)[1] != nw * nr * J || dim(Omega_init)[2] != nw * nr * J) {
        stop(paste("Omega initialized with dimensions", sQuote(dim(Omega_init)), "but", sQuote(nw * nr * J), "wanted"))
    }
    if (dim(gamR_init)[1] != nw * nr * J) {
        stop(paste("gamR initialized with", sQuote(dim(gamR_init)[1]), "rows but", sQuote(nw * nr * J), "wanted"))
    }
    if (dim(gamR_init)[2] != K) {
        stop(paste("gamR initialized with", sQuote(dim(gamR_init)[2]), "columns but", sQuote(K), "wanted"))
    }

    if (any(eigen(Omega_init, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
        stop(paste("The initial value for", sQuote("Omega"), "is not positive definite"))
    }

    if (scale_x) {
        XCovariate_ <- scale(XCovariate, center = FALSE, scale = TRUE)
    } else {
        XCovariate_ <- XCovariate
    }

    mcmctime <- system.time({
        if (second.exist) {
            fout <- .Call(`_metapack_fmodel_corr_modeling`,
                        as.matrix(Outcome),
                        as.matrix(SD),
                        as.matrix(XCovariate_),
                        as.matrix(WCovariate),
                        as.matrix(ZCovariate),
                        as.integer(Treat.n),
                        as.integer(Trial.n),
                        as.integer(group.n),
                        as.double(Npt),
                        as.double(c0),
                        as.double(dj0),
                        as.double(a0),
                        as.double(b0),
                        as.double(a1),
                        as.double(b1),
                        as.double(a2),
                        as.double(b2),
                        as.double(a3),
                        as.double(b3),
                        as.matrix(Omega0),
                        as.integer(K),
                        as.integer(T),
                        as.integer(ndiscard),
                        as.integer(nskip),
                        as.integer(nkeep),
                        as.integer(transform_type),
                        as.double(delta_stepsize),
                        as.double(R_stepsize),
                        as.double(TOL),
                        as.double(theta_init),
                        as.matrix(gamR_init),
                        as.matrix(Omega_init),
                        as.logical(verbose)
                    )
        }
    })
    xcc <- if (!is.null(colnames(XCovariate))) colnames(XCovariate) else paste0("beta", 1:ncol(XCovariate))
    wcc <- if (!is.null(colnames(WCovariate))) colnames(WCovariate) else paste0("gam", 1:ncol(WCovariate))
    if (is.null(group)) {
        rownames(fout$theta) <- c(paste0(rep(xcc, J), "_", rep(1:J, each = length(xcc))), paste0(rep(wcc, J), "_", rep(1:J, each = length(wcc))))
    } else {
        rownames(fout$theta) <- c(paste0(rep(xcc, J), "_", rep(1:J, each = length(xcc))),
        paste0(rep(wcc, 2 * J), rep(rep(c("*(1-2nd)", "*2nd"), each = length(wcc)), J), "_", rep(1:J, each = 2 * length(wcc))))
    }

    out <- list(
        Outcome = Outcome,
        SD = SD,
        Npt = Npt,
        XCovariate = XCovariate_,
        WCovariate = WCovariate,
        ZCovariate = ZCovariate,
        Treat = Treat.n,
        Trial = Trial.n,
        group = group.n,
        TrtLabels = Treat.order,
        TrialLabels = Trial.order,
        GroupLabels = group.order,
        K = K,
        T = T,
        fmodel = fmodel,
        scale_x = scale_x,
        prior = priorvals,
        control = ctrl,
        mcmctime = mcmctime,
        mcmc = mcvals,
        mcmc.draws = fout
    )
    class(out) <- "matrixcorr"
    out
}