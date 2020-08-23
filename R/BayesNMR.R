#' Fit Bayesian Network Meta-Regression Hierarchical Models Using Heavy-Tailed Multivariate Random Effects with Covariate-Dependent Variances
#' 
#' This is a function for running the Markov chain Monte Carlo algorithm for the BNMHHtMRe Model. The first six arguments are required.
#' @author Daeyoung Lim, \email{daeyoung.lim@uconn.edu}
#' @param Outcome aggregate mean of the responses for each arm of each study
#' @param SD standard deviation of the responses for each arm of each study
#' @param Covariate aggregate covariates for the mean component
#' @param Study study number in integers
#' @param Treat treatment applied; equivalent to the arm number of each study
#' @param Npt number of observations per trial
#' @param groupInfo list of grouping information; the control(baseline) group must start from 0; the aggregate covariates 'z' explaining the variance of the random effect of the t-th treatment will be construct based on this grouping information
#' @param prior list of hyperparameters; when not given, algorithm will run in default setting
#' @param mcmc list of MCMC-related parameters: number of burn-ins (ndiscard), number of thinning(nskip), and posterior sample size (nkeep)
#' @param add.z additional covariates other than the grouping vectors that should be column-concatenated to 'z'. This should have the same number of rows as 'y', and 'x'
#' @param scale.x logical variable for scaling x. Default to TRUE. If not scaled, the gamma[1] (different than gam in the function) cannot be interpreted as placebo
#' @param verbose logical variable for printing progress bar. Default to FALSE.
#' @param init initial values for beta (ns + nT dimensional) and phi. Dimensions must be conformant.
#' @return a dataframe with input arguments, posterior samples, Metropolis algorithm acceptance rates, etc
#' @examples
#' \dontrun{
#' data(df)
#' groupInfo <- list(c(0,1), c(2,3), c(4)) # define the variance structure
#' x <- df[,6:10]
#' fit <- bayes.nmr(df$y, df$sd, x, df$ids, df$iarm, df$npt, groupInfo,
#' 			prior = list(c01=1.0e05, c02=4, nu=3),
#' 			mcmc=list(ndiscard=2500,nskip=1,nkeep=10000))
#' }
#' @export
bayes.nmr <- function(Outcome, SD, Covariate, Study, Treat, Npt, groupInfo, prior = list(), mcmc = list(), add.z=list(), scale.x=TRUE, verbose=FALSE, init=list()) {
	mcvals <- list(ndiscard = 5000L, nskip = 1L, nkeep = 20000L)
	mcvals[names(mcmc)] = mcmc
	ndiscard <- mcvals$ndiscard
	nskip <- mcvals$nskip
	nkeep <- mcvals$nkeep
	if (length(groupInfo) == 0) {
		stop("groupInfo is missing.")
	}
	if (min(c(unlist(groupInfo))) != 0) {
		warning("Baseline treatment for groupInfo should start from 0.\nAssuming baseline is not ")
	}

	priorvals <- list(nu = 20, c01 = 1.0e05, c02 = 4)
	priorvals[names(prior)] <- prior
	nu <- priorvals$nu
	c01 <- priorvals$c01
	c02 <- priorvals$c02

	Treat.order <- sort(unique(Treat))
	Treat.n <- relabel.vec(Treat, Treat.order) - 1 # relabel the treatment numbers
	# if (min(unique(Treat)) != 0) {
	# 	stop("The treatments should start from 0. Please adjust accordingly.")
	# }
	nx <- ncol(Covariate)
	nz <- length(groupInfo)
	ns <- length(Outcome)
	K <- length(unique(Study))
	nT <- length(unique(Treat))
	z <- matrix(0, ns, nz)
	if (nz > 0) {
		for (j in 1:nz) {
			for (i in 1:ns) {
				if (Treat[i] %in% groupInfo[[j]]) {
					z[i,j] <- 1
				}
			}
		}
	}
	if (length(add.z) > 0) {
		z <- cbind(z, scale(add.z, center=TRUE, scale=TRUE))
	}
	if (scale.x) {
		Covariate <- scale(Covariate, center = TRUE, scale = TRUE)
	}
	init_final <- list(beta = numeric(nx+nT), phi = numeric(ncol(z)), sig2 = rep(1, ns))
	init_final[names(init)] <- init


	mcmctime <- system.time({
				fout <- .Call(`_metapack_BayesNMR`,
					  as.double(Outcome),
					  as.double(SD),
					  as.matrix(Covariate),
					  as.matrix(z),
					  as.integer(Study),
					  as.integer(Treat.n),
					  as.double(Npt),
					  as.double(nu),
					  as.double(1/c01),
					  as.double(1/c02),
					  as.integer(K),
					  as.integer(nT),
					  as.integer(ndiscard),
					  as.integer(nskip),
					  as.integer(nkeep),
					  as.logical(verbose),
					  as.double(init_final$beta),
					  as.double(init_final$phi),
					  as.double(init_final$sig2))
			})

	out <- list(Outcome = Outcome,
				SD = SD,
				Npt = Npt,
				Covariate = Covariate,
				z = z,
				Study = Study,
				Treat = Treat.n,
				TrtLabels = Treat.order,
				K = K,
				nT = nT,
				groupInfo = groupInfo,
				prior = priorvals,
				mcmctime = mcmctime,
				mcmc = mcvals,
				mcmc.draws = fout)
	class(out) <- "bayesnmr"
	out
}
