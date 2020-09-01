#' Fit Bayesian Inference for Multivariate Meta-Regression With a Partially Observed Within-Study Sample Covariance Matrix
#' 
#' This is a function for running the Markov chain Monte Carlo algorithm for the BMVMR_POCov Model. The first six arguments are required.
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
#' @param control list of parameters for localized Metropolis algorithm: the step sizes for R, Rho, delta, and Delta (R_stepsize, Rho_stepsize, delta_stepsize); they're all 0.2 by default
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
#' XCovariate <- cbind(cholesterol$bldlc, cholesterol$bhdlc, cholesterol$btg, cholesterol$age, cholesterol$durat, cholesterol$white, cholesterol$male, cholesterol$dm)
#' WCovariate <- cbind(1-cholesterol$onstat, cholesterol$trt * (1-cholesterol$onstat), cholesterol$onstat, cholesterol$trt * cholesterol$onstat)
#'
#' fmodel <- 3
#' fit <- bayes.parobs(Outcome, SD, scale(XCovariate, scale=TRUE, center=TRUE), WCovariate, Treat, Trial, Npt, fmodel,
#' 			mcmc=list(ndiscard=100000,nskip=1,nkeep=20000), verbose = TRUE)
#' }
#' @export
bayes.parobs <- function(Outcome, SD, XCovariate, WCovariate, Treat, Trial, Npt, fmodel = 1, prior = list(), mcmc = list(), control = list(), verbose=FALSE) {
	mcvals <- list(ndiscard = 5000L, nskip = 1L, nkeep = 20000L)
	mcvals[names(mcmc)] <- mcmc
	ndiscard <- mcvals$ndiscard
	nskip <- mcvals$nskip
	nkeep <- mcvals$nkeep

	ctrl <- list(R_stepsize = 0.2, Rho_stepsize = 0.2, delta_stepsize = 0.2, L_HMC = 7, eps_HMC = 0.00002)
	ctrl[names(control)] <- control
	R_stepsize <- ctrl$R_stepsize
	Rho_stepsize <- ctrl$Rho_stepsize
	delta_stepsize <- ctrl$delta_stepsize
	L_HMC <- ctrl$L_HMC
	eps_HMC <- ctrl$eps_HMC

	J = ncol(Outcome)
	nw = ncol(WCovariate)
 	priorvals <- list(c0 = 1.0e05, dj0 = 0.1 + nw, d0 = 0.1 + J, s0 = 0.1, a0 = 0.1, b0 = 0.1, Omega0 = diag(10,nw), Sigma0 = diag(10,J), nu0 = 10)
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
	

	Treat.order <- sort(unique(Treat))
	Treat.n <- relabel.vec(Treat, Treat.order) - 1 # relabel the treatment numbers

	Trial.order <- sort(unique(Trial))
	Trial.n <- relabel.vec(Trial, Trial.order) - 1 # relabel the trial numbers

	K <- length(unique(Trial))
	T <- length(unique(Treat))

	mcmctime <- system.time({
				fout <- .Call(`_metapack_BMVMR_POCovHMC`,
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
					  as.double(s0),
					  as.double(nu0),
					  as.double(a0),
					  as.double(b0),
					  as.matrix(Omega0),
					  as.matrix(Sigma0),
					  as.integer(K),
					  as.integer(T),
					  as.integer(fmodel),
					  as.integer(ndiscard),
					  as.integer(nskip),
					  as.integer(nkeep),
					  as.double(R_stepsize),
					  as.double(Rho_stepsize),
					  as.double(delta_stepsize),
					  as.integer(L_HMC),
					  as.double(eps_HMC),
					  as.logical(verbose))
			})
	if (!is.null(colnames(XCovariate)) && !is.null(colnames(WCovariate))) {
		rownames(fout$theta) <- c(rep(colnames(XCovariate), J), rep(colnames(WCovariate), J))
	}
	out <- list(Outcome = Outcome,
				SD = SD,
				Npt = Npt,
				XCovariate = XCovariate,
				WCovariate = WCovariate,
				Treat = Treat,
				Trial = Trial,
				TrtLabels = Treat.order,
				TrialLabels = Trial.order,
				K = K,
				T = T,
				fmodel = fmodel,
				prior = priorvals,
				mcmctime = mcmctime,
				mcmc = mcvals,
				mcmc.draws = fout)
	class(out) <- "bayes.parobs"
	out
}
