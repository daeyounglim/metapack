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
#' @param verbose logical variable for printing progress bar. Default to FALSE.
#' @return a dataframe with input arguments, posterior samples, Metropolis algorithm acceptance rates, etc
#' @examples
#' \dontrun{
#' data("cholesterol")
#' Outcome <- cbind(cholesterol$ldlcm, cholesterol$hdlcm, cholesterol$tgm)
#' SD <- cbind(cholesterol$ldlcsd, cholesterol$hdlcsd, cholesterol$tgsd)
#' Trial <- cholesterol$Trial
#' Treat <- cholesterol$Trt
#' Npt <- cholesterol$npt
#' XCovariate <- cbind(cholesterol$bl_ldlc, cholesterol$bl_hdlc, cholesterol$bl_tg, cholesterol$age, cholesterol$Dur, cholesterol$white, cholesterol$male, cholesterol$DM)
#' WCovariate <- cbind(1, cholesterol$Trt)
#'
#' fit <- bayes.parobs(Outcome, SD, XCovariate, WCovariate, Treat, Trial, Npt, 1,
#' 			mcmc=list(ndiscard=2500,nskip=1,nkeep=10000), verbose = TRUE)
#' }
#' @export
bayes.parobs <- function(Outcome, SD, XCovariate, WCovariate, Treat, Trial, Npt, fmodel = 1, prior = list(), mcmc = list(), verbose=FALSE) {
	mcvals <- list(ndiscard = 5000L, nskip = 1L, nkeep = 20000L)
	mcvals[names(mcmc)] = mcmc
	ndiscard <- mcvals$ndiscard
	nskip <- mcvals$nskip
	nkeep <- mcvals$nkeep


	J = ncol(Outcome)
	nw = ncol(WCovariate)
 	priorvals <- list(c0 = 1.0e05, dj0 = 0.1 + nw, d0 = 0.1 + J, s0 = 0.1, Omega0 = diag(10,nw), Sigma0 = diag(10,J), nu0 = 10)
	priorvals[names(prior)] <- prior
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
				fout <- .Call(`_metapack_BMVMR_POCov`,
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
					  as.matrix(Omega0),
					  as.matrix(Sigma0),
					  as.integer(K),
					  as.integer(T),
					  as.integer(fmodel),
					  as.integer(ndiscard),
					  as.integer(nskip),
					  as.integer(nkeep),
					  as.logical(verbose))
			})

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
