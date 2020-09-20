#' Fit test function
#' 
#' This is a function for testing the Metropolis-Hastings algorithm
#' @author Daeyoung Lim, \email{daeyoung.lim@uconn.edu}
#' @param X a matrix of data whose columns are X_i
#' @param ndiscard number of burn-in iterations
#' @param nskip number of thinning iterations
#' @param nkeep number of posterior samples
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
#' 			mcmc=list(ndiscard=100000,nskip=1,nkeep=20000), control(delta_rep=100, rho_rep=200), verbose = TRUE)
#' }
#' @export
test <- function(VSV,
				 vrtk,
				 j,
				 J,
				 iR,
				 iC,
				 ntk) {
	mcmctime <- system.time({
				fout <- .Call(`_metapack_testfun`,
					  as.matrix(VSV),
				   	  as.double(vrtk),
				      as.integer(j),
				      as.integer(J),
				      as.integer(iR),
				      as.integer(iC),
				      as.double(ntk))
			})
	out <- list(mcmc.draws = fout,
				mcmctime = mcmctime)
	out
}

# test <- function(a) {
# 	mcmctime <- system.time({
# 				fout <- .Call(`_metapack_testfun`,
# 					  as.double(a))
# 			})
# 	out <- list(mcmc.draws = fout,
# 				mcmctime = mcmctime)
# 	out
# }
# test <- function(X, ndiscard, nskip, nkeep) {
# 	if (missing(ndiscard)) ndiscard <- 5000L
# 	if (missing(nskip)) nskip <- 1L
# 	if (missing(nkeep)) nkeep <- 20000L

# 	mcmctime <- system.time({
# 				fout <- .Call(`_metapack_testfun`,
# 					  as.matrix(X),
# 					  as.integer(ndiscard),
# 					  as.integer(nskip),
# 					  as.integer(nkeep))
# 			})
# 	out <- list(mcmc.draws = fout,
# 				mcmctime = mcmctime,
# 				X = X,
# 				ndiscard = ndiscard,
# 				nskip = nskip,
# 				nkeep = nkeep)
# 	out
# }