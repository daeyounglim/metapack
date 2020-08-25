transform.data <- function(df){
  with(df, {
    #relabel Study names into to a numeric sequence (1:nstudy)
    na <- rle(Study)$lengths
    Study.order <- unique(Study)
    names(Study.order) <- 1:length(Study.order)
    Study <- rep(1:length(unique(Study)), times = na)

    #relabel the treatment according to treatment order specified
    if(is.null(Treat.order)){
      Treat.order <- sort(unique(Treat))
    }
    Treat <- relabel.vec(Treat, Treat.order)
    names(Treat.order) <- 1:length(Treat.order)

    data <- if(response == "normal"){
      cbind(Outcomes, SE, Study, Treat)
    } else {
      cbind(Outcomes, N, Study, Treat)
    }

    ordering <- order(Study, Treat)
    data <- data[ordering,]

    list(data = data, Treat.order = Treat.order, Study.order = Study.order)
  })
}

relabel.vec <- function(x, order)
{
  old.x <- x
  x <- rep(NA, length(old.x))
  for (i in seq(length(order))) x[old.x == order[i]] <- i #relabel studies in numerical order starting with one
  x
}


hpdarray <- function(A, conf.level = 0.95) {
  nR <- nrow(A)
  nC <- ncol(A)
  out <- array(0, dim = c(nR, nC, 2))
  dimnames(out)[[3]] <- c("lower", "upper")
  for (iC in 1:nC) {
    hpd_ic <- coda::HPDinterval(mcmc(t(A[,iC,]), end = dim(A)[3]), prob = conf.level)
    out[,iC,1] <- hpd_ic[,1]
    out[,iC,2] <- hpd_ic[,2]
  }
  out
}

