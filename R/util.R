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
  return(x)
}

vhpd <- function(x, level = 0.95) {
    n <- length(x)
    cl <- level
    gap <- max(1, min(n - 1, round(n * level))) / n
    if (level > gap) {
      warning("The desired level cannot be reached with the provided posterior sample size.\n The level should be smaller than ", gap, ". ", "Forcing the HPD level to ", gap-1.0e-4, "\n")
      cl <- gap-1.0e-4
    }
    alpha <- 1 - cl
    out <- .Call(`_metapack_vhpd`, as.vector(x), as.double(alpha))
    attr(out, "Empirical level") <- gap
    return(out)
}

mhpd <- function(x, level = 0.95) {
    n <- ncol(x)
    gap <- max(1, min(n - 1, round(n * level))) / n
    cl <- level
    if (level > gap) {
      warning("The desired level cannot be reached with the provided posterior sample size.\n The level should be smaller than ", gap, ". ", "Forcing the HPD level to ", gap-1.0e-4, "\n")
      cl <- gap-1.0e-4
    }
    alpha <- 1 - cl
    out <- .Call(`_metapack_mhpd`, as.matrix(x), as.double(alpha))
    attr(out, "Empirical level") <- gap
    return(out)
}

hpdarray <- function(A, level = 0.95) {
    nR <- nrow(A)
    nC <- ncol(A)
    n <- dim(A)[3]
    gap <- max(1, min(n - 1, round(n * level))) / n
    cl <- level
    if (level > gap) {
      warning("The desired level cannot be reached with the provided posterior sample size.\n The level should be smaller than ", gap, ". ", "Forcing the HPD level to ", gap-1.0e-4, "\n")
      cl <- gap-1.0e-4
    }
    alpha <- 1 - cl
  out <- array(0, dim = c(nR, nC, 2))
  dimnames(out)[[3]] <- c("lower", "upper")
  for (iC in 1:nC) {
        hpd_ic <- .Call(`_metapack_mhpd`, as.matrix(A[,iC,]), as.double(alpha))
    out[,iC,1] <- hpd_ic[,1]
    out[,iC,2] <- hpd_ic[,2]
  }
  attr(out, "Empirical level") <- gap 
  return(out)
}

ciarray <- function(A, level = 0.95) {
  nR <- nrow(A)
  nC <- ncol(A)
  out <- array(0, dim = c(nR, nC, 2))
  dimnames(out)[[3]] <- c("lower", "upper")
  sig.level <- 1 - level
  out[,,1] <- apply(A, c(1,2), function(xx) quantile(xx, prob = sig.level / 2))
  out[,,2] <- apply(A, c(1,2), function(xx) quantile(xx, prob = 1 - sig.level / 2))
  return(out)
}
