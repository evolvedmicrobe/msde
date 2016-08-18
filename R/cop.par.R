cop.par <-
function(X, dens.x, dens.y, Rho, mean, sd, from, to, n, ...,
                    zero.dens = FALSE, sd.infl = 1/5, debug = FALSE) {
  if(!missing(X)) {
    X <- as.matrix(X)
    d <- ncol(X)
    if(missing(mean)) mean <- colMeans(X)
    calc.mean <- FALSE
    if(missing(sd)) {
      sd <- 1
      sd <- apply(X, 2, sd)*sd.infl
    }
    calc.sd <- FALSE
    if(missing(Rho) & (missing(dens.x) | missing(dens.y))) {
      tmp <- apply(X, 2,
                   function(x) qqnorm(x, plot.it = FALSE)$x)
      Rho <- cor(tmp)
    }
    if(missing(dens.x) | missing(dens.y)) {      
      if(!missing(from)) {
        from <- rep(from, len = d)
        to <- rep(to, len = d)
        if(missing(n)) n <- 512
        n <- rep(n, len = d)
        dens.x <- vector("list", d)
        names(dens.x) <- colnames(X)
        dens.y <- dens.x
        for(i in 1:d) {
          dens <- density(X[,i], n = n[i], from = from[i], to = to[i], ...)
          dens.x[[i]] <- dens$x
          dens.y[[i]] <- dens$y
        }
      } else {
        if(missing(n)) n <- 512
        n <- rep(n, len = d)
        dens.x <- vector("list", d)
        names(dens.x) <- colnames(X)
        dens.y <- dens.x
        for(i in 1:d) {
          dens <- density(X[,i], n = n[i], ...)
          dens.x[[i]] <- dens$x
          dens.y[[i]] <- dens$y
        }
      }
    }
  }
  if(debug) browser()
  # format arguments
  # dens.y
  if(is.matrix(dens.y)) {
    dens.y <- unlist(apply(dens.y, 2, function(x) list(x)), recursive = FALSE)
  } else if(is.numeric(dens.y)) dens.y <- list(dens.y)
  if(!is.list(dens.y)) stop("Incorrect specification for dens.y.")
  if(!zero.dens) dens.y <- sapply(dens.y, function(x) pmax(x, min(x[x > 0])),
                                  simplify = FALSE)
  ldens.y <- sapply(dens.y, log, simplify = FALSE)
  nrv <- length(dens.y)
  ndens <- sapply(dens.y, length)
  # dens.x
  if(is.matrix(dens.x)) {
    dens.x <- unlist(apply(dens.x, 2, function(x) list(x)), recursive = FALSE)
  } else if(is.numeric(dens.x)) dens.x <- list(dens.x)
  if(!is.list(dens.x)) stop("Incorrect specification for dens.x.")
  # space between dens evaluations
  dx <- sapply(dens.x, function(x) x[2]-x[1])
  # range of density: outside is normal approx
  rx <- sapply(dens.x, range)
  rx[1,] <- rx[1,] - dx/2
  rx[2,] <- rx[2,] + dx/2
  if(length(dens.x) != nrv)
    stop("nrv specified by dens.x and dens.y do not agree.")
  if(any(sapply(dens.x, length) != ndens))
    stop("ndens specified by dens.x and dens.y do not agree.")
  # mean and sd
  if(missing(mean)) {
    calc.mean <- TRUE
    mean <- rep(NA, nrv)
  } else {
    if(length(mean) != nrv) stop("mean has wrong length.")
    calc.mean <- FALSE
  }
  if(missing(sd)) {
    calc.sd <- TRUE
    sd <- rep(NA, nrv)
  } else {
    if(length(sd) != nrv) stop("sd has wrong length.")
    calc.sd <- FALSE
  }
  if(missing(Rho)) Z <- matrix(NA, nrow(X), ncol(X))
  # renormalize dens.y to account for tails; calculate CDF, mean, and sd
  Dens.y <- vector("list", nrv)
  names(Dens.y) <- names(dens.y)
  for(i in 1:nrv) {
    # calculate mean and sd
    if(calc.mean) mean[i] <- sum(dens.x[[i]]*dens.y[[i]])*dx[i]
    if(calc.sd)
      sd[i] <- sqrt(sum(dens.x[[i]]^2*dens.y[[i]])*dx[i] - mean[i]^2)*sd.infl
    # tail probability not contained in dens estimate
    tail.probs <- pnorm(abs(rx[,i] - mean[i])/sd[i], lower.tail = FALSE)
    dens.y[[i]] <- dens.y[[i]]/sum(dens.y[[i]]*dx[i])*(1-sum(tail.probs))
    Dens.y[[i]] <- cumsum(c(tail.probs[1], dens.y[[i]]*dx[i]))
    # calculate Rho
    if(missing(Rho)) {
      ind <- (X[,i] - rx[1,i]) %/% dx[i] + 1
      p <- (X[,i] - rx[1,i]) %% dx[i]
      iind <- pmax(pmin(ind, ndens[i]), 1)
      Z[,i] <- ifelse(1 <= ind & ind <= ndens[i],
                      qnorm(Dens.y[[i]][iind] + p*dens.y[[i]][iind]),
                      (X[,i] - mean[i])/sd[i])
    }
  }
  # Rho
  if(missing(Rho)) Rho <- cor(Z)
  Rho <- cov2cor(Rho)
  list(dens.x = dens.x, dx = dx, rx = rx,
       dens.y = dens.y, ldens.y = ldens.y, Dens.y = Dens.y,
       Rho = Rho, mean = mean, sd = sd)
}
