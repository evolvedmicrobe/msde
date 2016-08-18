dcop <-
function(x, par, dens.x, dens.y, Rho, mean, sd, log = FALSE,
                 decomp = FALSE, max.Z = 10, debug = FALSE) {
  if(missing(par)) par <- cop.par(dens.x = dens.x, dens.y = dens.y,
                                  Rho = Rho, mean = mean, sd = sd)
  dens.x <- par$dens.x
  dens.y <- par$dens.y
  ldens.y <- par$ldens.y
  Dens.y <- par$Dens.y
  Rho <- par$Rho
  mean <- par$mean
  sd <- par$sd
  dx <- par$dx
  rx <- par$rx
  nrv <- length(dens.x)
  ndens <- sapply(dens.x, length)
  # x
  if(!is.matrix(x)) {
    x <- as.matrix(x)
    if(nrv != 1) x <- t(x)
  }
  if(ncol(x) != nrv) stop("nrv specified by x and dens do not agree.")
  nx <- nrow(x)
  # pdf and cdf calculations
  f <- matrix(NA, nx, nrv)
  F <- f
  if(debug) browser()
  for(i in 1:nrv) {
    ind <- (x[,i] - rx[1,i])%/%dx[i] + 1
    iind <- pmax(pmin(ind, ndens[i]), 1)
    p <- (x[,i] - rx[1,i])%%dx[i]
    f[,i] <- ifelse(1 <= ind & ind <= ndens[i],
                    ldens.y[[i]][iind],
                    dnorm(x[,i], mean = mean[i], sd = sd[i], log = TRUE))
    F[,i] <- ifelse(1 <= ind & ind <= ndens[i],
                    Dens.y[[i]][iind] + p*dens.y[[i]][iind],
                    pnorm(x[,i], mean = mean[i], sd = sd[i]))
  }
  max.Z <- abs(max.Z)
  Z <- qnorm(F)
  Z <- pmax(pmin(Z, max.Z), -max.Z)
  ldens.z <- dmvnorm(Z, sigma = Rho, log = TRUE)
  jac.z <- rowSums(f) - dmvnorm(Z, sigma = diag(nrv), log = TRUE)
  if(decomp) return(list(Z = Z, ldens.z = ldens.z, jac.z = jac.z))
  ans <- ldens.z + jac.z
  if(!log) ans <- exp(ans)
  ans
}
