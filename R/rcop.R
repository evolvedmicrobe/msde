rcop <-
function(n, par, dens.x, dens.y, Rho, mean, sd, debug = FALSE) {
  if(missing(par)) par <- cop.par(dens = dens.x, dens.y = dens.y,
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
  if(debug) browser()
  # Simulate Uniforms
  P <- pnorm(rmvnorm(n, sigma = Rho))
  X <- matrix(NA, n, nrv)
  nm <- which.min(sapply(par, function(l) is.null(names(l))))
  colnames(X) <- names(par[[nm]])
  for(i in 1:nrv) {
    tmp <- c(0, Dens.y[[i]], 1)
    tmpi <- !duplicated(tmp)
    ind <- as.numeric(cut(P[,i], breaks = tmp[tmpi])) - 1
    iind <- pmax(pmin(ind, sum(tmpi)-3), 1)
    Dy <- diff(Dens.y[[i]][tmpi])[iind]
    PDy <- P[,i] - Dens.y[[i]][tmpi][iind]
    X[,i] <- ifelse(1 <= ind & ind <= sum(tmpi)-3,
                    dens.x[[i]][tmpi][iind] + (PDy/Dy - 1/2)*dx[i],
                    qnorm(P[,i], mean = mean[i], sd = sd[i]))
  }
  X
}
