sde.init <-
function(data, dt, k, m, par.index, params, debug = FALSE) {
  nobs <- nrow(data)
  ndims <- ncol(data)
  if(missing(m)) m <- 2^k-1
  ncomp <- (nobs-1)*(m+1)+1
  init.data <- matrix(NA, ncomp, ndims)
  colnames(init.data) <- colnames(data)
  # interpolation to create missing data
  if(debug) browser()
  dt1 <- length(dt) == 1
  if(dt1) dt <- rep(dt, nobs-1)
  if(length(dt) != nobs-1) stop("Incorrect specification of dt.")
  dtnew <- rep(dt/(m+1), each = m+1)
  told <- cumsum(c(0, dt))
  tnew <- cumsum(c(0, dtnew))
  for(ii in 1:ndims) {
    init.data[,ii] <- approx(x = told, y = data[,ii],
                             xout = tnew)$y
  }
  if(dt1) dtnew <- dtnew[1]
  if(missing(par.index)) par.index <- ndims
  par.ind <- rep(0, ncomp)
  par.ind[seq(1, ncomp, len = nobs)] <- par.index
  ans <- list(data = init.data, dt = dtnew, par.index = par.ind)
  if(!missing(params)) ans$params <- params
  ans
}
