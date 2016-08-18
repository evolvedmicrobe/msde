#'@name sde.loglik
#'@title Loglikelihood Function
#'@description Applies the log-likelihood function for \code{model} to \code{x} using \code{theta} and \code{dt}
#'@param model An object of the sde.model class generated using \code{sde.make.model}
#'@param x A matrix of size (nreps x ncomp x ndims) containing data
#'@param dt A numeric vector of size 1 or ncomp-1 containing Delta-t
#'@param theta A matrix of size (nreps x nparams) containing parameters
#'@param debug A boolean (\code{FALSE} by default) if set to \code{TRUE}, will cause the function to open a browser mid-call
#'@return a matrix containing the result of the log-likelihood function
#'@examples
#'# Create the model
#'hest.model <- sde.make.model(list = hestList, model.name = "hest")
#'
#'theta <- c(alpha = .1, gamma = 5, beta = .8, sigma = .6, rho = -.7)
#'Y0 <- c(X = log(100), Z = .1)
#'
#'# simulate data
#'N <- 10
#'burn <- 10
#'dT <- 1/252
#'
#'hsim <- sde.sim(model = hest.model, init.data = Y0, params = theta, dt = dT, dt.sim = dT/100,
#'                N = N, burn = burn, nreps = 1)
#'
#'loglik <- sde.loglik(model = hest.model, x = hsim$data, dt = dT, theta = theta)
#'@export
sde.loglik <- function(model, x, dt, theta, debug = FALSE) {
  if(class(model) != "sde.model")
    stop("Expecting object of class sde.model.  Use sde.make.model to create.")
  # model constants
  ndims <- model$ndims
  data.names <- model$data.names
  nparams <- model$nparams
  param.names <- model$param.names
  # initialize
  if(debug) browser()
  if(is.matrix(x)) {
    ncomp <- nrow(x)
    x <- array(t(x), dim = c(ndims,ncomp,1))
  } else {
    ncomp <- dim(x)[2]
    x <- aperm(x, perm = 3:1)
  }
  if(!is.matrix(theta)) theta <- matrix(theta, ncol = 1) else theta <- t(theta)
  nreps <- max(dim(x)[3], ncol(theta))
  if(dim(x)[3] == 1) x <- array(x, dim = c(ndims,ncomp,nreps))
  if(ncol(theta) == 1) theta <- matrix(theta, nparams, nreps)
  if(!dim(x)[3] == ncol(theta)) {
    stop("x and theta must have the same number of samples.")
  }
  if(length(dt) == 1) dt <- rep(dt, ncomp-1)
  if(length(dt) != ncomp-1) stop("Incorrectly specified dt.")

  ans <- model$loglik(x = as.double(x),
                      theta = as.double(theta), deltaT = as.double(dt),
                      nCompData = as.integer(ncomp), nReps = as.integer(nreps))

  ans
}
