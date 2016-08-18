#'@name sde.drift
#'@title Drift Function
#'@description Applies the drift function for \code{model} to \code{x} using \code{theta}
#'@param model An object of the sde.model class generated using \code{sde.make.model}
#'@param x A matrix of data
#'@param theta A matrix of parameters
#'@return A matrix containing the result of the drift function
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
#'dr <- sde.drift(model = hest.model, x = hsim$data, theta = theta)
#'@export
sde.drift <- function(model, x, theta) {
  if(class(model) != "sde.model")
    stop("Expecting object of class sde.model.  Use sde.make.model to create.")
  # model constants
  ndims <- model$ndims
  data.names <- model$data.names
  nparams <- model$nparams
  param.names <- model$param.names
  # initialize
  if(!is.matrix(x)) x <- matrix(x, ncol = 1) else x <- t(x)
  if(!is.matrix(theta)) theta <- matrix(theta, ncol = 1) else theta <- t(theta)
  nreps <- max(ncol(x), ncol(theta))
  if(ncol(x) == 1) x <- matrix(x, ndims, nreps)
  if(ncol(theta) == 1) theta <- matrix(theta, nparams, nreps)
  if(ncol(x) != ncol(theta)) {
    stop("x and theta must have the same number of rows.")
  }

  ans <- model$drift(x = as.double(x), theta = as.double(theta),
                     nReps = as.integer(nreps))

  dr <- matrix(ans, nreps, ndims, byrow = TRUE)
  dr
}
