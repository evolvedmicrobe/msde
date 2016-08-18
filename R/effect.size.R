#'@name effect.size
#'@title effective sample size calculation
#'@description calculates the effective sample size for either the output of an MCMC sampler, or the weights of an importance sampler.
#'@param chain the output of an MCMC sampler for a single coordinate.
#'@param weights the weights for an importance sampler.
#'@param prop logical.  \code{prop = FALSE} returns the effective number of samples, whereas \code{prop = TRUE} returns it as a fraction.
#'@return A scalar.
#'@examples
#'# Gibbs sampler for the bivariate Normal
#'rho <- .9
#'nsamples <- 1e3
#'X <- matrix(NA, nsamples, 2)
#'colnames(X) <- c("X1", "X2")
#'x <- c(0,0)
#'for(ii in 1:nsamples) {
#'  x[1] <- rnorm(1, mean = rho*x[2], sd = sqrt(1-rho^2))
#'  x[2] <- rnorm(1, mean = rho*x[1], sd = sqrt(1-rho^2))
#'  X[ii,] <- x
#'}
#'
#'# effective number of samples
#'apply(X, 2, effect.size, prop = FALSE)
#'
#'# Importance sampler with independent Normals as the proposal distribution
#'Xprop <- rmvnorm(nsamples, sigma = diag(2))
#'# calculate weights
#'lprop <- dmvnorm(Xprop, sigma = diag(2), log = TRUE)
#'ltarg <- dmvnorm(Xprop, sigma = cbind(c(1,rho), c(rho,1)), log = TRUE)
#'wgt <- ltarg - lprop
#'wgt <- exp(wgt - max(wgt))
#'wgt <- wgt/sum(wgt)
#'# resampling of proposals to approximate the target
#'Xtarg <- Xprop[sample(x = 1:nsamples, size = nsamples, replace = TRUE,  prob = wgt),]
#'
#'# proportion of effective samples
#'effect.size(weights = wgt, prop = TRUE)
#'@export
effect.size <- function(chain = NULL, weights = NULL, prop = TRUE) {
  # inefficiency factor, basically the sum of the autocorrelations
  ineff <- function(x) {
    a <- invisible(acf(x, plot = FALSE))
    if(any(is.na(a$acf))) return(0)
    2*sum(a$acf)-1
  }
  if(!is.null(chain)) {
    N <- length(chain)
    # eff <- 1/(ineff(chain) + 1)
    eff <- 1/ineff(chain)
  } else {
    N <- length(weights)
    eff <- 1/(1+var(weights)/mean(weights)^2)
  }
  ifelse(prop, eff, N*eff)
}
