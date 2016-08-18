#' @name msde-package
#' @docType package
#' @aliases msde
#' @title Bayesian Inference for Multivariate Stochastic Differential Equations
#' @description MCMC sampling of the posterior distribution of the parameters of a multivariate SDE model given discretely observed data.  The user-defined models in msde are compiled in C++ which results in considerably faster posterior sampling than if the models were written in R.
#' @importFrom Rcpp sourceCpp
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @examples
#' # create the sde.model object for Heston's stochastic volatility model
#' hest.model <- sde.make.model(list = hestList, cpp.out = TRUE)
#'
#' # simulate data
#'
#' # model parameters
#' alpha <- .1
#' gamma <- 5.07
#' mu <- .05
#' sigma <- .48
#' rho <- -.77
#' beta <- 2*gamma*mu - sigma^2/2
#' theta <- c(alpha = alpha, gamma = gamma, beta = beta, sigma = sigma, rho = rho)
#' Y0 <- c(X = log(100), Z = 2*sqrt(mu)) # initial value
#' dT <- 1/252 # financial convention for one observation per day
#' ndays <- 501
#' # simulation
#' hest.sim <- sde.sim(model = hest.model, init.data = Y0, params = theta,
#'                     dt = dT, dt.sim = dT/1e2, N = ndays, burn = 100)
#'
#'
#' # parameter inference with Lebesgue prior
#' # basic Euler approximation with no missing data
#' init <- sde.init(data = hest.sim$data, dt = hest.sim$dt, m = 0,
#'                  par.index = 2, params = hest.sim$params)
#' nsamples <- 5e3
#' # standard deviations for random walk proposals
#' rw.jump.sd <- c(.2, 2, .2, .1, .1)
#' hest.post <- sde.post(model = hest.model, init = init, prior = NULL,
#'                       nsamples = nsamples, data.out.ind = data.out.ind,
#'                       rw.jump.sd = rw.jump.sd)
#'
#' # plot the posterior densities
#' theta.names <- expression(alpha, gamma, beta, sigma, rho)
#' par(mfrow = c(2,3), mar = c(4, 4.5, .5, .5))
#' for(ii in 1:5) {
#'   hist(hest.post$params[,ii], breaks = 100, freq = FALSE,
#'        xlab = theta.names[ii], main = "")
#'   abline(v = theta[ii], col = "red")
#' }
#' legend(x = "topright", legend = expression(hat(p)[0](theta*" | "*Y),
#'        theta), fill = c("black", "red"))
NULL
