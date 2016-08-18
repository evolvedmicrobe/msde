#'@name dgcir
#'@title gcir stationary distribution (regular scale)
#'@param x
#'@param gamma
#'@param mu
#'@param sigma
#'@param lambda
#'@param log
#'@param debug a boolean (\code{FALSE} by default) if set to \code{TRUE}, will cause the function to open a browser mid-call
#'@return a numeric vector
#'@export
dgcir <- function(x, gamma, mu, sigma, lambda, log = FALSE, debug = FALSE) {
  if(debug) browser()
  a2 <- 2*lambda
  b <- 2*gamma/sigma^2
  if(lambda != .5) {
    ans <- -2*log(sigma) -a2*log(x) - b * (x^(2-a2)/(2-a2) - mu*x^(1-a2)/(1-a2))
  } else {
    ans <- -2*log(sigma) + (b*mu-1)*log(x) - b*x
  }
  if(!log) ans <- exp(ans)
  ans
}
