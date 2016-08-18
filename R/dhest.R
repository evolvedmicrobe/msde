#'@name dhest
#'@title Stationary Distribution of Heston Model Log-Return Differences
#'@description these will roughly be centered at alpha*dT \cr
#'             requires Fourier inversion of the characteristic function, ie. numerical \cr
#'             complex integration.  currently done by grid method. \cr
#'             LConst allows fewer exp evaluations if alpha is held constant. \cr
#'@param x
#'@param alpha
#'@param gamma
#'@param mu
#'@param beta
#'@param sigma
#'@param rho
#'@param dT
#'@param inter
#'@param n
#'@param LConst
#'@return numeric
#'@export
dhest <- function(x, alpha, gamma, mu, beta, sigma, rho, dT,
                  inter = c(-5e3, 5e3), n = 1e3, LConst, debug = FALSE) {
  if(missing(mu)) mu <- (beta + sigma^2/2)/(2*gamma)
  y <- x - alpha*dT
  px <- seq(inter[1], inter[2], len = n)
  dx <- px[2]-px[1]
  if(debug) browser()

  G <- gamma + 1i*rho*sigma*px
  O <- sqrt(G^2 + sigma^2*(px^2 - 1i*px))
  K <- O*dT/2
  ans <- (1 + exp(-2*K))/2 + (O^2-G^2+2*gamma*G)/(2*gamma*O) * (1 - exp(-2*K))/2
  ans <- gamma*mu/sigma^2*G*dT - 2*gamma*mu/sigma^2 * (K + log(ans))
  ans <- if(missing(LConst)) exp(1i*(px %o% y) + ans) else LConst * exp(ans)
  ans <- colSums(Re(ans))*dx/2/pi
}
