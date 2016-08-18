effect.size <-
function(chain = NULL, weights = NULL, prop = TRUE) {
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
