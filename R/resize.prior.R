resize.prior <-
function(prior, new.size) {
  prior <- prior[sort(names(prior))]
  mv.names <- sort(c("Mu", "V"))
  gcop.names <- sort(c("dens.x", "dx", "rx", "dens.y", "ldens.y", "Dens.y",
                       "Rho", "mean", "sd"))
  if(identical(names(prior), mv.names)) prior.type <- "normal"
  if(identical(names(prior), gcop.names)) prior.type <- "gcop"
  if(!prior.type %in% c("normal", "gcop")) stop("prior must be either normal or gcop.")
  if(prior.type == "normal") {
    nold <- length(prior$Mu)
    if(nold > new.size) {
      prior$Mu <- prior$Mu[1:new.size]
      prior$V <- prior$V[1:nold,1:new.size]
    } else if(nold < new.size) {
      nnew <- new.size - nold
      prior$Mu <- c(prior$Mu, rep(0, nnew))
      V <- diag(new.size)
      V[1:nold,1:nold] <- prior$V
      prior$V <- V
    }
  }
  if(prior.type == "gcop") {
    nold <- length(gcop$mean)
  }
}
