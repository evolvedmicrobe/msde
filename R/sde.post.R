#'@name sde.post
#'@title Posterior Inference
#'@description Markov Chain Monte Carlo sampling of posterior distribution
#'@param model an object of the sde.model class generated using \code{sde.make.model}
#'@param init A list containing values for any combination of \code{init.data}, \code{init.params}, \code{par.index} and \code{dt}.
#'            Values provided in this list will be used instead of values provided as parameters to \code{sde.post}
#'@param init.data Data for the model
#'@param init.params Parameters for the model
#'@param par.index
#'@param dt Delta-t for the model
#'@param nsamples An integer giving the number of samples to generate with the Markov chain after the required
#'                number of samples have been burned
#'@param burn An integer giving the number of samples to burn before beginning to store samples
#'@param data.out.ind
#'@param prior A list containing the parameters for the model's prior function. \cr
#'             The parameters must have the correct names for the prior indicated by \code{prior.type} \cr
#'             normal: Mu, V \cr
#'             gcop: dens.x, dx, rx, dens.y, ldens.y, Dens.y, Rho, mean, sd \cr
#'             flat: flat prior does not accept parameters. Parameter list is ignored \cr
#'             custom: no argument checking for custom priors
#'@param prior.type A string denoting which prior should be used for Posterior Inference (\code{"default"} by default).\cr
#'                  Must be one of: default, flat, normal, gcop or custom\cr
#'                  if \code{"flat"} a flat prior will be used\cr
#'                  if \code{"normal"} a normal prior will be used\cr
#'                  if \code{"gcop"} a Gaussian Copula prior will be used\cr
#'                  if \code{"default"} the names of the parameters provided in \code{prior} will be used to determine
#'                      whether a normal or gcop prior should be used.
#'                      If the parameter list does not match either of these, a flat prior will be used\cr
#'                  if \code{"custom"} the custom prior that was passed to \code{sde.make.model} when \code{model} was created
#'                      will be used
#'@param rw.jump.sd A numeric vector containing random walk sd's (default \code{NULL})
#'@param update.data A boolean (default \code{TRUE})
#'@param update.params A boolean (default \code{TRUE})
#'@param update.multi A boolean (default \code{TRUE})
#'@param smp.old the result of a previous call to \code{sde.post}. Used for multi update
#'@param max.multi.tries An integer (default \code{5e3})
#'@param log.multi.acc A boolean. If \code{TRUE} and \code{update.multi == TRUE} log.multi.acc will be included in the output of this function (default \code{FALSE})
#'@param loglik.out A boolean. If \code{TRUE} loglik will be included in the output of this function (default \code{FALSE})
#'@param last.miss.out A boolean. If \code{TRUE} last.miss will be included in the output of this function (default \code{FALSE})
#'@param verbose A boolean that determines the amount of text printed out by \code{sde.post}. (default \code{TRUE})
#'@param debug a boolean (\code{FALSE} by default) if set to \code{TRUE}, will cause the function to open a browser mid-call
#'@return a list containing:\cr \cr
#'                          params, \cr
#'                          data, \cr
#'                          loglik (only if \code{loglik.out == TRUE}), \cr
#'                          dt, \cr
#'                          par.index, \cr
#'                          data.out.ind, \cr
#'                          init.data, \cr
#'                          init.params, \cr
#'                          rw.jump.sd, \cr
#'                          prior, \cr
#'                          prior.type, \cr
#'                          last.iter, \cr
#'                          last.miss (only if \code{last.miss.out == TRUE})\cr
#'                          multi.ind (only if \code{update.multi == TRUE})\cr
#'                          log.multi.acc (only if \code{update.multi == TRUE && log.multi.acc == TRUE})\cr
#'                          accept\cr
#'                          time
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
#'# posterior inference
#'
#'# parameters of normal prior
#'norm.prior <- c(theta, Y0[2])
#'norm.prior <- list(Mu = norm.prior,
#'                   V = diag(abs(norm.prior)/2) %*% f(length(norm.prior)) %*% diag(abs(norm.prior)/2))
#'
#'k <- 1
#'par.index <- 1
#'init <- sde.init(data = hsim$data, dt = dT, k = k, par.index = par.index, params = theta)
#'
#'# reduce the prior on N random variables to n <= N
#'n <- length(theta) + hest.model$ndim - init$par.index[1]
#'norm.prior$Mu <- norm.prior$Mu[1:n]
#'norm.prior$V <- norm.prior$V[1:n,1:n]
#'
#'nsamples <- 10
#'burn <- 2
#'rw.jump.sd <- rep(.1, 6)
#'
#'hpost <- sde.post(model = hest.model, init = init, nsamples = nsamples, burn = burn,
#'                  data.out.ind = nsamples, rw.jump.sd = rw.jump.sd,
#'                  prior = norm.prior)
#'
#'
#'# and a different model
#'gl05.model <- sde.make.model(list = gl05List, model.name = "gl05")
#'
#'theta <- log(c(.1,.7, .35, .2, .1, .9, .3, .1))
#'names(theta) <- paste("c", 1:8, sep = "")
#'
#'Y0 <- c(R = 5, P = 5, Q = 5, D = 5)
#'
#'# simulation
#'N <- 10
#'burn <- 10
#'dT <- 1
#'
#'gsim <- sde.sim(model = gl05.model, init.data = Y0, params = theta, dt = dT, dt.sim = dT/1e5,
#'                N = N, burn = burn, nreps = 1)
#'
#'
#'# inference
#'
#'# prior
#'norm.prior <- c(theta, Y0[-1])
#'norm.prior <- list(Mu = norm.prior,
#'                   V = diag(abs(norm.prior)/2) %*% f(length(norm.prior)) %*% diag(abs(norm.prior)/2))
#'psim <- rmvnorm(1e5, mean = norm.prior$Mu*1.1, sigma = norm.prior$V*2)
#'gcop.prior <- cop.par(X = psim, n = 512)
#'
#'k <- 2 # or 1, 2, 3, ...
#'par.index <- sample(1:gl05.model$ndims, N, replace = TRUE) # 1,2,3, or 4, but not 0
#'init <- sde.init(data = gsim$data, dt = dT, k = k, par.index = par.index, params = theta)
#'
#'n <- length(theta) + gl05.model$ndim - init$par.index[1]
#'gcop.prior$dens.x <- gcop.prior$dens.x[1:n]
#'gcop.prior$dens.y <- gcop.prior$dens.y[1:n]
#'gcop.prior$Dens.y <- gcop.prior$Dens.y[1:n]
#'gcop.prior$dx <- gcop.prior$dx[1:n]
#'gcop.prior$ldens.y <- gcop.prior$ldens.y[1:n]
#'gcop.prior$mean <- gcop.prior$mean[1:n]
#'gcop.prior$Rho <- gcop.prior$Rho[1:n,1:n]
#'gcop.prior$rx <- gcop.prior$rx[1:(2*n)]
#'gcop.prior$sd <- gcop.prior$sd[1:n]
#'
#'nsamples <- 10
#'burn <- 2
#'rw.jump.sd <- rep(.1, gl05.model$ndims+gl05.model$nparams-1)
#'
#'gpost <- sde.post(model = gl05.model, init = init, nsamples = nsamples, burn = burn,
#'                  data.out.ind = nsamples, rw.jump.sd = rw.jump.sd,
#'                  prior = gcop.prior)
#'@export
sde.post <- function(model,
                     init, init.data, init.params, par.index, dt, nsamples, burn,
                     data.out.ind, prior,
                     rw.jump.sd = NULL,
                     update.data = TRUE, update.params = TRUE,
                     update.multi = TRUE, smp.old,
                     max.multi.tries = 5e3, log.multi.acc = FALSE,
                     loglik.out = FALSE, last.miss.out = FALSE,
                     verbose = TRUE, debug = FALSE) {
  if(class(model) != "sde.model")
    stop("Expecting object of class sde.model.  Use sde.make.model to create.")
  # model constants
  ndims <- model$ndims
  data.names <- model$data.names
  nparams <- model$nparams
  param.names <- model$param.names
  if(!is.null(model$logCustomPrior)) {
    has.custom.prior <- TRUE
    custom.names <- model$custom.names
    if(!is.null(custom.names)) custom.names <- sort(custom.names)
  } else {
    has.custom.prior <- FALSE
    custom.names <- NULL
  }
  Prior.Types <- c("flat", "normal", "gcop", "custom") # allowed prior types
  # initial parameters
  if(!missing(init)) {
    if(!is.null(init$data)) init.data <- init$data
    if(!is.null(init$dt)) dt <- init$dt
    if(!is.null(init$par.index)) par.index <- init$par.index
    if(!is.null(init$params)) init.params <- init$params
  }
  # parse inputs
  ncomp <- nrow(init.data)
  nmiss0 <- ndims - par.index[1]
  nparams2 <- nparams+nmiss0
  if(ndims != ncol(init.data))
    stop("init.data does not have the right number of components.")
  if(!is.null(colnames(init.data))) {
    if(any(colnames(init.data) != data.names))
      stop("Incorrect data.names.")
  }
  if(nparams != length(init.params)) stop("init.params does not have the right length.")
  if(!is.null(colnames(init.params))) {
    if(any(colnames(init.params) != param.names))
      stop("Incorrect param.names.")
  }
  # time
  if(length(dt) == 1) dt <- rep(dt, ncomp-1)
  if(length(dt) != ncomp-1) stop("Incorrectly specified dt.")
  if(missing(burn)) burn <- .1
  if(burn < 1) burn <- nsamples*burn
  burn <- floor(burn)
  # storage
  if(missing(data.out.ind)) data.out.ind <- 2e3
  if(!is.list(data.out.ind)) {
    data.out.row <- data.out.ind
    data.out.col <- 1:ncomp
  } else {
    data.out.row <- data.out.ind$row
    data.out.col <- data.out.ind$col
  }
  if(length(data.out.row) == 1)
    data.out.row <- unique(floor(seq(1, nsamples, len = data.out.row)))
  if(is.logical(data.out.row)) data.out.row <- which(data.out.row)
  nsamples.out <- length(data.out.row)
  if(is.logical(data.out.col)) data.out.col <- which(data.out.col)
  ncomp.out <- length(data.out.col)
  data.out.ind <- list(row = data.out.row, col = data.out.col)
  if(all(par.index == ndims)) update.data <- FALSE
  if(update.data) {
    ndata.out <- nsamples.out*ndims*ncomp.out
  } else {
    ndata.out <- 1
    data.out.ind$row <- 1:nsamples
  }
  nparams.out <- ifelse(update.params, nsamples*nparams, 1)
  # prior specification
  mv.names <- sort(c("Mu", "V"))
  mv.dim <- c(Mu = nparams2, V = nparams2^2)
  gcop.names <- sort(c("dens.x", "dx", "rx", "dens.y", "ldens.y", "Dens.y",
                       "Rho", "mean", "sd"))
  #gcop.dim <- c(rep(nparams2, 6), nparams2^2, 2*nparams2, nparams2)
  gcop.dim <- c(dens.x = nparams2, dx = nparams2, rx = 2*nparams2,
                dens.y = nparams2, ldens.y = nparams2, Dens.y = nparams2,
                Rho = nparams2^2, mean = nparams2, sd = nparams2)
  prior.names <- names(prior)
  if(!is.null(prior.names)) prior.names <- sort(prior.names)
  prior2 <- prior
  if(debug) browser()
  if(is.null(prior)) {
    # flat prior
    prior.type <- "flat"
    prior <- list()
  } else if (identical(prior.names, mv.names)) {
    # mv prior
    prior.type <- "normal"
    if(any(sapply(prior, length)[prior.names] != mv.dim[prior.names]))
      stop("Incorrect prior specification (needs nparams + nmiss0 prior variables).")
    prior$V <- chol(prior$V)
  } else if(identical(prior.names, gcop.names)) {
    # gcop prior
    prior.type <- "gcop"
    if(any(sapply(prior, length)[prior.names] != gcop.dim[prior.names]))
      stop("Incorrect prior specification (needs nparams + nmiss0 prior variables).")
    prior$Rho <- chol(prior$Rho)
    prior$nbreaks <- sapply(prior$dens.x, length)
    prior <- sapply(prior, unlist)
  } else {
    # custom prior
    prior.type <- "custom"
    if(!has.custom.prior) stop("prior is not recognized default and no custom prior supplied.")
    if(!is.null(custom.names) && !identical(prior.names, custom.names))
      stop("Supplied custom prior names don't match those defined by sde.model.")
    prior <- sapply(prior, unlist, recursive = TRUE)
  }
  if(is.null(rw.jump.sd)) {
    rw.jump.sd <- abs(init.params)/4
    if(nmiss0 > 0) rw.jump.sd <- c(rw.jump.sd, abs(init.data[par.index[1]+(1:nmiss0)])/4)
  }
  if(length(rw.jump.sd) < nparams + nmiss0)
    stop("Incorrectly specified random walk jump sd's (need at least nparams + nmiss0).")
  # last missing output
  nmissN <- ndims - par.index[ncomp]
  if(nmissN == 0) last.miss.out <- FALSE
  # multiresolution
  if(missing(smp.old)) update.multi <- FALSE
  if(update.multi) {
    if(debug) browser()
    past.data <- smp.old$data
    past.params <- smp.old$params
    past.dt <- smp.old$dt
    past.prior <- smp.old$prior
    # past.prior.type <- smp.old$prior.type
    past.ncomp <- length(past.dt) + 1
    if(is.matrix(past.data)) {
      past.nsamples <- nrow(past.params)
      past.data <- array(past.data, dim = c(past.ncomp, ndims, past.nsamples))
      past.data <- aperm(past.data, perm = c(2,1,3))
    } else {
      past.nsamples <- length(smp.old$data.out.ind$row)
      past.data <- aperm(past.data, perm = c(3,2,1))
    }
    if(length(smp.old$data.out.ind$col) != past.ncomp)
      stop("Multiresolution requires complete data in smp.old.")
    if(!is.matrix(past.params)) {
      past.params <- matrix(past.params, nparams, past.nsamples)
    } else {
      past.params <- t(past.params[smp.old$data.out.ind$row,])
    }
    # multi.ind
    cdt <- cumsum(c(0, dt))
    past.cdt <- cumsum(c(0, past.dt))
    multi.ind <- sapply(cdt, function(x) {
      tmp <- which(past.cdt <= x)
      c(past.cdt[tmp[length(tmp)]], past.cdt[which(past.cdt >= x)[1]])
    })
    multi.ind <- which(multi.ind[1,] != multi.ind[2,])
    if(any(diff(multi.ind) <= 1))
      stop("Invalid multiresolution setting (max. 1 new data pt between old pts).")
    nmulti <- length(multi.ind)
    past.multi.ind <- (1:ncomp)[-multi.ind]
    nmulti.out <- ifelse(log.multi.acc, nmulti*nsamples.out, 1)
    # past prior
    past.prior.names <- names(past.prior)
    if(!is.null(past.prior)) past.prior.names <- sort(past.prior.names)
    if(is.null(past.prior)) {
      # flat prior
      past.prior.type <- "flat"
      past.prior <- list()
    }
    else if (identical(past.prior.names, mv.names)) {
      # mv prior
      past.prior.type <- "normal"
      if(any(sapply(past.prior, length)[past.prior.names] != mv.dim))
        stop("Incorrect prior specification (needs nparams + nmiss0 prior variables).")
      past.prior$V <- chol(past.prior$V)
    } else if(identical(past.prior.names, gcop.names)) {
      # gcop prior
      past.prior.type <- "gcop"
      if(any(sapply(past.prior, length)[past.prior.names] != gcop.dim))
        stop("Incorrect prior specification (needs nparams + nmiss0 prior variables).")
      past.prior$Rho <- chol(past.prior$Rho)
      past.prior$nbreaks <- sapply(past.prior$dens.x, length)
      past.prior <- sapply(past.prior, unlist)
    } else {
      # custom prior
      past.prior.type <- "custom"
      if(!has.custom.prior) stop("prior is not recognized default and no custom prior supplied.")
      if(!is.null(custom.names) && !identical(past.prior.names, custom.names))
        stop("Supplied custom prior names don't match those defined by sde.model.")
      past.prior <- sapply(past.prior, unlist, recursive = TRUE)
    }
  } else {
    log.multi.acc <- FALSE
    past.nsamples <- 1
    past.ncomp <- 1
    past.data <- 0
    past.params <- 0
    multi.ind <- 2
    past.multi.ind <- 0
    nmulti <- 1
    nmulti.out <- 1
    past.prior <- list(gcop = FALSE)
    past.prior.type <- "flat"
  }
  if(verbose) {
    message("Output size:")
    if(update.params) message("params = ", round(nparams.out, 2))
    if(update.data) message("data = ", round(ndata.out, 2))
    if(update.multi & log.multi.acc)
      message("multi.acc = ", round(nmulti.out, 2))
    if(verbose > 1) {
      ans <- readline("Press Q to quit, any other key to proceed: ")
      if(substr(ans, 1, 1) == "Q") {
        message("Ended by user.")
        return()
      }
    }
    message("Running posterior sampler...")
  }
  if(debug) browser()
  tm <- chrono()

  ans <- model$post(nParamsOut = as.integer(nparams.out),
                    nDataOut = as.integer(ndata.out),
                    initParams = as.double(init.params),
                    initData = as.double(t(init.data)),
                    deltaT = as.double(dt),
                    nObsDim = as.integer(par.index),
                    nCompData = as.integer(c(ncomp, ncomp.out)),
                    nSamples = as.integer(nsamples),
                    burn = as.integer(burn),
                    dataOutRow = as.integer(data.out.row-1),
                    dataOutCol = as.integer(data.out.col-1),
                    updateParams = as.double(update.params),
                    updateData = as.double(update.data),
                    updateMulti = as.double(update.multi),
                    rwJumpSd = as.double(rw.jump.sd),
                    priorType = which(prior.type == Prior.Types),
                    priorParams = prior,
                    multiIndex = as.integer(c(multi.ind, past.multi.ind)-1),
                    nMultiInd = as.integer(nmulti),
                    past_Params = as.double(past.params),
                    past_Data = as.double(past.data),
                    pastPriorType = which(past.prior.type == Prior.Types),
                    pastPriorParams = past.prior,
                    nMaxMultiTries = as.integer(max.multi.tries),
                    nPast_Samples = as.integer(past.nsamples),
                    nPastCompData = as.integer(past.ncomp),
                    updateLogMulti = as.integer(log.multi.acc),
                    nLogMultiOut = as.integer(nmulti.out),
                    updateLogLik = as.integer(loglik.out),
                    nLogLikOut = as.integer(ifelse(loglik.out, nsamples, 1)),
                    updateLastMiss = as.integer(last.miss.out),
                    nLastMissOut = as.integer(ifelse(last.miss.out, nsamples*nmissN, 1)))

  tm <- chrono(tm, display = verbose)
  names(ans) <- c("paramsOut", "dataOut", "logMultiOut",
                  "paramAccept", "gibbsAccept", "multiAccept",
                  "logLikOut", "lastMissOut", "lastIter")
  # acceptance rates
  if(debug) browser()
  accept <- NULL
  if(update.data) {
    accept <- c(accept, list(data = ans$gibbsAccept/(nsamples+burn)))
    if(verbose) {
      message("Gibbs accept: ", signif(mean(accept$data[par.index < ndims])*100,3), "%")
    }
  }
  if(update.params) {
    accept <- c(accept, list(params = ans$paramAccept[1:nparams]/(nsamples+burn)))
    if(verbose) {
      for(ii in 1:nparams)
        message(param.names[ii], " accept: ", signif(accept$params[ii]*100,3), "%")
    }
  }
  if(update.data && (nmiss0 > 0)) {
    accept <- c(accept,
                list(miss0 = ans$paramAccept[nparams+(1:nmiss0)]/(nsamples+burn)))
    for(ii in 1:nmiss0) {
      message("first ", data.names[par.index[1] + ii], " accept: ",
              signif(accept$miss0[ii]*100,3), "%")
    }
  }
  if(update.multi) {
    accept <- c(accept, list(multi = ans$multiAccept/(nsamples+burn)))
    if(verbose) message("multi accept: ", signif(accept$multi*100,3), "%")
  }
  out <- list()
  if(update.params) {
    out <- c(out, list(params = matrix(ans$paramsOut,
                         ncol = nparams, byrow = TRUE,
                         dimnames = list(NULL, param.names))))
  } else out <- c(out, list(params = init.params))
  if(update.data) {
    out <- c(out, list(data = aperm(array(ans$dataOut,
                         dim = c(ndims, ncomp.out, nsamples.out),
                         dimnames = list(data.names, NULL, NULL)), perm = 3:1)))
  } else out <- c(out, list(data = init.data))
  if(loglik.out) out <- c(out, list(loglik = ans$logLikOut))
  out <- c(out, list(dt = dt, par.index = par.index, data.out.ind = data.out.ind,
                     init.data = init.data, init.params = init.params,
                     rw.jump.sd = rw.jump.sd, prior = prior2))
  last.iter <- list(params = ans$lastIter[1:nparams],
                    data = matrix(ans$lastIter[nparams + 1:(ncomp*ndims)],
                      ncomp, ndims, byrow = TRUE))
  names(last.iter$params) <- param.names
  colnames(last.iter$data) <- data.names
  out <- c(out, list(last.iter = last.iter))
  if(last.miss.out) {
    last.miss <- matrix(ans$lastMissOut, nsamples, nmissN, byrow = TRUE)
    colnames(last.miss) <- data.names[par.index[ncomp]+(1:nmissN)]
    out <- c(out, list(last.miss = last.miss))
  }
  if(update.multi) {
    out <- c(out, list(multi.ind = multi.ind))
    if(log.multi.acc)
      out <- c(out,
               list(log.multi.acc = t(matrix(ans$logMultiOut, nrow = nmulti))))
  }
  out <- c(out, list(accept = accept))
  out
}
