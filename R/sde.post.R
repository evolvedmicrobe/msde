sde.post <-
function(model,
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
