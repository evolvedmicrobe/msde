#'@name sde.sim
#'@title Simulation
#'@description Runs a simulation using the supplied model, data and parameters.
#'@param model an object of the sde.model class generated using \code{sde.make.model}
#'@param init.data A matrix containing data for the model.
#'                 Must have as many columns as the model has dimensions.
#'                 If \code{init.data} has only 1 row, that row will be repeated for each of the \code{nreps} Markov chains, otherwise
#'                 Must have as many rows as \code{nreps}
#'@param params A matrix containing parameters for the model
#'              Must have as many columns as the model has parameters
#'              If \code{params} has only 1 row, that row will be repeated for each of the \code{nreps} Markov chains, otherwise
#'              Must have as many rows as \code{nreps}
#'@param dt A double. The ratio of \code{dt/dt.sim} determines
#'                    how many intermediate steps there will be between each stored value
#'@param dt.sim A double The ratio of \code{dt/dt.sim} determines
#'                       how many intermediate steps there will be between each stored value
#'@param N An integer that determines how many data values will be stored for each Markov chain
#'@param burn An integer representing the number of data values to burn for each Markov chain (default \code{0})
#'            OR a double 0 < burn < 1 representing the number of data values to burn for each chain as a fraction of N
#'@param nreps An integer giving the number of Markov chains to use
#'@param max.bad.draws An integer giving the maximum number of times that invalid data can be encountered in the simulation (default \code{5e3})
#'@param verbose A boolean that determines how much text is printed out by this function (default \code{TRUE})
#'@param debug a boolean (\code{FALSE} by default) if set to \code{TRUE}, will cause the function to open a browser mid-call
#'@return a list containing: data, params, dt, dt.sim, nbad
#'@examples
#'# Create the model
#'hest.model <- sde.make.model(list = hestList)
#'
#'theta <- c(alpha = .1, gamma = 5, beta = .8, sigma = .6, rho = -.7) # parameters
#'Y0 <- c(X = log(100), Z = .1) # initial values
#'
#'# simulate data
#'N <- 100
#'burn <- 0
#'dT <- 1/252
#'
#'hsim <- sde.sim(model = hest.model, init.data = Y0, params = theta, dt = dT, dt.sim = dT/100,
#'                N = N, burn = burn, nreps = 1)
#'
#'
#'# plot
#'par(mfrow = c(1,2))
#'plot(hsim$data[,"X"], type = "l", xlab = "Time (days)", ylab = expression(X[t]))
#'plot(hsim$data[,"Z"]^2/4, type = "l", xlab = "Time (days)", ylab = expression(V[t]))
#'@export
sde.sim <- function(model,
                    init.data, params, dt, dt.sim,
                    N, burn = 0, nreps = 1,
                    max.bad.draws = 5e3, verbose = TRUE,
                    debug = FALSE) {
  if(class(model) != "sde.model")
    stop("Expecting object of class sde.model.  Use sde.make.model to create.")
  # model constants
  ndims <- model$ndims
  data.names <- model$data.names
  nparams <- model$nparams
  param.names <- model$param.names
  # initialize
  if(!is.matrix(init.data)) init.data <- t(as.matrix(init.data))
  if(!is.matrix(params)) params <- t(as.matrix(params))
  if(ndims != ncol(init.data))
    stop("init.data does not have the right number of components.")
  if(nparams != ncol(params))
    stop("params does not have the right length/number of columns.")
  # data
  if(!is.null(colnames(init.data))) {
    if(any(colnames(init.data) != data.names))
      stop("Incorrect data.names.")
  }
  if(nrow(init.data) == 1) {
    init.data <- matrix(init.data, nrow = nreps, ncol = ndims, byrow = TRUE)
  }
  if(nrow(init.data) != nreps) stop("init.data does not have the right number of rows.")
  colnames(init.data) <- data.names
  # params
  init.params <- params
  if(!is.null(colnames(init.params))) {
    if(any(colnames(init.params) != param.names))
      stop("Incorrect param.names.")
  }
  if(nrow(params) == 1) {
    params <- matrix(params, nrow = nreps, ncol = nparams, byrow = TRUE)
  }
  if(nrow(params) != nreps) stop("params does not have the right number of rows.")
  colnames(init.params) <- param.names
  # time
  if(dt.sim <= dt) {
    r <- ceiling(dt/dt.sim)
    t <- dt/r
  } else {
    r <- 1
    t <- dt
  }
  if(burn < 1) burn <- N*burn
  burn <- floor(burn)
  if(verbose) {
    message("Normal draws required: ", round((N+burn)*r*nreps, 2))
    if(verbose > 1) {
      ans <- readline("Press Q to quit, any other key to proceed: ")
      if(substr(ans, 1, 1) == "Q") {
        message("Ended by user.")
        return()
      }
    }
    message("Running simulation...")
  }
  if(debug) browser()
  tm <- chrono()

  ans <- model$sim(nDataOut = as.integer((N+burn)*ndims*nreps),
                   N = as.integer(N+burn),
                   reps = as.integer(nreps),
                   r = as.integer(r),
                   delta = as.double(t),
                   MAXBAD = as.integer(max.bad.draws),
                   initData = as.double(t(init.data)),
                   params = as.double(t(params)))

  tm <- chrono(tm, display = verbose)
  names(ans) <- c("dataOut", "nBadDraws")
  if(verbose) message("Bad Draws = ", ans$nBadDraws)
  data <- aperm(array(ans$dataOut, dim = c(ndims, N+burn, nreps)), perm = 3:1)
  dimnames(data) <- list(NULL, NULL, data.names)
  out <- list(data = data[,burn+1:N,], params = init.params[,], dt = dt, dt.sim = t,
              nbad = ans$nBadDraws)
  out
}
