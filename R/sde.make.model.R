sde.make.model <-
function(list = NULL, ndims, nparams, data.names, param.names, custom.names,
                           sdeDr, sdeDf, isValidData, isValidParams, logCustomPrior,
                           cpp.out = TRUE, ..., debug = FALSE) {
  # list supplies values not otherwise specified
  nm <- c("ndims", "nparams", "data.names", "param.names", "custom.names",
          "sdeDr", "sdeDf", "isValidData", "isValidParams", "logCustomPrior")
  for(ii in nm) {
    eval(substitute(if(missing(ii) && !is.null(list$ii))
                    ii <- list$ii, list(ii = as.symbol(ii))))
  }
  # default data and parameter names
  if(missing(data.names)) data.names <- paste0("X", 1:ndims)
  if(missing(param.names)) param.names <- paste0("theta", 1:nparams)
  if(length(data.names) != ndims) stop("Incorrect data.names.")
  if(length(param.names) != nparams) stop("Incorrect param.names.")
  # check for custom prior
  if(missing(logCustomPrior)) logCustomPrior <- NULL
  if(missing(custom.names)) custom.names <- NULL
  sde.model <- list(ndims = ndims, nparams = nparams,
                    data.names = data.names, param.names = param.names,
                    sdeDr = sdeDr, sdeDf = sdeDf, isValidData = isValidData, isValidParams = isValidParams)
  if(!is.null(logCustomPrior)) {
    sde.model <- c(sde.model, list(logCustomPrior = logCustomPrior))
    if(!is.null(custom.names)) sde.model <- c(sde.model, list(custom.names = custom.names))
  } else {
    logCustomPrior <- "double CustomPrior::logPrior(double params[], double x[]) {
  return(0.0);
}"
  }
  # create c++ code
  h.nm <- c("Priors.h", "sdeCore.h", "ConvenienceFunctions.h", "sdeAPI-Rcpp.h")
  cpp.nm <- c("ConvenienceFunctions.cpp", "sdeCore.cpp", "Priors.cpp", "sdeAPI-Rcpp.cpp")
  hLines <- sapply(h.nm, function(nm) {
    con <- file(file.path(.msdeCppPath, nm), "r")
    zz <- readLines(con)
    close(con)
    zz
  })
  cppLines <- sapply(cpp.nm, function(nm) {
    con <- file(file.path(.msdeCppPath, nm), "r")
    zz <- readLines(con)
    close(con)
    zz
  })
  userLines <- c(sdeDr, sdeDf, isValidData, isValidParams, logCustomPrior)
  cpp.code <- c(paste0("const int nParams = ", nparams, ";"),
                paste0("const int nDims = ", ndims, ";"),
                unlist(hLines), userLines, unlist(cppLines),"\n")
  if(cpp.out) sde.model <- c(sde.model, list(cpp.code = cpp.code))

  # compile c++ code
  if(debug) browser()
  sourceCpp(code = paste(cpp.code, collapse = "\n"), env = environment(), ...)
  environment(sde.model$sim) <- globalenv()
  environment(sde.model$post) <- globalenv()
  environment(sde.model$drift) <- globalenv()
  environment(sde.model$diff) <- globalenv()
  environment(sde.model$loglik) <- globalenv()

  class(sde.model) <- "sde.model"
  sde.model
}
