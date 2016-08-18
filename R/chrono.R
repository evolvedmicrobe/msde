chrono <-
function(start, display = TRUE) {
  if(missing(start)) ans <- proc.time()
  else {
    ans <- proc.time() - start
    ans <- as.numeric(ans[3])
    if(display) message("Execution time: ", round(ans, 2), " seconds.")
  }
  ans
}
