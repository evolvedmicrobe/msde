sv.init <-
function(x, dt, vblock.size, debug = FALSE) {
  y <- diff(x)
  N <- length(y)
  ind <- cbind(floor(seq(1, N-vblock.size+1, len = N)),
               ceiling(seq(vblock.size, N, len = N)))
  if(debug) browser()
  tmp <- apply(ind, 1, function(i) {
    m <- mean(y[i[1]:i[2]])
    v <- var(y[i[1]:i[2]])
    c(alpha = (m+v/2)/dt, v = v/dt)
  })
  v <- tmp[2,]
  list(alpha = mean(tmp[1,]), v = c(v[1], v))
}
