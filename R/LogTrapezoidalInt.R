LogTrapezoidalInt <- function(f, lower, upper, n, ...){
  
  Xs <- seq(from = lower, to = upper, length.out = n)
  h <- (upper-lower)/(n-1)
  
  Ys <- f(Xs, ...)
  Ys[1] <- Ys[1] + log(5) - log(12)
  Ys[n] <- Ys[1] + log(5) - log(12)
  Ys[2] <- Ys[1] + log(13) - log(12)
  Ys[n-1] <- Ys[1] + log(13) - log(12)
  
  return(log(h)+logSums_MaxMethod(Ys))
  
}
