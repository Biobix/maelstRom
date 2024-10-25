logSums_MaxMethod <- function(logvec){
  
  # EXPECTED: a vector of logs
  logvec <- logvec[!is.na(logvec)]
  X1 <- max(logvec)
  
  H <- logvec - X1
  
  return(X1 + log(sum(exp(H))))
  
}
