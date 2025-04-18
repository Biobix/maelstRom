#' Returns an (approximate) log of a sum based on individual logs of the terms being summed
#'
#' \code{logSums_MaxMethod} calculates the (approximate) log of a sum of terms, based on the individual logs of these terms.
#' Simply taking the exponent of these individual terms, then summing them, then taking the log again, is usually not possible if taking the exponent of the
#' supplied values results in a numeric overflow. The rationale here is:
#' 
#' log(x + y + z) =(suppose x is the MAXIMUM amongst {x,y,z}) log(((x+y+z)/x) * x) = log(x) + log((x+y+z)/x) = log(x) + log(x/x + y/x + z/x) =
#' log(x) + log(exp(log(x/x)) + exp(log(y/x)) + exp(log(z/x))) = log(x) + log(exp(log(x)-log(x)) + exp(log(y)-log(x)) + exp(log(z)-log(x)))
#' 
#' Now, log(x)-log(x) cancels out, and the remaining (log(y)-log(x)) and (log(z)-log(x)) are negative (x is the maximum amontst {x,y,z}) meaning that,
#' worst case scenario, x is that much larger than y and z that these terms are high and negative, resulting in a calculated exp-value of zero;
#' this isn't too bad though, as the first term (log(x)) is still there, and if x is so much higher than the other terms, then the log of sums will be
#' mainly dictated by x anyway (when remaining on R's default numerical precision, and given a limited number of terms being summed; once again, this value is approximate).
#'
#' @param logvec Numeric vector. seperate log values, of which we want to calculate the log of the sum of their (non-log) values
#' @export
#' @return The log of the sum of the values for which the seperate logs were supplied as input.

logSums_MaxMethod <- function(logvec){
  
  # EXPECTED: a vector of logs
  logvec <- logvec[!is.na(logvec)]
  X1 <- max(logvec)
  
  H <- logvec - X1
  
  return(X1 + log(sum(exp(H))))
  
}
