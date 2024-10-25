# #' gradient of the beta-binomial log-likelihood function with respect to the pi-parameter
# #'
# #' \code{grad_pi} calculates the value of the gradient of the beta-binomial log-likelihood function with respect to pi at given observations & parameters.
# #' The function makes use of different implementations of this gradient,
# #' aiming to produce correct results even on extreme data, and this in a speedy fashion
# #'
# #' @param ms Numeric vector. Vector of number of successes
# #' @param ns Numeric vector. Vector of number of trials
# #' @param pi Number. Probability of success; \code{0 >= pi >= 1}
# #' @param theta Number. Dispersion parameter; \code{0 >= theta > +Inf}
# #' @param MemLim Integer. The memory limit, in bits, of numbers used during the calculation of the gradient.
# #'    For extreme data/parameter values, it may be necesssary to go beyond R's \code{double} memory, i.e. 53 bits,
# #'    in order to get correct results. However, this argument prevents boundless memory usage.
# #'    In case this limit is reached, \code{grad_pi} falls back on an implementation assuming a theta-parameter of zero,
# #'    as this should only happen for extremely close-to-zero theta-values anyway.
# #'    A warning message is generated if this occurs.
# #' @param Xtra Integer. An internal control parameter, determines the number of decimal places that are to be
# #'    stored correctly in memory during calculation of the gradient, upon which the number of
# #'    bits to be used per number depends. It's not recommended to change this from the default.
# #' @export
# #' @return A vector of values of the gradient, with the same length as \code{ms} and \code{ns}.

#grad_pi <- function(ms, ns, pi, theta, Xtra = 7, MemLim = 2048){
  
#  nsC <- ns
#  msC <- ms
#  nLow <- which(ns <= 0)
#  nHigh <- which(ns > 0)
  
#  if(length(nLow)>0){
#    ns <- nsC[nLow]
#    ms <- msC[nLow]
#    out1 <- grad_pi_old(ms, ns, pi, theta)
#  }
  
#  if(length(nHigh)>0){
    
#    ns <- nsC[nHigh]
#    ms <- msC[nHigh]
    
#    Hlp <- log10((1/theta)+max(ns))
#    LogParRat <- ceiling(Hlp + log10(Hlp*log(10)))
#    NecBits <- ceiling((LogParRat+Xtra)/log10(2))
    
#    if(NecBits <= 53 & theta !=0){
#      PV1 <- ifelse(ms>0, (1/theta) * (digamma(ms+(pi/theta)) - digamma(pi/theta)), 0)
#      PV2 <- ifelse((ns-ms)>0, (1/theta) * (digamma((ns-ms)+((1-pi)/theta)) - digamma((1-pi)/theta)), 0)
#      out2 <- PV1 - PV2
#    } else if(theta == 0 | NecBits > MemLim){
#      if(theta!=0){
#        warning("memory limit reached in grad_pi", call. = FALSE)
#      }
#      out2 <- ms/pi - (ns-ms)/(1-pi)
#    } else{
#      out2 <- grad_pi_MP(ms, ns, pi, theta, LogParRat+Xtra)
#    }
    
#  }
  
#  out <- rep(0, length(nsC))
#  if(length(nLow)>0){
#    out[nLow] <- out1
#  }
#  if(length(nHigh)>0){
#    out[nHigh] <- out2
#  }
#  return(out)
#}
