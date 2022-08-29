#' @name dpqrBetaBinom
#' @title The Beta-binomial distribution
#'
#' @description Density, distribution function, quantile function and random generation for the beta-binomial distribution with parameters \code{pi} and \code{theta}.
#' These can be interpreted as an expected probability of success and a dispersion parameter respectively,
#' and can respectively be expressed in terms of the common \code{alpha, beta} parameterization as \code{pi = alpha/(alpha+beta)} and \code{theta = 1/(alpha+beta)}.
#'
#' @param ms Numeric vector. Vector of number of successes
#' @param ns Numeric vector. Vector of number of trials
#' @param pi Number. Probability of success; \code{0 >= pi >= 1}
#' @param theta Number. Dispersion parameter; \code{0 >= theta > +Inf}
#' @param MemLim Integer. The memory limit, in bits, of numbers used during the calculation of the density.
#'    For extreme data/parameter values, it may be necesssary to go beyond R's \code{double} memory, i.e. 53 bits,
#'    in order to get correct results. However, this argument prevents boundless memory usage.
#'    In case this limit is reached, \code{dBetaBinom} falls back on the regular \code{dbinom} function,
#'    as this should only happen for extremely close-to-zero theta-values anyway.
#'    A warning message is generated if this occurs.
#' @param Xtra Integer. An internal control parameter, determines the number of decimal places that are to be
#'    stored correctly in memory during calculation of beta-binomial densities, upon which the number of
#'    bits to be used per number depends. It's not recommended to change this from the default.
#' @param LOG Logical. If TRUE, \code{dBetaBinom} returns log-densities.
#' @param log.p Logical. If TRUE, probabilities p returned and accepted by \code{pBetaBinom} and \code{qBetaBinom} respectively are given as log(p).
#' @return \code{dBetaBinom} gives the density, \code{pBetaBinom} gives the distribution function, \code{qBetaBinom} gives the quantile function and \code{rBetaBinom} generates random deviates. Where applicable, output has the same length as \code{ms} and/or \code{ns}.

NULL

#' @rdname dpqrBetaBinom
#' @export
dBetaBinom <- function(ms, ns, pi, theta, MemLim = 2048, Xtra = 7, LOG = FALSE) {
  
  nsC <- ns
  msC <- ms
  nLow <- which(ns <= 50)
  nHigh <- which(ns > 50)
  
  if(length(nLow)>0){ # We can use a specific beta-binomial implementation that is faster for low-count data
    ns <- nsC[nLow]
    ms <- msC[nLow]
    out1 <- dBetaBinom_cpp_old(ms, ns, pi, theta, LOG)
  }
  
  if(length(nHigh)>0){
    
    ns <- nsC[nHigh]
    ms <- msC[nHigh]
    
    Hlp <- log10((1/theta)+max(ns))
    LogParRat <- ceiling(Hlp + log10(Hlp*log(10)))
    Hlp2 <- ceiling(log10(abs(lbeta(a = max(ns), b = 1/theta))))
    NecBits <- ceiling((Hlp2+Xtra)/log10(2))
    
    NecBitsTRUE <- ceiling((LogParRat+Xtra)/log10(2)) # The previous necessary bits are calculated assuming beta-binomial calculation happen using lbeta-functions.
    # In the case of a large requirement in number of bits, however, we need to work in C++, which doesn't offer the lbeta function.
    # It does, however, offer an alternative computation using lgamma-functions,
    # in which case this expression is an approximation of the necessary number of bits.
    
    if(NecBits <= 53 & theta !=0){ # Necessary bits <= 53? Then an R double's default precision is enough
      PV1 <- ifelse(ms>0, -lbeta(a = ms, b = pi/theta) - log(ms), 0)
      PV2 <- ifelse((ns-ms)>0, -lbeta(a = (ns-ms), b = ((1-pi)/theta)) - log(ns-ms), 0)
      PV3 <- ifelse(ns>1, lbeta(a = ns, b = 1/theta) + log(ns), 0)
      out <- PV1 + PV2 + PV3
      if (LOG) {
        out2 <- out
      } else {
        out2 <- exp(out)
      }
    } else if ((theta == 0 | NecBitsTRUE > MemLim)){ # If theta==0, we're dealing with the regular binomial
      if(theta!=0){
        warning("memory limit reached in dBetaBinom\n", call. = FALSE)
      }
      out2 <- dbinom(ms, ns, pi, log = LOG)
    } else{
      out2 <- dBetaBinom_MP(ms, ns, pi, theta, LOG=LOG, NecPres=LogParRat+Xtra)
    }
    
  }
  
  out <- rep(0, length(nsC))
  if(length(nLow)>0){
    out[nLow] <- out1
  }
  if(length(nHigh)>0){
    out[nHigh] <- out2
  }
  return(out)
}


#' @rdname dpqrBetaBinom
#' @export
pBetaBinom <- function(ms, ns, pi, theta, lower.tail = TRUE, log.p = FALSE){
  out <- c()
  for(i in 1:length(ms)){
    out <- c(out, sum(dBetaBinom(ms=0:ms[i], ns=rep(ns[i], ms[i]+1), pi, theta, LOG=FALSE)))
  }
  
  if(log.p){
    if(lower.tail){
      return(log(out))
    }else{
      return(log(1-out))
    }
  }else{
    if(lower.tail){
      return(out)
    }else{
      return(1-out)
    }
  }
  
}


#' @rdname dpqrBetaBinom
#' @export
qBetaBinom <- function(p, ns, pi, theta, lower.tail = TRUE, log.p = FALSE){
  
  if(log.p){
    p <- exp(p)
  }
  
  if(!lower.tail){
    p <- 1-p
  }
  
  out <- c()
  for(i in 1:length(p)){
    CDF <- cumsum(dBetaBinom(ms=0:ns[i], ns=rep(ns[i], ns[i]+1), pi, theta, LOG=FALSE))
    CDF[length(CDF)] <- 1
    IND <- which(CDF >= min(p[i], 1))[1]-1
    out <- c(out, IND)
  }
  return(out)
}


#' @rdname dpqrBetaBinom
#' @export
rBetaBinom <- function(ns, pi, theta){
  out <- c()
  for(i in 1:length(ns)){
    P <- runif(1)
    CuVec <- cumsum(dBetaBinom(ms=0:ns[i], ns=rep(ns[i], ns[i]+1), pi, theta, LOG=FALSE))
    IND <- which(CuVec >= P)[1]-1
    out <- c(out, IND)
  }
  return(out)
}
