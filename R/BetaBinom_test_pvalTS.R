#' Performs an exact beta-binomial p-test
#'
#' \code{BetaBinom_test_pvalTS} performs an exact test,
#' returning a p-value that reflects the probability of the observed data \code{m}
#' assuming a beta-binomial PMF with given parameters \code{m}, \code{pi} and \code{theta}.
#' See \code{dBetaBinom} for more information on the used parameterization.
#'
#' @param m Number. Number of successes
#' @param n Number. Number of trials
#' @param pi Number. Probability of success; \code{0 >= pi >= 1}
#' @param theta Number. Overdispersion parameter; \code{0 >= theta > +Inf}
#' @export
#' @return The probability to make an observation equally or less likely than the input data

BetaBinom_test_pvalTS <- function(m, n, pi, theta){
  # returns two-sided pval only
  
  if (n == 0){
    return(1)
  } else{
    HV <- dBetaBinom(ms = 1:n, ns = rep(n, n), pi = pi, theta = theta, LOG = FALSE)
    return(sum(HV[HV<=dBetaBinom(m, n, pi = pi, theta = theta, LOG = FALSE)]))
  }
  
}
