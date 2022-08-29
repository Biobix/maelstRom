#' This function is a work in progress. As such, it is not exported yet; but it allows for some leniency on the sequencing error parameter, hence "SESlack"

pmf_betabinomMix_SEslack <- function(ref_counts, var_counts, probshift, piRR, piVV, pr, pv, prv, theta_RR, theta_VV, theta_het) {
  
  dmix <- pr * dBetaBinom(ref_counts, ref_counts + var_counts, pi = piRR, theta = theta_RR, LOG = FALSE) + pv * dBetaBinom(var_counts, ref_counts + var_counts, pi = piVV, theta = theta_VV, LOG = FALSE) + prv * dBetaBinom(ref_counts, ref_counts + var_counts, pi = probshift, theta = theta_het, LOG = FALSE)
  return(dmix)
}
