#' Probability Mass Function of the beta-binomial mixture distribution modeling population-level RNAseq data
#'
#' @description \code{pmf_betabinomMix} calculates the probability of observing given population-level RNAseq data (i.e. both reference- and variant counts of one or more samples)
#' assuming a beta-binomial mixture model with parameter values as determined by the input. More specifically, the formula used is (using MAGE's \code{dBetaBinom} function):
#'
#' @description pr * dBetaBinom(ref_counts, ref_counts + var_counts, pi = 1 - SE, theta = theta_hom, LOG = FALSE) + 
#' 
#' @description pv * dBetaBinom(var_counts, ref_counts + var_counts, pi = 1 - SE, theta = theta_hom, LOG = FALSE) + 
#' 
#' @description prv * dBetaBinom(ref_counts, ref_counts + var_counts, pi = probshift, theta = theta_het, LOG = FALSE)
#'
#' @param ref_counts Number or Numeric vector Reference count(s).
#' @param var_counts Number or Numeric vector. Variant count(s).
#' @param probshift Number. The reference allele fraction in heterozygotes, indicating allelic bias when deviating from 0.5
#' @param SE Number. Sequencing error rate.
#' @param pr Number. Reference homozygote genotype probability of the locus.
#' @param pv Number. Variant homozygote genotype probability of the locus.
#' @param prv Number. Heterozygote genotype probability of the locus.
#' @param theta_hom Number. The dispersion parameter of the homozygous peaks.
#' @param theta_het Number. The dispersion parameter of the heterozygous peak.
#' @return Probability of observing \code{ref_counts} and \code{var_counts}
#' @export

pmf_betabinomMix <- function(ref_counts, var_counts, probshift, SE, pr, pv, prv, theta_hom, theta_het) {
  
  dmix <- pr * dBetaBinom(ref_counts, ref_counts + var_counts, pi = 1 - SE, theta = theta_hom, LOG = FALSE) + pv * dBetaBinom(var_counts, ref_counts + var_counts, pi = 1 - SE, theta = theta_hom, LOG = FALSE) + prv * dBetaBinom(ref_counts, ref_counts + var_counts, pi = probshift, theta = theta_het, LOG = FALSE)
  return(dmix)
}
