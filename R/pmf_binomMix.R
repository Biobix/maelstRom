#' Probability Mass Function of the binomial mixture distribution modeling population-level RNAseq data
#'
#' @description \code{pmf_binomMix} calculates the probability of observing given population-level RNAseq data (i.e. both reference- and variant counts of one or more samples)
#' assuming a binomial mixture model with parameter values as determined by the input. More specifically, the formula used is:
#' 
#' @description pr * dbinom(ref_counts, ref_counts + var_counts, prob = 1 - SE) + 
#' 
#' @description pv * dbinom(var_counts, ref_counts + var_counts, prob = 1 - SE) + 
#' 
#' @description prv * dbinom(ref_counts, ref_counts + var_counts, prob = probshift)
#'
#' @param ref_counts Number or Numeric vector Reference count(s).
#' @param var_counts Number or Numeric vector. Variant count(s).
#' @param probshift Number. The reference allele fraction in heterozygotes, indicating allelic bias when deviating from 0.5
#' @param SE Number. Sequencing error rate.
#' @param pr Number. Reference homozygote genotype probability of the locus.
#' @param pv Number. Variant homozygote genotype probability of the locus.
#' @param prv Number. Heterozygote genotype probability of the locus.
#' @return Probability of observing \code{ref_counts} and \code{var_counts}
#' @export
#' @examples
#' pmf_binomMix(10, 10, 0.5, 0.002, 0.25, 0.25, 0.5)
#' pmf_binomMix(0, 10, 0.8, 0.002, 0.25, 0.25, 0.5)
#' pmf_binomMix(c(5, 8, 10, 3, 5, 6, 23), c(8, 8, 6, 4, 4, 10, 0), 0.8, 0.002, 0.25, 0.25, 0.5)
#' pmf_binomMix(c(5, 0, 0, 3, 5, 1, 23), c(1, 8, 6, 2, 0, 10, 0), 0.8, 0.002, 0.25, 0.25, 0.5)

pmf_binomMix <- function(ref_counts, var_counts, probshift, SE, pr, pv, prv) {
  dmix <- pr * dbinom(ref_counts, ref_counts + var_counts, prob = 1 - SE) + pv * dbinom(var_counts, ref_counts + var_counts, prob = 1 - SE) + prv * dbinom(ref_counts, ref_counts + var_counts, prob = probshift)
  return(dmix)
}
