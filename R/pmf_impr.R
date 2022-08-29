#' Probability Mass Function for an imprinted SNP.
#'
#' @description \code{pmf_impr} calculates the probability of observing population-level RNAseq data (reference- and variant allele counts per sample)
#' assuming an imprinted binomial mixture model with the mean reference allele fraction in heterozygotes fixed at 0.5 (no allelic bias).
#' More specifically, the formula used is (though under-the-hood not using R's built-in dbinom function but our own mathematical implementation):
#' 
#' @description pr * dbinom(x = ref_counts, size = ref_counts + var_counts, prob = 1-SE) +
#' 
#' @description 0.5 * prv * dbinom(x = ref_counts, size = ref_counts + var_counts, prob = (0.5-(impr/2)) / (1-(impr/2)) * (1-SE) + (0.5/1-(impr/2)) * SE) +
#' 
#' @description 0.5 * prv * dbinom(x = ref_counts, size = ref_counts + var_counts, prob = (0.5-(impr/2)) / (1-(impr/2)) * SE + (0.5/1-(impr/2)) * (1-SE)) +
#' 
#' @description pv * dbinom(x = ref_counts, size = ref_counts + var_counts, prob = SE)
#' 
#' @description pr, pv and prv are calculated assuming Hardy-Weinberg-Equilibrium with given input parameters (inbreeding coefficient \code{inbr} and reference allele frequency \code{allelefreq})
#'
#' @param ref_counts Number or Numeric vector Reference count(s).
#' @param var_counts Number or Numeric vector. Variant count(s).
#' @param allelefreq Number. Allele frequency.
#' @param impr Number. Degree of imprinting.
#' @param SE Number. Sequencing error rate.
#' @param inbr Number. Degree of inbreeding (default = 0).
#' @return Probability of observing \code{ref_counts} and \code{var_counts}
#' @export
#' @examples
#' pmf_impr(10, 10, 0.5, 0, 0.002, 0.12)
#' pmf_impr(10, 10, 0.5, 1, 0.002)
#' pmf_impr(0, 10, 0.5, 1, 0.002, 0.12)
#' pmf_impr(c(5, 8, 10, 3, 5, 6, 23), c(8, 8, 6, 4, 4, 10, 0), 0.5, 1, 0.002, 0.12)
#' pmf_impr(c(5, 0, 0, 3, 5, 1, 23), c(1, 8, 6, 2, 0, 10, 0), 0.5, 1, 0.002)

pmf_impr <- function(ref_counts, var_counts, allelefreq, impr, SE, inbr = 0) {
  pr <- allelefreq ^ 2 + (inbr * allelefreq * (1 - allelefreq))
  pv <- (1 - allelefreq) ^ 2 + (inbr * allelefreq * (1 - allelefreq))
  prv <- allelefreq * (1 - allelefreq) * (1 - inbr)

  x_m <- MAGE::multinomial_coeff(ref_counts, var_counts)

  c1 <- ((0.5 - impr / 2) / (1 - impr / 2)) * (1 - SE) + (0.5 / (1 - impr / 2)) * SE
  c2 <- (0.5 / (1 - impr / 2)) * (1 - SE) + ((0.5 - impr / 2) / (1 - impr / 2)) * SE

  #homo
  m_RR <- x_m * (1 - SE) ^ ref_counts * SE ^ var_counts
  m_VV <- x_m * SE ^ ref_counts * (1 - SE) ^ var_counts
  #hetero
  m_RV <- x_m * c1 ^ ref_counts * c2 ^ var_counts
  m_VR <- x_m * c2 ^ ref_counts * c1 ^ var_counts

  PMF <- as.numeric(pr * m_RR + pv * m_VV + prv * m_RV + prv * m_VR)
  return(PMF)
}
