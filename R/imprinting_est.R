#' Estimate the degree of imprinting.
#'
#' \code{imprinting_est} estimates the degree of imprinting using iterative likelihood ratio tests for a specific SNP,
#' going from imprinting = 0 to imprinting = 1 in steps of 0.01, performing LRTs along the way (against imprinting = 0), and finally retaining the best fit.
#' See the function \code{pmf_impr} for what the degree of imprinting entails, mathematically
#'
#' @param ref_counts Numeric list. Reference counts.
#' @param var_counts Numeric list. Variant counts.
#' @param allelefreq Number. Allele frequency.
#' @param SE Number. Sequencing error rate.
#' @param inbr Number. Degree of inbreeding (default = 0).
#' @export
#' @return A list containing the following components:
#' \item{est_i}{The estimated degree of imprinting.}
#' \item{LRT}{The test statistic of the likelihood ratio test against no imprinting.}
#' \item{p_value}{The p-value of the likelihood ratio test against no imprinting.}
#' \item{GOF_likelihood}{A goodness-of-fit value based on count-corrected likelihood, i.e. mean log PMF-value according to \code{pmf_impr} multiplied by each sample's coverage + 1 across samples of a locus.}
#' @examples
#' imprinting_est(c(5, 8, 10, 3, 5, 6, 23), c(8, 8, 6, 4, 4, 10, 0), 0.5, 0.002)
#' imprinting_est(c(5, 0, 0, 3, 5, 1, 23), c(1, 8, 6, 2, 0, 10, 0), 0.5, 0.002, 0.12)

imprinting_est <- function(ref_counts, var_counts, allelefreq, SE, inbr = 0) {
  #PMF1 => i = 0
  PMF1 <- maelstRom::pmf_impr(ref_counts, var_counts, allelefreq, 0, SE, inbr)

  #PMF2 => i = î
  i_est <- 0
  product_max <- -Inf
  PMF2 <- rep(0, length(ref_counts))
  for (i in seq(0, 1, 0.01)) {
    PMF2_iest <- maelstRom::pmf_impr(ref_counts, var_counts, allelefreq, i, SE, inbr)
    product_est <- sum(log(PMF2_iest))
    if (product_est > product_max) {
      product_max <- product_est
      i_est <- i
      PMF2 <- PMF2_iest
    }
  }

  #LRT: DETECTION OF IMPRINTING
  LRT <- - 2 * sum(log(PMF1)) + 2 * sum(log(PMF2))

  #DETERMINE P WITH CORRECTED NULL DISTRIBUTION
  p_impr <- 1 / 2 * pchisq(LRT, 1, lower.tail = FALSE) + 1 / 2 * pchisq(LRT, 0, lower.tail = FALSE)

  #GOF BASED ON CORRECTED LIKELIHOODS
  PMF_corrected <- maelstRom::pmf_impr(ref_counts, var_counts, allelefreq, i_est, SE, inbr) * (ref_counts + var_counts + 1)
  logLikelihood <- mean(log(PMF_corrected))

  results <- list("est_i" = i_est, "LRT" = LRT, "p_value" = p_impr, "GOF_likelihood" = logLikelihood)
  return(results)
}
