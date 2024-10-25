#' Goodness-of-fit test comparing data following a uniform distribution to a Hardy-Weinberg distribution.
#'
#' \code{GOF_uniform} determines the GOF of the heterozygous samples to the HWE model.
#'
#' @param ref_counts Numeric list. Reference counts.
#' @param var_counts Numeric list. Variant counts.
#' @param allelefreq Number. Allele frequency.
#' @param SE Number. Sequencing error rate.
#' @param inbr Number. Degree of inbreeding (default = 0).
#' @return Geometric mean of the ratios of the putative heterozygous samples as a measure of the GOF.
#' @export
#' @examples
#' GOF_uniform(c(5, 8, 10, 3, 5, 6, 23), c(8, 8, 6, 4, 4, 10, 0), 0.5, 0.002, 0.12)
#' GOF_uniform(c(5, 0, 0, 3, 5, 1, 23), c(1, 8, 6, 2, 0, 10, 0), 0.5, 0.002)

GOF_uniform <- function(ref_counts, var_counts, allelefreq, SE, inbr = 0) {
  prv <- 2 * allelefreq * (1 - allelefreq) * (1 - inbr)
  pr <- allelefreq ^ 2 + inbr * allelefreq * (1 - allelefreq)
  pv <- (1 - allelefreq) ^ 2 + inbr * allelefreq * (1 - allelefreq)

  PMF_uniform <- pr * dbinom(ref_counts, ref_counts + var_counts, prob = 1 - SE) +
    pv * dbinom(var_counts, ref_counts + var_counts, prob = 1 - SE) +
    (prv / (ref_counts + var_counts + 1))

  PMF_HWE <- maelstRom::pmf_impr(ref_counts, var_counts, allelefreq = allelefreq, impr = 0, SE = SE, inbr = inbr)

  hetero_frac <- round(prv * length(ref_counts))
  ratio <- PMF_uniform / PMF_HWE
  ratio_sorted <- sort(ratio)
  RT_hetero <- prod(ratio_sorted[1:hetero_frac]) ^ (1 / hetero_frac)

  return(RT_hetero)
}
