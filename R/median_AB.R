#' Calculates the median allelic bias
#'
#' \code{median_AB} calculates median allelic bias of a given locus, i.e. the heterozygous reference allele ratio.
#' It does this by ordering the samples according to degree of heterozygosity (i.e. minor over major allele count), 
#' then calculating the location of the "median heterozyogous sample" assuming HWE given input parameters
#' (\code{allelefreq}, \code{inbr}) as \code{round(number of samples \* allelefreq \* (1 - allelefreq) \* (1 - inbr))}, 
#' 2 times this unrounded number being the expected number of heterozygotes. Of this median heterozygous observation,
#' the allelic bias is returned (reference over total allele ratio).
#' 
#' @param ref_counts Numeric vector. Reference counts.
#' @param var_counts Numeric vector. Variant counts.
#' @param allelefreq Number. Allele frequency.
#' @param inbr Number. Inbreeding coefficient (default = 0).
#' @return Median heterozygous allelic bias of \code{ref_counts} and \code{var_counts} taking into account
#' \code{allelefreq} and \code{inbr}.
#' @export
#' @examples
#' median_AB(c(5, 8, 10, 3, 5, 6, 23), c(8, 8, 6, 4, 4, 10, 0), 0.5, 0.12)
#' median_AB(c(5, 0, 0, 3, 5, 1, 23), c(1, 8, 6, 2, 0, 10, 0), 0.5, 0.12)
#' median_AB(c(5, 0, 0, 3, 5, 1, 23), c(1, 8, 6, 2, 0, 10, 0), 0.5)

median_AB <- function(ref_counts, var_counts, allelefreq, inbr = 0) {
  pqmin <- min(c(allelefreq, (1 - allelefreq)))
  F_alt <- -pqmin / (1 - pqmin)
  if (inbr < F_alt) {
    inbr <- F_alt
  }

  median_hetero <- round(length(ref_counts) * allelefreq * (1 - allelefreq) * (1 - inbr))
  ratios_hetero <- apply(data.frame(ref_counts, var_counts), 1, min) / apply(data.frame(ref_counts, var_counts), 1, max)
  df_ratios <- data.frame("ref" = ref_counts, "var" = var_counts, "hetero" = ratios_hetero, "ratio" = ref_counts / (var_counts + ref_counts))
  df_ratios <- df_ratios[order(-df_ratios$hetero), ]

  #med_ase <- (df_ratios[median_hetero, "ratio"] - 0.5) / 0.5
  med_ase <- df_ratios[median_hetero, "ratio"]
  if (length(med_ase) == 0) med_ase <- 0

  return(med_ase)
}
