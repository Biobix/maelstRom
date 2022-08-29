#' Calculate median imprinting.
#'
#' \code{median_imprinting} calculates median imprinting of allelic counts for a locus. 
#' It does this by ordering the samples according to degree of heterozygosity (i.e. minor over major allele count), 
#' then calculating the location of the "median heterozyogous sample" assuming HWE given input parameters
#' (\code{allelefreq}, \code{inbr}) as \code{round(number of samples \* allelefreq \* (1 - allelefreq) \* (1 - inbr))}, 
#' 2 times this unrounded number being the expected number of heterozygotes. Of this median heterozygous observation,
#' the degree of imprining is returned as \code{2 \* (0.5 - minor over major allele count)}.
#'
#' @param ref_counts Numeric vector. Reference counts.
#' @param var_counts Numeric vector. Variant counts.
#' @param allelefreq Number. Allele frequency.
#' @param inbr Number. Inbreeding coefficient (default = 0).
#' @return Median degree of imprinting of \code{ref_counts} and \code{var_counts} taking into account \code{allelefreq} and \code{inbr}.
#' @export
#' @examples
#' median_imprinting(c(5, 8, 10, 3, 5, 6, 23), c(8, 8, 6, 4, 4, 10, 0), 0.5, 0.12)
#' median_imprinting(c(5, 0, 0, 3, 5, 1, 23), c(1, 8, 6, 2, 0, 10, 0), 0.5, 0.12)
#' median_imprinting(c(5, 0, 0, 3, 5, 1, 23), c(1, 8, 6, 2, 0, 10, 0), 0.5)

median_imprinting <- function(ref_counts, var_counts, allelefreq, inbr = 0) {
  median_hetero <- round(length(ref_counts) * allelefreq * (1 - allelefreq) * (1 - inbr))
  ratios <- sort((apply(data.frame(ref_counts, var_counts), 1, min) / apply(data.frame(ref_counts, var_counts), 1, max)), TRUE)

  med_impr <- 2 * (0.5 - ratios[median_hetero])
  if (length(med_impr) == 0) med_impr <- 0
  return(med_impr)
}
