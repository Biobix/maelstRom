#' Goodness-of-fit test for symmetry.
#'
#' \code{symmetry_gof} calculates the chi-squared test for symmetry of allelic counts.
#' This tests whether the number of samples expressing more than 50% of the reference allele AND the amount of samples expressing less than 50% of the reference allele
#' correspond to the amounts we expect when assuming HWE AS WELL AS no allelic shift (on average equal expression of reference- and variant alleles in heterozygotes).
#' These are the expected number of homozygote ref-samples plus half the expected number of heterozygous samples,
#' and the expected number of homozygote var-samples plus half the expected number of heterozygous samples, both under HWE, respectively.
#' This function is used during imprinting analyses ONLY, as these are the assumptions made during those analyses.
#'
#' @param ref_counts Numeric list. Reference counts.
#' @param var_counts Numeric list. Variant counts.
#' @param allelefreq Number. Allele frequency.
#' @return The p-value of chi-squared test for \code{ref_counts} and \code{var_counts}.
#' @export
#' @examples
#' symmetry_gof(c(5, 8, 10, 3, 5, 6, 23, 4, 12, 10, 9, 6, 7, 25),
#' c(8, 8, 6, 4, 4, 10, 0, 7, 8, 4, 2, 7, 15, 13), 0.5)
#' symmetry_gof(c(5, 0, 0, 3, 5, 1, 23, 0, 12, 0, 9, 0, 1, 0),
#' c(1, 8, 6, 2, 0, 10, 0, 7, 0, 4, 0, 7, 0, 13), 0.5)

symmetry_gof <- function(ref_counts, var_counts, allelefreq) {
  fracA <- round(ref_counts / (ref_counts + var_counts), 2)
  obs <- c(length(which(fracA < 0.5)), length(which(fracA > 0.5)))
  exp <- c((1-allelefreq), allelefreq) * (length(ref_counts) - length(which(fracA == 0.5)))
  dat <- data.frame(obs, exp)
  chi2 <- chisq.test(dat)$p.value
  return(chi2)
}
