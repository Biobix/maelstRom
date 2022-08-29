#' Logistic regression on the degree of heterozygosity between cases and controls
#'
#' \code{LOItest_logreg} performs logistic regression on the degree of heterozygosity in cases compared to control samples.
#'     This tests whether there is a significant difference (p-value) between the fraction of heterozygous samples in
#'     cases versus controls, which (in the case of a higher fraction in cases) can indicate loss of imprinting.
#'
#' @param ref_counts_ctrl,var_counts_ctrl Numeric vectors. Reference and variant counts of control samples.
#' @param ref_counts_case,var_counts_case Numeric vectors. Reference and variant counts of case samples.
#' @export
#' @return Beta-value (beta) and p-value (p.value) of the generalised linear model comparing the fraction of heterozygous samples
#'     in cases versus controls.
#' @examples
#' LOItest_logreg(c(5, 8, 10, 3, 5, 6, 23), c(8, 8, 6, 4, 4, 10, 0),
#'     c(5, 0, 0, 3, 5, 1, 23), c(1, 8, 6, 2, 0, 10, 0))

LOItest_logreg <- function(ref_counts_ctrl, var_counts_ctrl, ref_counts_case, var_counts_case) {
  counts_c <- data.frame("ref" = ref_counts_ctrl, "var" = var_counts_ctrl)
  counts_c$least <- apply(counts_c[, 1:2], 1, min)
  counts_c$most <- apply(counts_c[, 1:2], 1, max)
  counts_t <- data.frame("ref" = ref_counts_case, "var" = var_counts_case)
  counts_t$least <- apply(counts_t[, 1:2], 1, min)
  counts_t$most <- apply(counts_t[, 1:2], 1, max)
  
  tumc <- counts_t[, c("least", "most")]
  conc <- counts_c[, c("least", "most")]
  resp <- as.matrix(rbind(tumc, conc))
  status <- as.factor(c(rep("T", nrow(tumc)), rep("C", nrow(conc))))
  tempres <- stats::glm(resp ~ status, family = binomial)
  results <- list("beta" = coef(summary(tempres))[2, 1], "p.value" = anova(tempres, test="LRT")$`Pr(>Chi)`[2])
  
  return(results)
}

