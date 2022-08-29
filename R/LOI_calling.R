#' Call LOI samples over al imprinted loci.
#'
#' \code{LOI_calling} calls LOI samples based on binomial p-value.
#'
#' @param data_hash Hash. Hash of SNP positions with a data frame for every position with at least columns of allelic
#'     counts ("A", "T", "C", "G"), reference and variants counts ("ref_count" and "var_count", respectively) and
#'     sample ids ("samples").
#' @param samples_all Character list. List of all sample names.
#' @param SE Number. Sequencing error rate.
#' @return Data frame with LOI calls per samples. The data frame contains sample ids ("samples"), number of least and
#'    most expressed alleles ("least" and "most", respectively), a quality check ("warning"), binomial p-values ("binom"),
#'    FDR adjusted p-values ("FDR") and LOI vs not LOI calls ("geno"). A "!" in the warning column indicates that more
#'    non-standard alleles than standard alleles were present for a sample (with a count of at least 2).

LOI_calling <- function(data_hash, samples_all, SE) {
  data_all <- data.frame("samples" = samples_all, "least" = rep(0, length(samples_all)), "most" = rep(0, length(samples_all)), "warning" = rep("", length(samples_all)))
  for (z in hash::keys(data_hash)) {
    ref <- data_hash[[z]]$ref_count
    var <- data_hash[[z]]$var_count
    errors <- data_hash[[z]]$A + data_hash[[z]]$T + data_hash[[z]]$C + data_hash[[z]]$G - data_hash[[z]]$ref_count - data_hash[[z]]$var_count
    warn_errors <- ifelse((errors > data_hash[[z]]$ref_count & errors > data_hash[[z]]$var_count & errors > 2), "!", "")

    counts <- data.frame("ref" = data_hash[[z]]$ref_count, "var" = data_hash[[z]]$var_count, "warning" = warn_errors, "samples" = data_hash[[z]]$sample, stringsAsFactors = FALSE)
    counts <- counts[match(data_all$samples, counts$samples), ]
    counts$samples <- data_all$samples
    counts[is.na(counts)] <- 0
    counts[which(counts$warning == 0), "warning"] <- ""
    counts$least <- apply(counts[, 1:2], 1, min)
    counts$most <- apply(counts[, 1:2], 1, max)
    data_all$least <- data_all$least + counts$least
    data_all$most <- data_all$most + counts$most
    data_all$warning <- ifelse((counts$warning == "!" | data_all$warning == "!"), "!", "")
  }

  data_zero <- data_all[which(data_all[, "most"] == 0 & data_all[, "least"] == 0), ]
  data_zero$binom <- rep("NA", nrow(data_zero))
  data_zero$geno <- rep("NA", nrow(data_zero))
  data_zero$FDR <- rep("NA", nrow(data_zero))

  if (nrow(data_zero) > 0) data_all <- data_all[-which(data_all[, "most"] == 0 & data_all[, "least"] == 0), ]

  data_all$binom <- pbinom(data_all[, "least"] - 1, data_all[, "least"] + data_all[, "most"], prob = SE, lower.tail = FALSE)
  p_values <- data_all$binom
  data_all$FDR <- p.adjust(p_values, method = "BH")
  data_all$geno <- ifelse(data_all$FDR <= 0.05, "LOI", "not LOI")
  if (nrow(data_zero) > 0) data_all <- rbind(data_all, data_zero)

  return(data_all)
}
