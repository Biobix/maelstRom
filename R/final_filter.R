#' Final filtering of MAGE analysis results and writing to output files
#'
#' \code{final_filter} only retains significantly imprinted SNPs (after adjusting for multiple testing) and
#' SNPs of interest (with suitable GOF and degree of (median) imprinting) over all chromosomes. Results and
#' allelic count files are generated. When both file_all_counts and file_impr_counts are set to FALSE,
#' this function can be used to simply filter the results_df input.
#'
#' @param data_hash Hash. Hash of SNP positions with a data frame for every SNP position.
#' @param results_df Data frame. Results data frame with columns: "position", "gene", "LRT", "p", "estimated.i",
#'     "allele.frequency", "dbSNP", "reference", "variant", "est_SE", "coverage", "nr_samples", "GOF", "symmetry",
#'     "med_impr", est_inbreeding", "tot_inbreeding".
#' @param results_wd String. Directory where results files are written to.
#' @param gof_filt Number. Minimal Goodness of Fit, which is the mean(log(sample likelihood under imprinted model * sample coverage + 1)) across samples of a locus. A good (and default) cutoff is 0.8.
#' @param adj_p_filt Number. The FDR adjusted singnificance level filter (default is 0.05).
#' @param med_impr_filt Number. Minimal median imprinting (default is 0.8).
#' @param i_filt Number. Minimal degree of imprinting (default is 0.6).
#' @param file_all Logical. Should a file with all SNP information (imprinted and non-imprinted SNPs) be made (default is TRUE).
#' @param file_impr Logical. Should a file with imprinted SNP information be made (default is TRUE).
#' @param file_all_counts Logical. Should a file with all SNP counts (imprinted and non-imprinted SNPs) be made (default is FALSE).
#' @param file_impr_counts Logical. Should a file with imprinted SNP counts be made (default is TRUE).
#' @export
#' @return Data frame with results filtered on adjusted p-value, GOF, median imprinting and degree of imprinting.

final_filter <- function(data_hash, results_df, results_wd, gof_filt = 1.2, adj_p_filt = 0.05, med_impr_filt = 0.8, i_filt = 0.6, file_all = TRUE, file_impr = TRUE, file_all_counts = FALSE, file_impr_counts = TRUE) {
  #WRITE ALL RESULTS TO FILE CHR_IMPRINTING_ALL.TXT
  if (file_all) {
    for (chr in unique(results_df$chromosome)) {
      f <- file(paste(results_wd, "chr", chr, "_results_imprinting_all.txt", sep = ""), "w")
      out <- write.table(results_df[which(results_df$chromosome == chr), ], f, row.names = FALSE, sep = "\t")
      close(f)
    }
  }

  #WRITE ALL ALLELECOUNTS TO FILE CHR_ALLELECOUNTS_ALL.TXT
  if (file_all_counts) {
    for (chr in unique(results_df$chromosome)) {
      results_chr <- results_df[which(results_df$chromosome == chr), ]
      f2 <- file(paste(results_wd, "chr", chr, "_allelecounts_all.txt", sep = ""), "w")
      writeLines(paste("position", "dbSNP", "ref", "var", "sample_id", "sample_nr", sep = "\t"), f2)
      for (z in as.character(results_chr$dbSNP)) {
        out <- writeLines(paste(data_hash[[z]]$position, z, data_hash[[z]]$ref_count, data_hash[[z]]$var_count, data_hash[[z]]$sample, data_hash[[z]]$sample_nr, sep = "\t"), f2)
      }
      close(f2)
    }
  }

  results_df <- results_df[which(results_df$GOF > gof_filt), ]
  if(nrow(results_df)==0){
    warning("No results remain after filtering")
    return(results_df)
  }

  #ADJUST FOR MULTIPLE TESTING WITH BEJAMINI_HOCHBERG FDR
  p_values <- results_df$p
  adjusted_p <- p.adjust(p_values, method = "BH")
  results_df$adj_p <- adjusted_p

  #FILTER OUT SAMPLES ON ADJUSTED P-VALUE, MEDIAN IMPRINTING, ESTIMATED DEGREE OF IMPRINTING
  results_df <- results_df[which(results_df$adj_p <= adj_p_filt & results_df$med_impr >= med_impr_filt & results_df$estimated.i >= i_filt), ]
  if(nrow(results_df)==0){
    warning("No results remain after filtering")
    return(results_df)
  }

  #WRITE ALL RESULTS TO FILE IMPRINTED_GENES.TXT
  if (file_impr) {
    f <- file(paste(results_wd, "imprinted_genes.txt", sep = ""), "w")
    out <- write.table(results_df, f, row.names = FALSE, sep = "\t")
    close(f)
  }

  #WRITE IMPRINTED ALLELECOUNTS TO FILE CHR_ALLELECOUNTS_IMPRINTED.TXT
  if (file_impr_counts) {
    for (chr in unique(results_df$chromosome)) {
      results_chr <- results_df[which(results_df$chromosome == chr), ]
      f2 <- file(paste(results_wd, "chr", chr, "_allelecounts_imprinted.txt", sep = ""), "w")
      writeLines(paste("position", "dbSNP", "ref", "var", "sample_id", "sample_nr", sep = "\t"), f2)
      for (z in as.character(results_chr$dbSNP)) {
        out <- writeLines(paste(data_hash[[z]]$position, z, data_hash[[z]]$ref_count, data_hash[[z]]$var_count, data_hash[[z]]$sample, data_hash[[z]]$sample_nr, sep = "\t"), f2)
      }
      close(f2)
    }
  }

  return(results_df)
}
