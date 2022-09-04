#' multi-locus wrapper function for imprinting detection and loss-of-imprinting analysis
#'
#' \code{impr_LOI_analysis} is a wrapper function of \code{symmetry_gof}, \code{imprinting_est}, \code{median_imprinting}, \code{final_filter}, and
#' \code{LOItest_logreg} performing MAGE's entire (loss of) imprinting pipeline. It goes over the following steps:
#' \enumerate{
#'   \item Prior filtering on basis of the \code{allelefreq_prel} column that should be present in the control dataframe,
#'   which it should be after running \code{AllelicMeta_est}.
#'   \item Prior filtering via \code{symmetry_gof}, see its documentation for more details.
#'   \item Running \code{imprinting_est} to detect imprinted loci and write results to \code{impr_res}; this includes a measure of median
#'   imprinting as returned by the function \code{median_imprinting} as well.
#'   \item Filter \code{impr_res} down to significantly and sufficiently imprinted loci using \code{final_filter}, then perform loss-of-imprinting
#'   detection using \code{LOItest_logreg}. LOI-detection results are written to \code{LOI_res}.
#' }
#' 
#' @param DataList A list containing two lists of dataframes, which in turn contain control- and case-data (in that order). Subsequent entries in both lists
#' correspond to subsequent loci (lists must be the same size and named per-locus). Each of the dataframes should at least contain both a \code{ref_count} and
#' \code{var_count} column; control dataframes should also contain the columns \code{allelefreq_prel} (from running \code{AllelicMeta_est}), \code{locus_id}, 
#' \code{ref}, \code{var}, and \code{coverage} (from running \code{prior_filter}).
#' @param SE Number. Sequencing error rate, population metaparameter.
#' @param inbr Number. Inbreeding coefficient, population metaparameter.
#' @param MinMinorAllelefreq Number. Minimal allowed minor allele frequency as determinde via the estimated reference allele frequency by \code{AllelicMeta_est} to consider the locus for imprinting analysis.
#' @param sym_filt Number. Minimum required symmetry statistic to consider the locus for imprinting analysis, see \code{symmetry_gof}.
#' @param adj_p_filt Number. Required FDR-corrected significance level of imprinting for loci to be considered for loss-of-imprinting detection, see the function \code{final_filter}. 
#' @param gof_filt Number. Required minimum goodness-of-fit heuristic of imprinted loci to be considered for loss-of-imprinting detection, see the function \code{final_filter}.
#' @param med_i_filt Number. Required minimum median imprinting of imprinted loci to be considered for loss-of-imprinting detection, see the functions \code{median_imprinting} and \code{final_filter}.
#' @param i_filt Number. Required minimum estimated imprinting of imprinted loci to be considered for loss-of-imprinting detection, see the functions \code{imprinting_est} and \code{final_filter}.
#' @export
#' @return A list containing the following components:
#' \item{impr_res}{Results of the imprinting detection analysis, containing, for every locus:}
#' \itemize{
#'      \item{position}{The locus' name, according to \code{names(DataList)}}
#'      \item{LRT}{The test statistic of the likelihood ratio test against no imprinting (see \code{imprinting_est}).}
#'      \item{p}{The p-value of the likelihood ratio test against no imprinting (see \code{imprinting_est}).}
#'      \item{estimated.i}{The estimated degree of imprinting (see \code{imprinting_est}).}
#'      \item{allele.frequency}{Reference allele frequency as estimated by \code{AllelicMeta_est}, already given as input.}
#'      \item{reference}{Reference allele nucleotide.}
#'      \item{variant}{Variant allele nucleotide.}
#'      \item{med_cov}{Median coverage (reference + variant) across the locus.}
#'      \item{nr_samples}{Number of samples (with a total reference + variant read count of at least 1) covering the locus.}
#'      \item{GOF}{Goodness-of-fit heuristic as determined by \code{imprinting_est}.}
#'      \item{symmetry}{Symmetry statistic as determined by \code{symmetry_gof}.}
#'      \item{med_impr}{Median imprinting as determined by \code{median_imprinting.}}
#'    }
#' \item{LOI_res}{A dataframe containing, only for loci fitting well to MAGE's imprinting model (see \code{imprinting_est}, determined by \code{gof_filt}) 
#' that are significantly (determined by the \code{adj_p_filt input}) and suffiently (determined by \code{i_filt} and \code{med_i_filt}) imprinted,
#' the results of the differential imprinting analysis (which equals loss of imprinting if there's less imprinting in cases). This contains the same
#' columns as \code{impr_res}, but supplemented with:}
#'    \itemize{
#'      \item{DI_pval}{The p-value testing for differential imprinting using logistic regression, see the documentation for \code{LOItest_logreg} for more details.}
#'    }

impr_LOI_analysis <- function(DataList, SE, inbr, MinMinorAllelefreq=0.15, sym_filt=0.05, adj_p_filt=0.05, gof_filt=0.8, med_i_filt=0.8, i_filt=0.6){
  
  ImprData <- DataList[[1]]
  caseList <- DataList[[2]]

  for(LOC in names(ImprData)){
    if (ImprData[[LOC]]$allelefreq_prel[1] <= MinMinorAllelefreq || 
        ImprData[[LOC]]$allelefreq_prel[1] >= (1 - MinMinorAllelefreq)) {
      ImprData[[LOC]] <- NULL
    } else {
      ImprData[[LOC]]$sym <- MAGE::symmetry_gof(ImprData[[LOC]]$ref_count, 
                                                ImprData[[LOC]]$var_count, ImprData[[LOC]]$allelefreq_prel[1])
      if (ImprData[[LOC]]$sym[1] <= sym_filt) {
        ImprData[[LOC]] <- NULL
      }
    }
  }
  
  
  # Detect imprinted control loci
  impr_res <- data.frame()
  for(LOC in names(ImprData)){
    i_results <- MAGE::imprinting_est(ImprData[[LOC]]$ref_count, ImprData[[LOC]]$var_count, 
                                      allelefreq = ImprData[[LOC]]$allelefreq_prel[1], 
                                      SE = SE, inbr = inbr)
    # An additional robustified "median imprinting" across samples to be used as possible 
    # additional filter criterion:
    med_imp <- MAGE::median_imprinting(ImprData[[LOC]]$ref_count, ImprData[[LOC]]$var_count, 
                                       allelefreq = ImprData[[LOC]]$allelefreq_prel[1], inbr = inbr)
    # Write various results to a dataframe:
    results_z <- data.frame("position" = ImprData[[LOC]]$locus_id[1], "LRT" = i_results$LRT, 
                            "p" = i_results$p_value, "estimated.i" = i_results$est_i, "allele.frequency" = 
                              ImprData[[LOC]]$allelefreq_prel[1], "reference" = ImprData[[LOC]]$ref[1], "variant" = 
                              ImprData[[LOC]]$var[1], "med_cov" = ImprData[[LOC]]$coverage[1], "nr_samples" = 
                              nrow(ImprData[[LOC]]), "GOF" = i_results$GOF_likelihood, "symmetry" = 
                              ImprData[[LOC]]$sym[1], "med_impr" = med_imp, stringsAsFactors = FALSE)
    impr_res <- rbind(impr_res, results_z)
  }
  
  
  # Retain significantly imprinted loci (5% FDR) utilizing some additional filters, amongst
  # which a custom  Goodness-Of-Fit which more or less corresponds to a locus' likelihood of
  # the imprinted model*coverage; 0.8 is a good cutoff. Other filter criteria are imprinting 
  # (0.6) and median imprinting (0.8)
  impr_res_FIN <- MAGE::final_filter(data_hash=NULL, impr_res, results_wd=NULL, gof_filt = gof_filt, 
                                     med_impr_filt = med_i_filt, i_filt = i_filt, adj_p_filt = adj_p_filt, file_all = FALSE, file_impr = FALSE, 
                                     file_all_counts = FALSE, file_impr_counts = FALSE)
  
  
  # From among actually imprinted loci, detect differential expression in case data
  pos_impr <- as.character(impr_res_FIN$position)
  LOI_res <- impr_res_FIN
  LOI_res$DI_pval <- 1
  for(LOC in pos_impr){
    CData <- ImprData[[LOC]]; TData <- caseList[[LOC]]
    p_DI <- MAGE::LOItest_logreg(CData$ref_count, CData$var_count, 
                                 TData$ref_count, TData$var_count)$p.value
    LOI_res$DI_pval[LOI_res$position == LOC] <- p_DI
  }
  
  return(list("impr_res"=impr_res, "LOI_res"=LOI_res))
  
}
