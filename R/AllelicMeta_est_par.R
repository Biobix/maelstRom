#' multi-locus wrapper function of AllelicMeta_est
#'
#' \code{AllelicMeta_est_par} is a wrapper function of \code{AllelicMeta_est} allowing multiple loci to be given as input (as lists)
#' and focussing solely on obtaining a reliable sequencing error rate and inbreeding coefficient estimate (the population metaparameters),
#' this via extra inputs specifying filter criteria so as to only retain high-quality loci.
#'
#' @param DataList List of dataframes. Each dataframes should at least contain a column of reference- and variant-allele counts (named "ref_count" and "var_count" respectively).
#' @param deltaF Number. Expectation-Maximisation threshold, minimal difference between
#'     two consecutive iterations (default is 1e-08).
#' @param maxIT Number. Maximum number of iterations of the Expectation-Maximisation algorithm (default is 100).
#' @param SE_prior Number. Initial estimate of the sequencing error rate (default is 0.002).
#' @param F_inbr_prior Number. Initial estimate of the inbreeding coefficient used for calculating initial genotype frequencies (default is NULL, in which case the initial genotype frequences all get set to 1/3).
#' @param HetProb Number. Allelic bias in heterozygotes (expected reference over total allele count in RNAseq data; default is 0.5
#' @param pA_filt Number. Minimum allowed estimated minor allele frequency for a locus to be considered in the metaparameter estimates.
#' @param SE_filt Number. Maximum allowed estimated sequencing error rate for a locus to be considered in the metaparameter estimates.
#' @param NumSamp_filt Number. Minimum allowed number of samples for a locus to be considered in the metaparameter estimates.
#' @param MedianCov_filt Number. Minimum allowed median coverage across samples (reference plus variant allele count) for a locus to be considered in the metaparameter estimates.
#' @export
#' @return A list containing the following components:
#' \item{DataList_out}{An updated DataList, now containing each locus' sequencing error estimate (est_SE), inbreeding coefficient estimate (est_inbr), reference allele frequency estimate (allelefreq_prel), preliminary genotype probabilities (prr_prel, prv_prel, pvv_prel) and preliminary genotype (genotype_prel)}
#' \item{SE_vec}{Vector of estimated sequencing error rates of reliable loci.}
#' \item{F_vec}{Vector of estimated inbreeding coefficients of reliable loci.}

AllelicMeta_est_par <- function(DataList, deltaF = 10^-8, maxIT = 100, SE_prior = 0.002, F_inbr_prior = NULL, HetProb = 0.5,
                                pA_filt = 0.15, SE_filt = 0.035, NumSamp_filt = 20, MedianCov_filt = 4) {
  
  SE_vec <- c()
  F_vec <- c()
  
  positions <- names(DataList)
  results <- data.frame()
  
  for(n in names(DataList)){
    MetaEst_res <- MAGE::AllelicMeta_est(ref_counts = DataList[[n]]$ref_count, var_counts = DataList[[n]]$var_count,
      deltaF = deltaF, maxIT = maxIT, SE_prior = SE_prior, F_inbr_prior = F_inbr_prior, HetProb = HetProb)
    # You can store each locus' Sequencing Error (SE) and inbreeding coefficient (F) estimate 
    # in its dataframe, if you want:
    DataList[[n]]$est_SE <- MetaEst_res$SE
    DataList[[n]]$est_inbr <- MetaEst_res$F_inbr
    DataList[[n]]$allelefreq_prel <- MetaEst_res$allelefreq
    DataList[[n]]$prr_prel <- (MetaEst_res$genoprobs)$p.rr.
    DataList[[n]]$prv_prel <- (MetaEst_res$genoprobs)$p.rv.
    DataList[[n]]$pvv_prel <- (MetaEst_res$genoprobs)$p.vv.
    DataList[[n]]$genotype_prel <- MetaEst_res$genotypes
    # Only take a locus' estimates into account if the locus is high-quality:
    if (!(MetaEst_res$allelefreq <= pA_filt || MetaEst_res$allelefreq >= (1 - pA_filt) 
          # allelefreq returns allele frequency of the ref-allele,
          # so just in case this is the minor allele on population level 
          # (even though it is the most expresse one across all RNAseq data),
          # we perform the filtering like this
          || MetaEst_res$SE > SE_filt) & nrow(DataList[[n]])>=NumSamp_filt & 
        median(DataList[[n]]$ref_count + DataList[[n]]$var_count) >= MedianCov_filt) {
      SE_vec <- c(SE_vec, MetaEst_res$SE)
      F_vec <- c(F_vec, MetaEst_res$F_inbr)
    }
  }
  
  returnOBJ <- list("DataList_out" = DataList, "SE_vec" = SE_vec, "F_vec" = F_vec)
  return(returnOBJ)
}
