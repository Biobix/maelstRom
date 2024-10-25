#' Prior filtering of loci.
#'
#' \code{prior_filter} filters loci prior to analysis using maelstRom: 
#' \itemize{
#'   \item a filter throwing out samples showcasing relevant allele counts for alleles not present in the ref_alleles column (if checkref_filter == TRUE)
#'   \item a filter throwing out loci showcasing a low minor allele fraction (if prior_allelefreq_filter == TRUE, governed by min_PrioAlleleFreq)
#'   \item a filter for median coverage (governed by min_median_cov, only looking at ref_count and var_count)
#'   \item a filter for the number of samples (governed by min_nr_samples)
#' }
#' Also, samples having zero-counts in both "ref_count" and "var_count" are discarded, as they're useless in further analyses.
#'
#' @param data_pos Data frame. Data frame of a SNP position with columns:
#'    "ref_alleles", "A", "T", "C", "G", "ref", "var", "ref_count" and "var_count".
#' @param min_median_cov Number. Minimal median coverage necessary to retain locus (default is 5).
#' @param min_nr_samples Number. Minimal number of samples necessary to retain locus (default is 30).
#' @param checkref_filter Logical. Should samples with a high allele count of alleles not present in the ref_alleles column of data_pos be filtered out?
#' With "high", we mean either the highest allele count in that sample, OR the second highest if there's heuristic evidence that the sample is not homozygous for its most common allele.
#' The particular check involves, for each sample, extracting its first and second most common allele, then using the entire population data in \code{data_pos} to get
#' genotype probabilities for homozygocity in the most common allele and heterozygosity in these alleles, using these allele's observed population frequencies and assuming
#' Hardy-Weinberg Equilibrium (inbreeding parameter = 0). These are then multiplied with that particular sample's corresponding evidence for homozygocity of the most common allele
#' (using a multinomial distribution assuming equal observed frequencies for all other alleles, which should be low and are dictated by e.g. sequencing errors)
#' and evidence for heterozygosity if its most common alleles (using a multinomial distribution assuming equal expression of these two alleles which should be high,
#' and equal expression of the two remaining nucleotides which should be low, and dictated by e.g. sequencing error rate).
#' @param prior_allelefreq_filter Logical. Should loci be filtered using a minimal prior minor allele frequency,
#' which is simply determined as the percent occurrence of the less common (variant) allele across nucleotides. (default is FALSE).
#' @param min_PriorAlleleFreq Number. Prior allele frequency to filter SNPs on if \code{prior_allelefreq_filter} is TRUE (default is 0.1).
#' @export
#' @return The data as data frame with total and median coverage or NULL if the SNP positions was filtered.

prior_filter <- function(data_pos, min_median_cov = 5, min_nr_samples = 30, checkref_filter = FALSE,
                         prior_allelefreq_filter = FALSE, min_PriorAlleleFreq = 0.1) {
  
  if(checkref_filter){
    ref_alleles <- unlist(strsplit(data_pos$ref_alleles[1], "/"))
  }
  
  if(checkref_filter | prior_allelefreq_filter){
    allelecount_samples <- lapply(1:nrow(data_pos), function(y) c(data_pos$A[y], data_pos$T[y], data_pos$C[y], data_pos$G[y]))
    allelecount_samples_percent <- lapply(1:nrow(data_pos), function(y) allelecount_samples[[y]] / sum(allelecount_samples[[y]]) * 100)
    allelecount_samples_table <- lapply(1:nrow(data_pos), function(y) setNames(allelecount_samples_percent[[y]], c("A", "T", "C", "G")))
    
    allelecount_loci <- colSums(matrix(unlist(allelecount_samples), ncol = 4, byrow = T))
    allelecount_loci_percent <- allelecount_loci / sum(allelecount_loci) * 100
    allelecount_loci_table <- setNames(allelecount_loci_percent, c("A", "T", "C", "G"))
  }
  
  #FILTERING
  if(checkref_filter){
    alleles_loci_ordered <- sort(allelecount_loci_table, TRUE)
    #whichsample <- hash::hash(data_pos$sample, 1:nrow(data_pos))
    marker <- rep(TRUE, nrow(data_pos))
    for (s in 1:nrow(data_pos)) {
      alleles_samples_ordered <- sort(allelecount_samples_table[[s]], TRUE)
      highest_sample_alleles <- names(alleles_samples_ordered)
      if (highest_sample_alleles[1] %in% ref_alleles && highest_sample_alleles[2] %in% ref_alleles) {
        #data_pos <- data_pos
        next
      } else if (!(highest_sample_alleles[1] %in% ref_alleles)) {
        #data_pos <- data_pos[-which(data_pos$sample == s), ]
        marker[s] <- FALSE
      } else if (!(highest_sample_alleles[2] %in% ref_alleles)) {
        p1 <- (as.numeric(alleles_samples_ordered[1]) + as.numeric(alleles_samples_ordered[2])) / 2
        p3 <- (as.numeric(alleles_samples_ordered[3]) + as.numeric(alleles_samples_ordered[4])) / 2
        p_data_hetero <- dmultinom(as.numeric(alleles_samples_ordered), prob = c(p1, p1, p3, p3))
        p1 <- as.numeric(alleles_samples_ordered[1])
        p2 <- (as.numeric(alleles_samples_ordered[2]) + as.numeric(alleles_samples_ordered[3]) + as.numeric(alleles_samples_ordered[4])) / 3
        p_data_homo <- dmultinom(as.numeric(alleles_samples_ordered), prob = c(p1, p2, p2, p2))
        p_allele <- ref_alleles[which(ref_alleles %in% highest_sample_alleles[1])]
        q_allele <- highest_sample_alleles[2]
        p <- as.numeric(allelecount_loci_table[which(names(allelecount_loci_table) == p_allele)])
        q <- as.numeric(allelecount_loci_table[which(names(allelecount_loci_table) == q_allele)])
        p_homo <- p ^ 2
        p_hetero <- 2 * p * q
        if (p_data_hetero * p_hetero >= p_data_homo * p_homo) {
          #data_pos <- data_pos[-which(data_pos$sample == s), ]
          marker[s] <- FALSE
        } else {
          #data_pos <- data_pos
          next
        }
      }
    }
    data_pos <- data_pos[marker,]
    if (nrow(data_pos) == 0) return(NULL)
  }
  

  if (prior_allelefreq_filter &&  allelecount_loci_percent[which(is.element(c("A", "T", "C", "G"), data_pos$var[1]))] / 100 < min_PriorAlleleFreq) return(NULL)
  
  
  #DELETE NON_STANDARD ALLELES FROM SEQUENCES
  if (any(data_pos$ref_count + data_pos$var_count == 0)) data_pos <- data_pos[- which(data_pos$ref_count + data_pos$var_count == 0), ]
  if (nrow(data_pos) == 0) return(NULL)
  
  
  #FILTER ON COVERAGE AND ON NUMBER OF SAMPLES
  data_pos$total <- data_pos$ref_count + data_pos$var_count
  data_pos$coverage <- median(data_pos$total)
  if (median(data_pos$total) < min_median_cov || nrow(data_pos) < min_nr_samples) return(NULL)
  
  
  return(data_pos)
}
