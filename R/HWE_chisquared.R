#' Chi-squared test assessing the Hardy-Weinberg assumption.
#'
#' \code{HWE_chisquared} performs a per-SNP chi-squared test comparing observed genotype
#'    frequencies to those expected when assuming Hardy-Weinberg Equilibrium (HWE).
#'    A per-SNP corrected inbreeding coefficient is used if \code{Fmedian} is lower than the
#'    mathematically minimal possible value of -q/(1-q), with q the minor allele frequency (i.e. this minimum is used instead).
#'
#' @param data List or dataframe. A named list of per-SNP data frames (or one such dataframe) for a single chromosome, with names
#'    corresponding to the position on the chromosome. The list contains genotype probabilities
#'    ("prr", "prv" and "pvv") for each sample ("sample").
#' @param inbr Number. Estimated or assumed degree of inbreeding within the population
#'    under study (default = 0).
#' @param results Data frame. Data frame with SNP results in each row, containing
#'    at least the chromosomal position of the SNP (results$position). If NULL, one dataframe is expected as the \code{data} input argument.
#' @return The \code{results} input data frame appended with the following columns:
#' \item{HWEpval}{The p-value of the chi-squared test.}
#' \item{HWE_TestStat}{The test statistic of the chi-squared test.}
#' @export

HWE_chisquared <- function(data, inbr = 0, results = NULL){
  
  if(is.null(results)){
    prrOBS <- sum(data$prr)
    prvOBS <- sum(data$prv)
    pvvOBS <- sum(data$pvv)
    Rall <- 2 * prrOBS + prvOBS
    Vall <- 2 * pvvOBS + prvOBS
    TOTall <- Rall + Vall
    p <- Rall / TOTall
    q <- Vall / TOTall
    
    # Fail-safe for negative values of inbreeding:
    # For a given locus with multiple samples, the mathematically lowest possible inbreeding value
    # actually equals -q/(1-q), with q the minor allele frequency.
    # So, if inbr is lower than the -q/(1-q) of a specific locus,
    # it needs to be put at -q/(1-q) to avoid negative "rr/rv/vv" frequencies.
    pqmin <- min(c(p, q))
    F_alt <- - pqmin / (1 - pqmin)
    if (inbr < F_alt) {
      inbr <- F_alt
    }
    
    # Avoid precision errors that would return negative values through max(0, ...)
    prrEXP <- max(0, (p ^ 2 + inbr * p * q)) * nrow(data)
    prvEXP <- max(0, (2 * p * q * (1 - inbr))) * nrow(data)
    pvvEXP <- max(0, (q ^ 2 + inbr * p * q)) * nrow(data)
    
    PVAL <- chisq.test(rbind(c(prrOBS, prvOBS, pvvOBS), c(prrEXP, prvEXP, pvvEXP)))$p.value
    STAT <- chisq.test(rbind(c(prrOBS, prvOBS, pvvOBS), c(prrEXP, prvEXP, pvvEXP)))$statistic
    
    return(list("HWEpval" = PVAL, "HWE_TestStat" = STAT))
  }else{
    
    inbrC <- inbr

    counter <- 0
    
    results$HWEpval <- 0
    results$HWE_TestStat <- 0
    
    for (ii in 1:nrow(results)) {
      inbr <- inbrC
      counter <- counter + 1
      DataNow <- data[[names(data)[ii]]]
      prrOBS <- sum(DataNow$prr)
      prvOBS <- sum(DataNow$prv)
      pvvOBS <- sum(DataNow$pvv)
      Rall <- 2 * prrOBS + prvOBS
      Vall <- 2 * pvvOBS + prvOBS
      TOTall <- Rall + Vall
      p <- Rall / TOTall
      q <- Vall / TOTall
      
      # Fail-safe for negative values of inbreeding:
      # For a given locus with multiple samples, the mathematically lowest possible inbreeding value
      # actually equals -q/(1-q), with q the minor allele frequency.
      # So, if inbr is lower than the -q/(1-q) of a specific locus,
      # it needs to be put at -q/(1-q) to avoid negative "rr/rv/vv" frequencies.
      pqmin <- min(c(p, q))
      F_alt <- - pqmin / (1 - pqmin)
      if (inbr < F_alt) {
        inbr <- F_alt
      }
      
      # Avoid precision errors that would return negative values through max(0, ...)
      prrEXP <- max(0, (p ^ 2 + inbr * p * q)) * nrow(DataNow)
      prvEXP <- max(0, (2 * p * q * (1 - inbr))) * nrow(DataNow)
      pvvEXP <- max(0, (q ^ 2 + inbr * p * q)) * nrow(DataNow)
      
      PVAL <- chisq.test(rbind(c(prrOBS, prvOBS, pvvOBS), c(prrEXP, prvEXP, pvvEXP)))$p.value
      STAT <- chisq.test(rbind(c(prrOBS, prvOBS, pvvOBS), c(prrEXP, prvEXP, pvvEXP)))$statistic
      results$HWEpval[counter] <- PVAL
      results$HWE_TestStat[counter] <- STAT
    }
    
    return(results)
    
  }
  
}
