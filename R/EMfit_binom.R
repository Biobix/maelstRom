#' Models allelic RNAseq counts using an expectation-maximization fit of a binomial mixture distribution 
#'
#' \code{EMfit_binom} estimates, per locus, parameters of an assumed binomial mixture model via expectation maximization (genotype
#' frequencies, allelic bias i.e. heterozygous p-parameter).
#' Through this fit, per-samples genotype probabilities are obtained as well. A fit assuming no allelic bias (heterozygous p-parameter = 0.5) is
#' performed as well, for the purpose of significant allelic bias detection via likelihood ratio test.
#'
#' @param data_counts Data frame. Data frame of a SNP with reference and variant counts ("ref_count" and
#'     "var_count", respectively) for each sample ("sample").
#' @param SE Number. Sequencing error rate.
#' @param allelefreq Number. Allele frequency. Only used when pInitEst is TRUE (default = 0.5)
#' @param inbr Number. Degree of inbreeding (default = 0).
#' @param dltaco Number. Minimal difference between 2 iterations (default = 0.001).
#' @param HWE Logical. Should HWE be used for allele frequency estimation, not recommended (default = FALSE).
#' @param p_InitEst Logical. Calculate initial estimates of pr, pv and prv from allelefreq, not recommended (default = FALSE).
#' @export
#' @return A list containing the following components:
#' \item{AB}{The estimated allelic shift.}
#' \item{AB_lrt}{The test statistic of the likelihood ratio test.}
#' \item{AB_p}{The p-value of the likelihood ratio test.}
#' \item{GOF}{The goodness-of-fit value based on the corrected likelihood.}
#' \item{nrep}{The number of iterations.}
#' \item{quality}{Indicates the quality of a locus. An "!" indicates bad data or a bad fit due to no apparent heterozygosity.}
#' \item{rho_vv}{The variant (homozygote) genotype probability for the SNP.}
#' \item{rho_rr}{The reference (homozygote) genotype probability for the SNP.}
#' \item{rho_rv}{The heterozyous genotype probability for the SNP.}
#' \item{rho_vv_H0}{The variant (homozygote) genotype probability with Allelic Bias = 0.5.}
#' \item{rho_rr_H0}{The reference (homozygote) genotype probability with Allelic Bias = 0.5.}
#' \item{rho_rv_H0}{The heterozyous genotype probability with Allelic Bias = 0.5.}
#' \item{data_hash}{Data frame. Input data frame with extra columns: allelefreq, genotype and genotype probabilities (prr, prv, pvv) per sample.}

EMfit_binom <- function(data_counts, SE, allelefreq = 0.5, inbr = 0, dltaco = 0.001, HWE = FALSE, p_InitEst = FALSE) {

  probshift <- 0.5 # Initial estimate position heterozygous peak
  dlta <- 1 # Track change in probshift between iterations
  nrep <- 0 # Track number of iterations
  if (p_InitEst) {
    prv <- 2 * allelefreq * (1 - allelefreq) * (1 - inbr) # Corrected for inbreeding
    pr <- allelefreq ^ 2 + inbr * allelefreq * (1 - allelefreq)
    pv <- (1 - allelefreq) ^ 2 + inbr * allelefreq * (1 - allelefreq)
  } else{
    pr <- 1 / 3
    prv <- 1 / 3
    pv <- 1 / 3
  }

  pr_H0 <- pr # Will be used later when constructing the LRT testing for detection of significant heterozygous peak shift
  pv_H0 <- pv
  prv_H0 <- prv

  quality <- ""

  while (dlta > dltaco & nrep < 100) {
    nrep <- nrep + 1
    probshiftold <- probshift

    spv <- pv * dbinom(data_counts$var_count, data_counts$ref_count + data_counts$var_count, prob = 1 - SE)
    sprv <- prv * dbinom(data_counts$ref_count, data_counts$ref_count + data_counts$var_count, prob = probshift)
    spr <- pr * dbinom(data_counts$ref_count, data_counts$ref_count + data_counts$var_count, prob = 1 - SE)

    pdata <- rowSums(cbind(spv, sprv, spr))

    if (any(pdata==0)) {
      # This usually only occurs at very high-count data
      # If this happens, the log-binomial densities are calculated. The highest one is chosen (rr, rv, vv)
      # and set to 1 (the other two are set to zero).
      # This may seem like a rough approximation, but is close enough to the truth for high-count data,
      # for which the binomial probabilities are expected to have a several order-of-magnitudes difference
      ProblemCases <- which(pdata == 0)
      for (case in ProblemCases) {
        var_part <- data_counts$var_count[case]
        ref_part <- data_counts$ref_count[case]
        spv_part <- dbinom(var_part, ref_part + var_part, 1 - SE, log = TRUE)
        sprv_part <- dbinom(ref_part, ref_part + var_part, probshift, log = TRUE)
        spr_part <- dbinom(ref_part, ref_part + var_part, 1 - SE, log = TRUE)
        spvec <- c(spr_part, sprv_part, spv_part)
        if (spr_part == max(spvec)) {
          spr[case] <- 1
          pdata[case] <- 1
        } else if (sprv_part == max(spvec)) {
          sprv[case] <- 1
          pdata[case] <- 1
        } else {
          spv[case] <- 1
          pdata[case] <- 1
        }
      }
    }

    # Normalise probabilities as to sum to 1
    spv <- spv / pdata
    sprv <- sprv / pdata
    spr <- spr / pdata

    if (sum(sprv) == 0) {
      # No evidence for heterozygosity in the data, the EM-model cannot be fitted
      # The quality-flag is assigned ("!") and functions skips to the output-phase
      quality <- "!"
      break
    }

    # start EM-algorithms

    # E-step: z-terms are equivalent to spr, sprv and spv

    # M-step:
    # Calculate each peak's (rr, rv, vv) relative weight in the mixture distribution,
    # either imposing HWE or not
    allelefreq <- mean(spr) + mean(sprv) / 2
    if (HWE) { # Calculate weights while imposing HWE (not recommended)
      prv <- 2 * allelefreq * (1 - allelefreq) * (1 - inbr)
      pr <- allelefreq ^ 2 + inbr * allelefreq * (1 - allelefreq)
      pv <- (1 - allelefreq) ^ 2 + inbr * allelefreq * (1 - allelefreq)
    }
    else { # Calculate weights through regular EM (recommended)
      pv <- mean(spv)
      prv <- mean(sprv)
      pr <- mean(spr)
    }

    # MLE of binomial p-parameter of heterozygous peak
    probshift <- sum(sprv * (data_counts$ref_count / (data_counts$ref_count + data_counts$var_count))) / sum(sprv)
    # Small correction so this shift does not get more extreme than the sequencing error (p-parameter for the homozygous peaks):
    probshift[probshift > (1 - 5 * SE)] <- (1 - 5 * SE)
    probshift[probshift < 5 * SE] <- 5 * SE
    # Convergence criterion:
    dlta <- abs(probshiftold - probshift)

    allelefreq <- mean(spr) + mean(sprv) / 2
    data_counts$allelefreq <- allelefreq
  }

  data_counts$prr <- spr; data_counts$prv <- sprv; data_counts$pvv <- spv
  data_counts$genotypeN <- c("rr", "rv", "vv")[apply(data_counts[, c("prr", "prv", "pvv")], 1, function(x) which(x == max(x))[1])]

  # A GOF-heuristic
  dmixase_corrected <- MAGE::pmf_binomMix(data_counts$ref_count, data_counts$var_count, probshift, SE, mean(spr), mean(spv), mean(sprv)) * (data_counts$ref_count + data_counts$var_count + 1)
  logLikelihood <- mean(log(dmixase_corrected))
  #probshift_adj <- -(as.numeric(probshift) - 0.5) / 0.5


  # Fit the model without shift to perform an LRT testing for a significant heterozygous peak shift away from 0.5
  dlta <- 1
  nrep2 <- 0
  MyFlag <- ""
  Q <- 1000 # Q-criterion is used here to check for convergence (as probshift is now set)

  while (dlta > dltaco & nrep2 < 100) {
    Qold <- Q
    nrep2 <- nrep2 + 1
    spv_H0 <- pv_H0 * dbinom(data_counts$var_count, data_counts$ref_count + data_counts$var_count, prob = 1 - SE)
    sprv_H0 <- prv_H0 * dbinom(data_counts$ref_count, data_counts$ref_count + data_counts$var_count, prob = 0.5) # probshift set to 0.5
    spr_H0 <- pr_H0 * dbinom(data_counts$ref_count, data_counts$ref_count + data_counts$var_count, prob = 1 - SE)

    pdata <- rowSums(cbind(spr_H0, spv_H0, sprv_H0))

    if (any(pdata == 0)) {
      ProblemCases <- which(pdata == 0)
      for (case in ProblemCases) {
        var_part <- data_counts$var_count[case]
        ref_part <- data_counts$ref_count[case]
        spv_part <- dbinom(var_part, ref_part + var_part, 1 - SE, log = TRUE)
        sprv_part <- dbinom(ref_part, ref_part + var_part, 0.5, log = TRUE)
        spr_part <- dbinom(ref_part, ref_part + var_part, 1 - SE, log = TRUE)
        spvec <- c(spr_part, sprv_part, spv_part)
        if (spr_part == max(spvec)) {
          spr_H0[case] <- 1
          pdata[case] <- 1
        } else if (sprv_part == max(spvec)) {
          sprv_H0[case] <- 1
          pdata[case] <- 1
        } else {
          spv_H0[case] <- 1
          pdata[case] <- 1
        }
      }
    }

    spv_H0 <- spv_H0 / pdata
    sprv_H0 <- sprv_H0 / pdata
    spr_H0 <- spr_H0 / pdata

    if (sum(sprv_H0) == 0) {
      # Set a flag, signifying the LRT can't proceed
      MyFlag <- "!"
      break
    }


    allelefreq_H0 <- mean(spr_H0) + mean(sprv_H0) / 2
    if (HWE) {
      prv_H0 <- 2 * allelefreq_H0 * (1 - allelefreq_H0) * (1 - inbr)
      pr_H0 <- allelefreq_H0 ^ 2 + inbr * allelefreq_H0 * (1 - allelefreq_H0)
      pv_H0 <- (1 - allelefreq_H0) ^ 2 + inbr * allelefreq_H0 * (1 - allelefreq_H0)
    } else {
      pv_H0 <- mean(spv_H0)
      prv_H0 <- mean(sprv_H0)
      pr_H0 <- mean(spr_H0)
    }

    Q <- sum(ifelse(pr_H0 > 0, spr_H0 * log(pr_H0), 0) + spr_H0 * dbinom(data_counts$ref_count, data_counts$ref_count + data_counts$var_count, prob = 1 - SE, log = TRUE) +
               ifelse(prv_H0 > 0, sprv_H0 * log(prv_H0), 0) + sprv_H0 * dbinom(data_counts$ref_count, data_counts$ref_count + data_counts$var_count, prob = 0.5, log = TRUE) +
               ifelse(pv_H0 > 0, spv_H0 * log(pv_H0), 0) + spv_H0 * dbinom(data_counts$var_count, data_counts$ref_count + data_counts$var_count, prob = 1 - SE, log = TRUE))
    dlta <- abs(Qold - Q)
  }

  if (MyFlag == "!") {
    lrtstat <- NA
    pval <- NA
  } else {
    dmixh0 <- MAGE::pmf_binomMix(data_counts$ref_count, data_counts$var_count, probshift = 0.5, SE,  mean(spr_H0), mean(spv_H0), mean(sprv_H0))
    dmixase <- MAGE::pmf_binomMix(data_counts$ref_count, data_counts$var_count, probshift = probshift, SE, mean(spr), mean(spv), mean(sprv))
    lrtstat <- -2 * (sum(log(dmixh0)) - sum(log(dmixase)))
    pval <- pchisq(lrtstat, df = 1, lower.tail = F)
  }

  results <- list(AB = probshift, AB_lrt = lrtstat, AB_p = pval, GOF = logLikelihood, nrep = nrep, quality = quality, rho_vv = pv ,rho_rr = pr, rho_rv = prv, rho_vv_H0 = pv_H0 ,rho_rr_H0 = pr_H0, rho_rv_H0 = prv_H0, data_hash = data_counts)
  return(results)
}




