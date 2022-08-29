#' Estimate sequencing error rate and inbreeding using an Expectation-Maximisation algorithm
#'
#' \code{AllelicMeta_est} calculates sequencing error rates, inbreeding
#'     coefficients, allele frequencies, genotypes and genotype probabilities from \code{ref_counts} and \code{var_counts}
#'     for a specific locus. To this end it fits a (regular) binomial mixture model with no allelic bias in heterozygotes (both alleles equally likely).
#'     This somewhat simplistic model (no overdispersion parameter, no allelic bias) makes for an inferior but very fast fit,
#'     which is ideal to get rough estimates of the allelic population metaparameters (sequencing error rate and inbreeding coefficient) by e.g. taking their mean or median across a large amount of loci.
#'     The per-sample allele frequencies and genotype (probabilities) are less reliable and, as such, not this function's main use.
#'
#' @param ref_counts Numeric vector. Reference counts.
#' @param var_counts Numeric vector. Variant counts.
#' @param deltaF Number. Expectation-Maximisation threshold, minimal difference between
#'     two consecutive iterations (default is 1e-08).
#' @param maxIT Number. Maximum number of iterations of the Expectation-Maximisation algorithm (default is 100).
#' @param SE_prior Number. Initial estimate of the sequencing error rate (default is 0.002).
#' @param F_inbr_prior Number. Initial estimate of the inbreeding coefficient used for calculating initial genotype frequencies (default is NULL, in which case the initial genotype frequences all get set to 1/3).
#' @param HetProb Number. Allelic bias in heterozygotes (expected reference over total allele count in RNAseq data; default is 0.5
#' @export
#' @return A list containing the following components:
#' \item{allelefreq}{The estimated allele frequency.}
#' \item{SE}{The estimated sequencing error rate.}
#' \item{F_inbr}{The estimated inbreeding coefficient.}
#' \item{genotypes}{The most likely genotype (rr, rv, vv) of each sample.}
#' \item{genoprobs}{The genotype probabilties (p(rr), p(rv), p(vv)) for each sample.}
#' \item{nrep}{The number of iterations}
#' @examples
#' AllelicMeta_est(c(5, 8, 10, 3, 5, 6, 23), c(8, 8, 6, 4, 4, 10, 0))
#' AllelicMeta_est(c(5, 0, 0, 3, 5, 1, 23), c(1, 8, 6, 2, 0, 10, 0),
#'     SE_prior = 0.2, F_inbr_prior = 0.1)

AllelicMeta_est <- function(ref_counts, var_counts, deltaF = 10^-8, maxIT = 100, SE_prior = 0.002, F_inbr_prior = NULL, HetProb = 0.5) {

  delta <- 1
  nrep <- 0

  SE <- SE_prior
  if(is.null(F_inbr_prior)) {
    pr <- 1/3
    pv <- 1/3
    prv <- 1/3
  } else {
    tots <- ref_counts + var_counts
    rf <- mean(ref_counts/tots)
    vf <- 1-rf
    rf <- (rf - SE) / (1 - 2 * SE)
    vf <- (vf - SE) / (1 - 2 * SE)
    pr <- rf^2 + F_inbr_prior * rf * vf
    pv <- vf^2 + F_inbr_prior * rf * vf
    prv <- 1 - pr - pv
  }

  while(delta > deltaF & nrep < maxIT) {

    SEold <- SE
    prold <- pr
    pvold <- pv
    prvold <- prv

    nrep <- nrep + 1
    spr <- pr * dbinom(ref_counts, ref_counts + var_counts, prob = 1 - SE)
    spv <- pv * dbinom(var_counts,ref_counts + var_counts, prob = 1 - SE)
    sprv <- prv * dbinom(ref_counts, ref_counts + var_counts, prob = HetProb)

    pdata <- rowSums(cbind(spv, sprv, spr))

    if (any(pdata == 0)) {
      ProblemCases <- which(pdata == 0)
      for(case in ProblemCases){
        var_part <- var_counts[case]
        ref_part <- ref_counts[case]
        spv_part <- dbinom(var_part, ref_part + var_part, 1 - SE, log = TRUE)
        sprv_part <- dbinom(ref_part, ref_part + var_part, HetProb, log = TRUE)
        spr_part <- dbinom(ref_part, ref_part + var_part, 1 - SE, log = TRUE)
        spvec <- c(spr_part, sprv_part, spv_part)
        if(spr_part == max(spvec)) {
          spr[case] <- 1
          pdata[case] <- 1
        } else if(sprv_part == max(spvec)) {
          sprv[case] <- 1
          pdata[case] <- 1
        } else {
          spv[case] <- 1
          pdata[case] <- 1
        }
      }
    }

    if(sum(spr + spv) == 0) {
      allelefreq <- 1
      SE <- 1
      F_inbr <- 1
      genoprobs <- data.frame("p(rr)" = spr, "p(rv)" = sprv, "p(vv)" = spv)
      genotypes <- rep("NR", length(ref_counts))
      returnOBJ <- list(allelefreq, SE, F_inbr, genotypes, genoprobs)
      names(returnOBJ) <- c("allelefreq", "SE", "F_inbr", "genotypes", "genoprobs")
      return(returnOBJ)
    }

    # pdata toepassen;
    spr <- spr / pdata
    spv <- spv / pdata
    sprv <- sprv / pdata

    pv <- mean(spv)
    prv <- mean(sprv)
    pr <- mean(spr)

    tot_counts <- ref_counts + var_counts
    Nrr <- sum(tot_counts * spr)
    Nvv <- sum(tot_counts * spv)
    Xrr <- sum(var_counts * spr)
    Xvv <- sum(ref_counts*spv)
    SE <- (Xrr + Xvv) / (Nrr + Nvv)

    SEdif <- abs(SEold - SE)
    prdif <- abs(prold - pr)
    pvdif <- abs(pvold - pv)
    prvdif <- abs(prvold - prv)
    delta <- max(c(SEdif, prdif, pvdif, prvdif))

  }

  # allelefreq = fraction reference alleles
  allelefreq <- mean(spr) + mean(sprv)/2
  genoprobs <- data.frame("p(rr)" = spr, "p(rv)" = sprv, "p(vv)" = spv)
  genotypes <- c("rr", "rv", "vv")[apply(genoprobs, 1, function(x) which(x == max(x))[1])]

  F_inbr <- 1 - (sum(sprv) / (2 * allelefreq * (1 - allelefreq) * length(ref_counts)))

  returnOBJ <- list(allelefreq, SE, F_inbr, genotypes, genoprobs, nrep)
  names(returnOBJ) <- c("allelefreq", "SE", "F_inbr", "genotypes", "genoprobs", "nrep")
  return(returnOBJ)
}
