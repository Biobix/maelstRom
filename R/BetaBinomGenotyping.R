#' multi-locus wrapper function of EMfit_betabinom_robust for the purpose of genotyping
#'
#' \code{BetaBinomGenotyping} is a wrapper function of \code{EMfit_betabinom_robust} allowing multiple loci to be given as input (as lists).
#' Besides calling the latter to perform a EM-fit and update the inputted per-locus dataframes as listed in \code{EMfit_betabinom_robust} output,
#' it also returns a results-dataframe with one row for every locus, containing genotyping- and Allelic Bias-detection results.
#' This function is mainly useful when following MAGE's vignette step-by-step, as its input requirement are rather stictly dependent on performing all previous steps in that pipeline.
#' See the help page of \code{EMfit_betabinom_robust} for more information on these outputs
#' 
#' @param DataList List of dataframes. Each dataframes should at least contain the columns "ref", ref_count", "var", "var_count" and "est_SE",
#' which should be the case if following the vignette up to this function's appearance.
#' @param allelefreq,SE,inbr,dltaco,HWE,p_InitEst,ThetaInits,ReEstThetas,NoSplitHom,NoSplitHet,ResetThetaMin,ResetThetaMax,DistRob,CookMargin,LikEmpNum,LikMargin,NumHetMin,MaxOutFrac,thetaTRY,fitH0 All remaining parameters of \code{EMfit_betabinom_robust}
#' @export
#' @return A list containing the following components:
#' \item{DataList_out}{The updated DataList, corresponding to the \code{data_hash} output of \code{EMfit_betabinom_robust}}
#' \item{Geno_AB_res}{A dataframe containing per-locus results of \code{EMfit_betabinom_robust} EM-fit as well as some other metrics, namely:}
#'    \itemize{
#'      \item{position}{The locus' name, according to \code{names(DataList)}}
#'      \item{probshift}{Fitted reference allele fraction in RNAseq reads, indicating allelic bias when different from 0.5}
#'      \item{LRT}{The likelihood ratio test statistic, testing for significant allelic bias}
#'      \item{p}{The likelihood ratio test p-value, testing for significant allelic bias}
#'      \item{quality}{Equals "!" if the sample contains no fitted heterozygotes, otherwise ""}
#'      \item{allele_frequency}{estimated reference allele frequency in the population}
#'      \item{reference}{reference allele nucleotide}
#'      \item{variant}{variant allele nucleotide}
#'      \item{est_SE}{per-locus sequencing error rate estimate, as outputted earlier by \code{AllelicMeta_est}}
#'      \item{coverage}{median coverage across samples}
#'      \item{nr_samples}{number of samples covering this locus with at least one reference- or variant-count}
#'      \item{median_AB}{median allelic bias, as outputted by \code{median_AB}}
#'      \item{rho_rr}{fitted reference homozygous fraction in the population}
#'      \item{rho_rv}{fitted heterozygous fraction in the population}
#'      \item{rho_vv}{fitted variant homozygous fraction in the population}
#'      \item{theta_hom}{fitted overdispersion parameter for the homozygous PMFs}
#'      \item{theta_het}{fitted overdispersion parameter for the heterozygous PMF}
#'      \item{theta_hom_NoShift}{fitted overdispersion parameter for the homozygous PMFs assuming no allelic bias}
#'      \item{theta_het_NoShift}{fitted overdispersion parameter for the heterozygous PMF assuming no allelic bias}
#'      \item{Chi2PVAL}{p-value of a chi square test assessing Hardy-Weinberg-Equilibrium on the locus given the inbreeding coefficient metaparameter; see \code{HWE_chisquared}}
#'      \item{Chi2STAT}{test statistic of a chi square test assessing Hardy-Weinberg-Equilibrium on the locus given the inbreeding coefficient metaparameter; see \code{HWE_chisquared}}
#'    }

BetaBinomGenotyping <- function(DataList, allelefreq=0.5, SE, inbr = 0, dltaco = 10^-6, HWE = FALSE, p_InitEst = FALSE, 
                                ThetaInits = "moment", ReEstThetas = "moment", NoSplitHom = TRUE, NoSplitHet = TRUE,
                                ResetThetaMin = 10^-10, ResetThetaMax = 10^-1, DistRob = "Cook", CookMargin = 5, LikEmpNum = 1000, LikMargin = 0,
                                NumHetMin = 5, MaxOutFrac = 0.5, thetaTRY = c(10^-1, 10^-3, 10^-7), fitH0 = TRUE){
  positions <- names(DataList)
  results <- data.frame()
  for (z in positions) {
    MAGEres <- MAGE::EMfit_betabinom_robust(data_counts = DataList[[z]], 
                                            SE = SE, inbr = inbr,
                                            allelefreq=allelefreq, dltaco = dltaco, HWE = HWE, p_InitEst = p_InitEst, 
                                            ThetaInits = ThetaInits, ReEstThetas = ReEstThetas, NoSplitHom = NoSplitHom, NoSplitHet = NoSplitHet,
                                            ResetThetaMin = ResetThetaMin, ResetThetaMax = ResetThetaMax, DistRob = DistRob, CookMargin = CookMargin, LikEmpNum = LikEmpNum, LikMargin = LikMargin,
                                            NumHetMin = NumHetMin, MaxOutFrac = MaxOutFrac, thetaTRY = thetaTRY, fitH0 = fitH0)
    DataList[[z]] <- MAGEres$data_hash
    # median_AB also calculates a robust median AB,
    # besides the AB as determined during EMfit_betabinom_robust's fitting procedure.
    # This can be used as additional filter when detecting significant AB:
    med_AB <- MAGE::median_AB(DataList[[z]]$ref_count, DataList[[z]]$var_count, 
                              DataList[[z]]$allelefreq[1], inbr) 
    res_loc <- data.frame("position" = z, "probshift" = as.numeric(MAGEres$AB), 
                          "LRT" = as.numeric(MAGEres$AB_lrt), "p" = as.numeric(MAGEres$AB_p),
                          "quality" = MAGEres$quality, "allele.frequency" = DataList[[z]]$allelefreq[1], 
                          "reference" = DataList[[z]]$ref[1], "variant" = DataList[[z]]$var[1], 
                          "est_SE" = DataList[[z]]$est_SE[1], "coverage" =  median(DataList[[z]]$ref_count+DataList[[z]]$var_count),
                          "nr_samples" = nrow(DataList[[z]]), "median_AB" = med_AB, 
                          "rho_rr" = MAGEres$rho_rr, "rho_rv" = MAGEres$rho_rv, "rho_vv" = MAGEres$rho_vv, 
                          "theta_hom" = MAGEres$theta_hom, "theta_het" = MAGEres$theta_het,
                          "theta_hom_NoShift" = MAGEres$theta_hom_NoShift, 
                          "theta_het_NoShift" = MAGEres$theta_het_NoShift, stringsAsFactors = FALSE)
    results <- rbind(results, res_loc) # results; one position per line
  }
  results <- MAGE::HWE_chisquared(data = DataList, inbr, results = results)
  results$Chi2PVAL[is.na(results$Chi2PVAL)] <- 
    results$Chi2STAT[is.na(results$Chi2STAT)] <- -1
  return(list("DataList_out" = DataList,"Geno_AB_res" = results))
}
