#' Performs a differential Allelic Divergence analysis on a control- and case-dataset
#'
#' \code{dAD_analysis} is a wrapper function bundling \code{EMfit_betabinom_robust}, \code{EMfit_betabinom}, \code{EMfit_betabinom_popcomb},
#' \code{pmf_betabinomMix}, and \code{HWE_chisquared} to run the entire differential Allelic Divergence detection analysis as discussed in the package vignette.
#' To this end, it goes over the following steps:
#' \enumerate{
#'   \item Calling \code{EMfit_betabinom_robust} on both the control- and case-dataset separately for the sole purpose of outlier detection;
#'   detected outliers in both datasets play no part in future Expectation-Maximization-fits.
#'   \item Perform one fit on the joint (case+control) data using \code{EMfit_betabinom}, i.e. assuming all population parameters are shared including
#'   the heterozygous overdispersion parameter, i.e. assuming no differential Allelic Divergence.
#'   \item Perform a fit on the control- and case-data using \code{EMfit_betabinom_popcomb}, which assumes all population parameters are shared
#'   except for the heteroyzous overdispersion parameter, i.e. a fit accomodating differential Allelic Divergence
#'   \item Statistically test for differential Allelic Divergence using a likelihood ratio test with one degree of freedom comparing the
#'   fits from the two previous points.
#'   \item Filling out the results dataframe \code{dAD_res}.
#' }
#' 
#' @param Datalist A list containing two lists of dataframes, which in turn contain control- and case-data (in that order). Subsequent entries in both lists
#' correspond to subsequent loci (lists must be the same size and named per-locus). Each of the dataframes should at least contain both a \code{ref_count} and
#' \code{var_count} column.
#' @param SE Number. Sequencing error rate population metaparameter.
#' @param inbr Number. Inbreeding coefficient population metaparameter.
#' @param allelefreq,dltaco,HWE,p_InitEst,ThetaInits,ReEstThetas,NoSplitHom,NoSplitHet,ResetThetaMin,ResetThetaMax,DistRob,CookMargin,LikEmpNum,LikMargin,NumHetMin,MaxOutFrac,thetaTRY Remaining
#' control parameters to be passed to \code{EMfit_betabinom_robust}, \code{EMfit_betabinom}, and \code{EMfit_betabinom_popcomb}. See their respective documentation
#' for more information. Usually, these parameter's default values are fine.
#' @export
#' @return A dataframe with one row per locus, containing for that locus in its columns:
#' \item{LocName}{Locus name.}
#' \item{PiFitH0}{Heterozygous pi parameter (reference allele fraction) for the joint fit sharing all parameters.}
#' \item{PiFitH1}{Heterozygous pi parameter (reference allele fraction) for the fit allowing differential Allelic Divergence.}
#' \item{ThetaHetH0}{Heterozygous theta parameter (overdispersion) for the joint fit sharing all parameters.}
#' \item{ThetaHetCTRL}{Heterozygous theta parameter (overdispersion) of controls for the fit allowing differential Allelic Divergence.}
#' \item{ThetaHetCASE}{Heterozygous theta parameter (overdispersion) of cases for the fit allowing differential Allelic Divergence.}
#' \item{RhoHetH0}{Heterozygous rho parameter (overdispersion, ranging from 0 to 1) for the joint fit sharing all parameters.}
#' \item{RhoHetCTRL}{Heterozygous rho parameter (overdispersion, ranging from 0 to 1) of controls for the fit allowing differential Allelic Divergence.}
#' \item{RhoHetCASE}{Heterozygous rho parameter (overdispersion, ranging from 0 to 1) of cases for the fit allowing differential Allelic Divergence.}
#' \item{NumHetCTRL}{Estimate number of heterozygotes in controls, according to its separate outlier-detection fit from step 1 above. Outliers are themselves genotyped and counted towards this number, just not used in the actual fitting procedure.}
#' \item{NumHetCASE}{Estimate number of heterozygotes in cases, according to its separate outlier-detection fit from step 1 above. Outliers are themselves genotyped and counted towards this number, just not used in the actual fitting procedure.}
#' \item{RobFlagCTRL}{The \code{RobFlag} output of the separate outlier-detection fit from step 1 on controls, see \code{EMfit_betabinom_robust}'s documentation.}
#' \item{RobFlagCASE}{The \code{RobFlag} output of the separate outlier-detection fit from step 1 on cases, see \code{EMfit_betabinom_robust}'s documentation.}
#' \item{HWECTRL}{HWE chi squared p-value for controls, see \code{HWE_chisquared}. Performed on genotyping data according to the outlier-detection fit of step 1.}
#' \item{HWECASE}{HWE chi squared p-value for cases, see \code{HWE_chisquared}. Performed on genotyping data according to the outlier-detection fit of step 1.}
#' \item{CovCTRL_mean}{Mean coverage (reference + variant count) across control data.}
#' \item{CovCASE_mean}{Mean coverage (reference + variant count) across case data.}
#' \item{CovCTRL_med}{Median coverage (reference + variant count) across control data.}
#' \item{CovCASE_med}{Median coverage (reference + variant count) across case data.}
#' \item{NumOutCTRL}{Number of outliers detected in control data.}
#' \item{NumOutCASE}{Number of outliers detected in case data.}
#' \item{QualityCTRL}{Quality flag returned by \code{EMfit_betabinom_robust} (see documentation) of the outlier-detection fit on control data (step 1 above).}
#' \item{QualityCASE}{Quality flag returned by \code{EMfit_betabinom_robust} (see documentation) of the outlier-detection fit on case data (step 1 above).}
#' \item{pr}{Reference homozygote fraction according to the differential Allelic Divergence accomodating fit (different heterozygous overdispersion parameters).}
#' \item{prv}{Heterozygote fraction according to the differential Allelic Divergence accomodating fit (different heterozygous overdispersion parameters).}
#' \item{pv}{Variant homozygote fraction according to the differential Allelic Divergence accomodating fit (different heterozygous overdispersion parameters).}
#' \item{ThetaHomH0}{Homozygous theta parameter (overdispersion) for the joint fit sharing all parameters.}
#' \item{ThetaHomH1}{Homozygous theta parameter (overdispersion) for the fit allowing differential Allelic Divergence.}

dAD_analysis <- function(DataList, SE, inbr = 0, allelefreq=0.5, dltaco = 10^-6, HWE = FALSE, p_InitEst = FALSE, 
                         ThetaInits = "moment", ReEstThetas = "moment", NoSplitHom = TRUE, NoSplitHet = TRUE,
                         ResetThetaMin = 10^-10, ResetThetaMax = 10^-1, DistRob = "Cook", CookMargin = 5, LikEmpNum = 1000, LikMargin = 0,
                         NumHetMin = 5, MaxOutFrac = 0.5, thetaTRY = c(10^-1, 10^-3, 10^-7)){
  
  controlListP <- DataList[[1]]
  caseListP <- DataList[[2]]
  positions <- names(controlListP)
  
  dAD_res <- data.frame(LocName = names(controlListP), PiFitH0 = 0, PiFitH1 = 0, ThetaHetH0 = 0,
                        ThetaHetCTRL = 0, ThetaHetCASE = 0, RhoHetH0 = 0, RhoHetCTRL = 0, RhoHetCASE = 0, LRTpval = 0, 
                        NumHetCTRL = 0, NumHetCASE = 0, RobFlagCTRL = "", RobFlagCASE = "", HWECTRL = 0,
                        HWECASE = 0, CovCTRL_mean = 0, CovCASE_mean = 0, CovCTRL_med = 0, CovCASE_med = 0,
                        NumOutCTRL = 0, NumOutCASE = 0, QualityCTRL = "N", QualityCASE = "N", pr = 0, 
                        prv = 0, pv = 0, ThetaHomH0 = 0, ThetaHomH1 = 0)
  
  for (LOC in positions) {
    
    CTRL_DF <- data.frame("ref_count" = controlListP[[LOC]]$ref_count, 
                          "var_count" = controlListP[[LOC]]$var_count, "isCase" = 0)
    CASE_DF <- data.frame("ref_count" = caseListP[[LOC]]$ref_count, 
                          "var_count" = caseListP[[LOC]]$var_count, "isCase" = 1)
    
    # Previously, we used EMfit_betabinom_robust() to ROBUSTLY fit our models by removing
    # outliers via Cook's distance. Now however, we're working with data from two different
    # sources (control- and case tissue) and have no way of knowing in advance how similar
    # these are. As such, outlier detections should happen ON BOTH SETS SEPARATELY yet the
    # two hypotheses we'll be fitting share all or some parameters between the two. 
    # For this, we can use the EMfit_betabinom_robust() function with its "fitH0" set to 
    # FALSE, which will not complete the entire eqtl-detection pipeline, but will cut it 
    # short before fitting the unshifted model. By running this function on both our 
    # datasets separately in advance of the dAD-relevant fits, we ensure correct outlier
    # detection ánd the use of the same dataset in our upcoming likelihood ratio tests.
    
    # 1. Detect and extract outliers using EMfit_betabinom_robust()
    OUTfitCTRL <- MAGE::EMfit_betabinom_robust(data_counts = CTRL_DF, SE = SE, inbr = inbr, 
      allelefreq=allelefreq, dltaco = dltaco, HWE = HWE, p_InitEst = p_InitEst, ThetaInits = ThetaInits, ReEstThetas = ReEstThetas, 
      NoSplitHom = NoSplitHom, NoSplitHet = NoSplitHet, ResetThetaMin = ResetThetaMin, ResetThetaMax = ResetThetaMax, DistRob = DistRob,
      CookMargin = CookMargin, LikEmpNum = LikEmpNum, LikMargin = LikMargin, NumHetMin = NumHetMin, MaxOutFrac = MaxOutFrac, 
      thetaTRY = thetaTRY, 
      fitH0 = FALSE)
    OUTfitCASE <- MAGE::EMfit_betabinom_robust(data_counts = CASE_DF, SE = SE, inbr = inbr,
      allelefreq=allelefreq, dltaco = dltaco, HWE = HWE, p_InitEst = p_InitEst, ThetaInits = ThetaInits, ReEstThetas = ReEstThetas, 
      NoSplitHom = NoSplitHom, NoSplitHet = NoSplitHet, ResetThetaMin = ResetThetaMin, ResetThetaMax = ResetThetaMax, DistRob = DistRob,
      CookMargin = CookMargin, LikEmpNum = LikEmpNum, LikMargin = LikMargin, NumHetMin = NumHetMin, MaxOutFrac = MaxOutFrac, 
      thetaTRY = thetaTRY,
      fitH0 = FALSE)
    OUTfitCTRL_DH <- OUTfitCTRL$data_hash; OUTfitCASE_DH <- OUTfitCASE$data_hash
    CTRL_DF$Outlier <- OUTfitCTRL_DH$Outlier; CASE_DF$Outlier <- OUTfitCASE_DH$Outlier
    CurDF <- rbind(CTRL_DF, CASE_DF)
    # These results are also handy for an estimation of the number of heterozygotes in 
    # controls and tumors AND the number of outliers, both good filter criteria.
    # We can also include a "RobFlag" which gives more information about outlier detection 
    # (e.g. none detected, or so unreasonably many that none were removed)
    dAD_res$NumHetCTRL[dAD_res$LocName == LOC] <- sum(OUTfitCTRL_DH$prv)
    dAD_res$NumHetCASE[dAD_res$LocName == LOC] <- sum(OUTfitCASE_DH$prv)
    dAD_res$NumOutCTRL[dAD_res$LocName == LOC] <- sum(CTRL_DF$Outlier)
    dAD_res$NumOutCASE[dAD_res$LocName == LOC] <- sum(CASE_DF$Outlier)
    dAD_res$RobFlagCTRL[dAD_res$LocName == LOC] <- OUTfitCTRL$RobFlag
    dAD_res$RobFlagCASE[dAD_res$LocName == LOC] <- OUTfitCASE$RobFlag
    
    
    # 2. Perform one fit on the entire (non-outlying) data, i.e. all parameters shared,
    # i.e. assuming no dAD; the null hypothesis in dAD detection.
    # Remark this uses the non-robust fitting function, since outliers were already detected.
    NOdAD_fit <- MAGE::EMfit_betabinom(data_counts = CurDF[CurDF$Outlier == 0,], SE = SE, inbr = inbr, 
                                       allelefreq=allelefreq, dltaco = dltaco, HWE = HWE, p_InitEst = p_InitEst, 
                                       ThetaInits = ThetaInits, ReEstThetas = ReEstThetas, NoSplitHom = NoSplitHom, NoSplitHet = NoSplitHet,
                                       ResetThetaMin = ResetThetaMin, ResetThetaMax = ResetThetaMax, thetaTRY = thetaTRY,                                 
                                       fitH0 = FALSE)
    NOdAD_fit_DF <- NOdAD_fit$data_hash
    PiH0 <- NOdAD_fit$AB
    rho_rr <- NOdAD_fit$rho_rr; rho_rv <- NOdAD_fit$rho_rv; rho_vv <- NOdAD_fit$rho_vv
    ThetaHomH0 <- NOdAD_fit$theta_hom; ThetaHetH0 <- NOdAD_fit$theta_het
    
    dAD_res$PiFitH0[dAD_res$LocName == LOC] <- PiH0
    dAD_res$ThetaHetH0[dAD_res$LocName == LOC] <- ThetaHetH0
    
    
    # 3. Perform a fit on the (non-outlying) data allowing separate theta_het parameters 
    # for control- and case-data
    FullFit <- MAGE::EMfit_betabinom_popcomb(data_counts = CurDF[CurDF$Outlier == 0,], SE = SE, inbr = inbr, 
                                             allelefreq=allelefreq, dltaco = dltaco, HWE = HWE, p_InitEst = p_InitEst, 
                                             ThetaInits = ThetaInits, ReEstThetas = ReEstThetas, NoSplitHom = NoSplitHom, NoSplitHet = NoSplitHet,
                                             ResetThetaMin = ResetThetaMin, ResetThetaMax = ResetThetaMax, thetaTRY = thetaTRY,
                                             probshift_init = PiH0)
    ParamVec <- FullFit$ParamVec
    dAD_res$PiFitH1[dAD_res$LocName == LOC] <- ParamVec["probshift"]
    dAD_res$ThetaHetCTRL[dAD_res$LocName == LOC] <- ParamVec["theta_het_control"]
    dAD_res$ThetaHetCASE[dAD_res$LocName == LOC] <- ParamVec["theta_het_case"]
    
    
    # 4. Perform the Likelihood Ratio Test for dAD detection
    # Likelihood of fit with all parameters shared:
    LikTot <- MAGE::pmf_betabinomMix(CurDF[CurDF$Outlier==0,]$ref_count, 
                                     CurDF[CurDF$Outlier==0,]$var_count, probshift = PiH0, SE, rho_rr, rho_vv, rho_rv, 
                                     theta_hom = ThetaHomH0, theta_het = ThetaHetH0)
    # Likelihood of fit with separate theta_het (calculated in two steps because of the 
    # different theta)
    LikCTRL <- MAGE::pmf_betabinomMix(CTRL_DF[CTRL_DF$Outlier==0,]$ref_count, 
                                      CTRL_DF[CTRL_DF$Outlier==0,]$var_count, probshift = ParamVec["probshift"], SE,
                                      ParamVec["pr"], ParamVec["pv"], ParamVec["prv"], theta_hom = ParamVec["theta_hom"],
                                      theta_het = ParamVec["theta_het_control"])
    LikCASE <- MAGE::pmf_betabinomMix(CASE_DF[CASE_DF$Outlier==0,]$ref_count, 
                                      CASE_DF[CASE_DF$Outlier==0,]$var_count, probshift = ParamVec["probshift"], SE,
                                      ParamVec["pr"], ParamVec["pv"], ParamVec["prv"], theta_hom = ParamVec["theta_hom"], 
                                      theta_het = ParamVec["theta_het_case"])
    lrtstat <- -2 * (sum(log(LikTot)) - sum(log(c(LikCTRL, LikCASE))))
    LRTpval <- pchisq(lrtstat, df = 1, lower.tail = F)
    dAD_res$LRTpval[dAD_res$LocName == LOC] <- LRTpval
    
    
    # 5. Fill out the results dataframe
    dAD_res$QualityCTRL[dAD_res$LocName == LOC] <- # spot bad quality data
      OUTfitCTRL$quality; dAD_res$QualityCASE[dAD_res$LocName == LOC] <-OUTfitCASE$quality 
    # Test HWE on both the control and tumor data:
    HWEtest_CTRL <- MAGE::HWE_chisquared(inbr = inbr, data = OUTfitCTRL_DH)
    dAD_res$HWECTRL[dAD_res$LocName == LOC] <- HWEtest_CTRL$PVAL
    HWEtest_CASE <- MAGE::HWE_chisquared(inbr = inbr, data = OUTfitCASE_DH)
    dAD_res$HWECASE[dAD_res$LocName == LOC] <- HWEtest_CASE$PVAL
    # Mean and median coverages:
    dAD_res$CovCTRL_mean[dAD_res$LocName == LOC] <- 
      mean(CTRL_DF[CTRL_DF$Outlier==0,]$ref_count + CTRL_DF[CTRL_DF$Outlier==0,]$var_count)
    dAD_res$CovCASE_mean[dAD_res$LocName == LOC] <- 
      mean(CASE_DF[CASE_DF$Outlier==0,]$ref_count + CASE_DF[CASE_DF$Outlier==0,]$var_count)
    dAD_res$CovCTRL_med[dAD_res$LocName == LOC] <- 
      median(CTRL_DF[CTRL_DF$Outlier==0,]$ref_count + CTRL_DF[CTRL_DF$Outlier==0,]$var_count)
    dAD_res$CovCASE_med[dAD_res$LocName == LOC] <- 
      median(CASE_DF[CASE_DF$Outlier==0,]$ref_count + CASE_DF[CASE_DF$Outlier==0,]$var_count)
    dAD_res$pr[dAD_res$LocName == LOC] <- ParamVec["pr"]
    dAD_res$prv[dAD_res$LocName == LOC] <- ParamVec["prv"]
    dAD_res$pv[dAD_res$LocName == LOC] <- ParamVec["pv"]
    
    dAD_res$ThetaHomH0[dAD_res$LocName == LOC] <- ThetaHomH0
    dAD_res$ThetaHomH1[dAD_res$LocName == LOC] <- ParamVec["theta_hom"]
  }
  
  dAD_res$HWECTRL[is.na(dAD_res$HWECTRL)] <- 
    dAD_res$HWECTRL[is.na(dAD_res$HWECTRL)] <- -1
  
  dAD_res$RhoHetH0 <- 1/((1/dAD_res$ThetaHetH0)+1)
  dAD_res$RhoHetCTRL <- 1/((1/dAD_res$ThetaHetCTRL)+1)
  dAD_res$RhoHetCASE <- 1/((1/dAD_res$ThetaHetCASE)+1)
  
  return(dAD_res)
}
