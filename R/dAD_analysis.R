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

dAD_analysis <- function(datas,
                         SE, inbr,
                         
                         # FitH0-parameter (regarding AB) doesn't need to be TRUE here
                         
                         allelefreq=0.5, dltaco = 10^-6, HWE = FALSE, p_InitEst = FALSE, p_inits = c(1/3, 1/3, 1/3), 
                         ThetaInits = "moment", ReEstThetas = "moment", NoSplitHom = TRUE, NoSplitHet = TRUE,
                         ResetThetaMin = 10^-10, ResetThetaMax = 10^-1, DistRob = "Cook", CookMargin = 5, 
                         LikEmpNum = 1000, LikMargin = 0, NumHetMin = 5, MaxOutFrac = 0.5, 
                         thetaTRY = c(10^-1, 10^-3, 10^-7), probshift_init = 0.5, FirstFewFixed = NULL, 
                         MemLim = 2048, Xtra = 5, epsabs = 0.001, MaxIt = 100
                         
){
  
  data <- datas[[1]]
  data_t <- datas[[2]]
  
  LocNames <- names(data)
  
  dAD_ResDF <- data.frame(Locus = LocNames, Gene = "", dAD_pval = 0,
                          
                          phi_rr = 0, phi_rv = 0, phi_vv = 0, Pi = 0, ThetaHetC = 0, ThetaHetT = 0, RhoC = 0, RhoT = 0, ThetaHom = 0,
                          phi_rr_H0 = 0, phi_rv_H0 = 0, phi_vv_H0 = 0, PiH0 = 0, ThetaHetH0 = 0, RhoH0 = 0, ThetaHomH0 = 0,
                          phi_rr_OnlyC = 0, phi_rv_OnlyC = 0, phi_vv_OnlyC = 0, PiOnlyC = 0, ThetaHetOnlyC = 0, ThetaHomOnlyC = 0,
                          phi_rr_OnlyT = 0, phi_rv_OnlyT = 0, phi_vv_OnlyT = 0, PiOnlyT = 0, ThetaHetOnlyT = 0, ThetaHomOnlyT = 0,
                          
                          NumHetC = 0, NumHetT = 0, RobFlagC = "", RobFlagT = "", HWEC = 0, HWET = 0, CovC = 0, CovT = 0, 
                          CovC_Med = 0, CovT_Med = 0, NumOutC = 0, NumOutT = 0, "nrep_H0" = 101, "nrep_H1" = 101, QualityC = "N", QualityT = "N"
                          
  )
  
  CData_list <- list(); TData_list <- list(); TData_list2 <- list(); CTData_listnames <- c()
  
  for(N in LocNames){
    
    CData <- data[[N]]
    if(nrow(CData)<10){
      next
    }
    TData <- data_t[[N]]
    if(is.null(TData)){
      next
    }
    if(nrow(TData)==0){
      next
    }
    if(median(CData$ref_count+CData$var_count) <= 3 | median(TData$ref_count+TData$var_count) <= 3){
      next
    }
    
    dAD_ResDF$Gene[dAD_ResDF$Locus == N] <- CData$Gene[1]
    C_DF <- data.frame("ref_count" = CData$ref_count, "var_count" = CData$var_count, "isCase" = 0)
    T_DF <- data.frame("ref_count" = TData$ref_count, "var_count" = TData$var_count, "isCase" = 1)
    
    # We FIRST remove outliers in both control- and tumordata, for which we perform separate beta-binomial fits on both
    OUTfitC <- maelstRom::EMfit_betabinom_robust(data_counts = C_DF, SE = SE, inbr = inbr, fitH0 = FALSE,
                                                 allelefreq=allelefreq, dltaco = dltaco, HWE = HWE, p_InitEst = p_InitEst, ThetaInits = ThetaInits, ReEstThetas = ReEstThetas, 
                                                 NoSplitHom = NoSplitHom, NoSplitHet = NoSplitHet, ResetThetaMin = ResetThetaMin, ResetThetaMax = ResetThetaMax, 
                                                 DistRob = DistRob, CookMargin = CookMargin, LikEmpNum = LikEmpNum, LikMargin = LikMargin, NumHetMin = NumHetMin, 
                                                 MaxOutFrac = MaxOutFrac, thetaTRY = thetaTRY, FirstFewFixed = FirstFewFixed, epsabs = epsabs)
    
    OUTfitT <- maelstRom::EMfit_betabinom_robust(data_counts = T_DF, SE = SE, inbr = inbr, fitH0 = FALSE,
                                                 allelefreq=allelefreq, dltaco = dltaco, HWE = HWE, p_InitEst = p_InitEst, ThetaInits = ThetaInits, ReEstThetas = ReEstThetas, 
                                                 NoSplitHom = NoSplitHom, NoSplitHet = NoSplitHet, ResetThetaMin = ResetThetaMin, ResetThetaMax = ResetThetaMax, 
                                                 DistRob = DistRob, CookMargin = CookMargin, LikEmpNum = LikEmpNum, LikMargin = LikMargin, NumHetMin = NumHetMin, 
                                                 MaxOutFrac = MaxOutFrac, thetaTRY = thetaTRY, FirstFewFixed = FirstFewFixed, epsabs = epsabs)
    
    
    OUTfitC_DH <- OUTfitC$data_hash
    OUTfitT_DH <- OUTfitT$data_hash
    C_DF$Outlier <- OUTfitC_DH$Outlier
    T_DF$Outlier <- OUTfitT_DH$Outlier
    
    dAD_ResDF$phi_rr_OnlyC[dAD_ResDF$Locus == N] <- OUTfitC$rho_rr
    dAD_ResDF$phi_rv_OnlyC[dAD_ResDF$Locus == N] <- OUTfitC$rho_rv
    dAD_ResDF$phi_vv_OnlyC[dAD_ResDF$Locus == N] <- OUTfitC$rho_vv
    dAD_ResDF$PiOnlyC[dAD_ResDF$Locus == N] <- OUTfitC$AB
    dAD_ResDF$ThetaHetOnlyC[dAD_ResDF$Locus == N] <- OUTfitC$theta_het
    dAD_ResDF$ThetaHomOnlyC[dAD_ResDF$Locus == N] <- OUTfitC$theta_hom
    
    dAD_ResDF$phi_rr_OnlyT[dAD_ResDF$Locus == N] <- OUTfitT$rho_rr
    dAD_ResDF$phi_rv_OnlyT[dAD_ResDF$Locus == N] <- OUTfitT$rho_rv
    dAD_ResDF$phi_vv_OnlyT[dAD_ResDF$Locus == N] <- OUTfitT$rho_vv
    dAD_ResDF$PiOnlyT[dAD_ResDF$Locus == N] <- OUTfitT$AB
    dAD_ResDF$ThetaHetOnlyT[dAD_ResDF$Locus == N] <- OUTfitT$theta_het
    dAD_ResDF$ThetaHomOnlyT[dAD_ResDF$Locus == N] <- OUTfitT$theta_hom
    
    #####
    CData$Outlier <- OUTfitC_DH$Outlier
    TData$Outlier <- OUTfitT_DH$Outlier
    #####
    
    CurDF <- rbind(C_DF, T_DF)
    dAD_ResDF$NumOutC[dAD_ResDF$Locus == N] <- sum(C_DF$Outlier)
    dAD_ResDF$NumOutT[dAD_ResDF$Locus == N] <- sum(T_DF$Outlier)
    
    
    
    # Perform ONE BETABINOMIAL FIT on the entire data, i.e. assuming the same theta-parameter, i.e. assuming NO dAD_Res_AllChr
    PiEstRes <- maelstRom::EMfit_betabinom(data_counts = CurDF[CurDF$Outlier == 0,], SE = SE, inbr = inbr, fitH0 = FALSE,
                                           
                                           allelefreq=allelefreq, dltaco = dltaco, HWE = HWE, p_InitEst = p_InitEst, probshift_init = probshift_init, 
                                           ThetaInits = ThetaInits, ReEstThetas = ReEstThetas, NoSplitHom = NoSplitHom, NoSplitHet = NoSplitHet,
                                           ResetThetaMin = ResetThetaMin, ResetThetaMax = ResetThetaMax, thetaTRY = thetaTRY,
                                           FirstFewFixed = FirstFewFixed, epsabs = 0.001, MaxIt = 100
                                           
    )
    
    PiEstRes_DF <- PiEstRes$data_hash
    if(sum(PiEstRes_DF$prv) < 10){
      next
    }
    PiH0 <- PiEstRes$AB
    phi_rr_H0 <- PiEstRes$rho_rr
    phi_rv_H0 <- PiEstRes$rho_rv
    phi_vv_H0 <- PiEstRes$rho_vv
    ThetaHomH0 <- PiEstRes$theta_hom
    ThetaHetH0 <- PiEstRes$theta_het
    
    dAD_ResDF$PiH0[dAD_ResDF$Locus == N] <- PiH0
    dAD_ResDF$ThetaHetH0[dAD_ResDF$Locus == N] <- ThetaHetH0
    dAD_ResDF$ThetaHomH0[dAD_ResDF$Locus == N] <- ThetaHomH0
    dAD_ResDF$phi_rr_H0[dAD_ResDF$Locus == N] <- phi_rr_H0
    dAD_ResDF$phi_rv_H0[dAD_ResDF$Locus == N] <- phi_rv_H0
    dAD_ResDF$phi_vv_H0[dAD_ResDF$Locus == N] <- phi_vv_H0
    
    # Perform a fit assuming different theta between control and tumor samples
    FullFit <- maelstRom::EMfit_betabinom_popcomb(data_counts = CurDF[CurDF$Outlier == 0,], SE = SE, inbr = inbr,
                                                  allelefreq=allelefreq, dltaco = dltaco, HWE = HWE, p_InitEst = p_InitEst, p_inits = p_inits, 
                                                  ThetaInits = ThetaInits, ReEstThetas = ReEstThetas, NoSplitHom = NoSplitHom, NoSplitHet = NoSplitHet,
                                                  ResetThetaMin = ResetThetaMin, ResetThetaMax = ResetThetaMax, thetaTRY = thetaTRY,
                                                  probshift_init = probshift_init, FirstFewFixed = FirstFewFixed, MemLim = MemLim, Xtra = Xtra, epsabs = epsabs, MaxIt = MaxIt)
    
    ParamVec <- FullFit$ParamVec
    # Nog prv en age-kolommen nodig...
    T_DF_2 <- TData
    #T_DF$sample <- TData$TCGA_sample_ID
    T_DF_2 <- T_DF_2[T_DF_2$Outlier==0,]
    T_DF_2$prr <- (FullFit$GenoDF)$prr[(nrow(FullFit$GenoDF)-nrow(T_DF_2)+1):nrow(FullFit$GenoDF)]
    T_DF_2$prv <- (FullFit$GenoDF)$prv[(nrow(FullFit$GenoDF)-nrow(T_DF_2)+1):nrow(FullFit$GenoDF)]
    T_DF_2$pvv <- (FullFit$GenoDF)$pvv[(nrow(FullFit$GenoDF)-nrow(T_DF_2)+1):nrow(FullFit$GenoDF)]
    
    dAD_ResDF$phi_rr[dAD_ResDF$Locus == N] <- ParamVec["pr"]
    dAD_ResDF$phi_rv[dAD_ResDF$Locus == N] <- ParamVec["prv"]
    dAD_ResDF$phi_vv[dAD_ResDF$Locus == N] <- ParamVec["pv"]
    
    dAD_ResDF$Pi[dAD_ResDF$Locus == N] <- ParamVec["probshift"]
    dAD_ResDF$ThetaHetC[dAD_ResDF$Locus == N] <- ParamVec["theta_het_control"]
    dAD_ResDF$ThetaHetT[dAD_ResDF$Locus == N] <- ParamVec["theta_het_case"]
    dAD_ResDF$ThetaHom[dAD_ResDF$Locus == N] <- ParamVec["theta_hom"]
    
    
    # Determine the log-likelihood for both models
    
    LikTot <- maelstRom::pmf_betabinomMix(CurDF[CurDF$Outlier==0,]$ref_count, CurDF[CurDF$Outlier==0,]$var_count, probshift = PiH0, SE, 
                                          phi_rr_H0, phi_vv_H0, phi_rv_H0, theta_hom = ThetaHomH0, theta_het = ThetaHetH0)
    LikC <- maelstRom::pmf_betabinomMix(C_DF[C_DF$Outlier==0,]$ref_count, C_DF[C_DF$Outlier==0,]$var_count, probshift = ParamVec["probshift"], SE, 
                                        ParamVec["pr"], ParamVec["pv"], ParamVec["prv"], theta_hom = ParamVec["theta_hom"], theta_het = ParamVec["theta_het_control"])
    LikT <- maelstRom::pmf_betabinomMix(T_DF[T_DF$Outlier==0,]$ref_count, T_DF[T_DF$Outlier==0,]$var_count, probshift = ParamVec["probshift"], SE, 
                                        ParamVec["pr"], ParamVec["pv"], ParamVec["prv"], theta_hom = ParamVec["theta_hom"], theta_het = ParamVec["theta_het_case"])
    lrtstat <- -2 * (sum(log(LikTot)) - sum(log(c(LikC, LikT))))
    LRTpval <- pchisq(lrtstat, df = 1, lower.tail = F)
    
    dAD_ResDF$dAD_pval[dAD_ResDF$Locus == N] <- LRTpval
    
    dAD_ResDF$NumHetC[dAD_ResDF$Locus == N] <- sum(OUTfitC_DH$prv)
    dAD_ResDF$NumHetT[dAD_ResDF$Locus == N] <- sum(OUTfitT_DH$prv)
    dAD_ResDF$RobFlagC[dAD_ResDF$Locus == N] <- OUTfitC$RobFlag
    dAD_ResDF$RobFlagT[dAD_ResDF$Locus == N] <- OUTfitT$RobFlag
    dAD_ResDF$QualityC[dAD_ResDF$Locus == N] <- OUTfitC$quality
    dAD_ResDF$QualityT[dAD_ResDF$Locus == N] <- OUTfitT$quality
    HWEtest_C <- maelstRom::HWE_chisquared(inbr = inbr, data = OUTfitC_DH)
    dAD_ResDF$HWEC[dAD_ResDF$Locus == N] <- HWEtest_C$HWEpval
    HWEtest_T <- maelstRom::HWE_chisquared(inbr = inbr, data = OUTfitT_DH)
    dAD_ResDF$HWET[dAD_ResDF$Locus == N] <- HWEtest_T$HWEpval
    dAD_ResDF$CovC[dAD_ResDF$Locus == N] <- mean(C_DF[C_DF$Outlier==0,]$ref_count + C_DF[C_DF$Outlier==0,]$var_count)
    dAD_ResDF$CovT[dAD_ResDF$Locus == N] <- mean(T_DF[T_DF$Outlier==0,]$ref_count + T_DF[T_DF$Outlier==0,]$var_count)
    dAD_ResDF$CovC_Med[dAD_ResDF$Locus == N] <- median(C_DF[C_DF$Outlier==0,]$ref_count + C_DF[C_DF$Outlier==0,]$var_count)
    dAD_ResDF$CovT_Med[dAD_ResDF$Locus == N] <- median(T_DF[T_DF$Outlier==0,]$ref_count + T_DF[T_DF$Outlier==0,]$var_count)
    
    dAD_ResDF$nrep_H0[dAD_ResDF$Locus == N] <- PiEstRes$nrep
    dAD_ResDF$nrep_H1[dAD_ResDF$Locus == N] <- FullFit$nrep
    
    CData_list <- c(CData_list, list(CData))
    TData_list <- c(TData_list, list(TData))
    TData_list2 <- c(TData_list2, list(T_DF_2))
    
    CTData_listnames <- c(CTData_listnames, N)
  }
  
  dAD_ResDF <- dAD_ResDF[dAD_ResDF$CovC!=0,]
  dAD_ResDF$RhoH0 <- 1/((1/dAD_ResDF$ThetaHetH0)+1)
  dAD_ResDF$RhoC <- 1/((1/dAD_ResDF$ThetaHetC)+1)
  dAD_ResDF$RhoT <- 1/((1/dAD_ResDF$ThetaHetT)+1)
  
  names(CData_list) <- CTData_listnames
  names(TData_list) <- CTData_listnames
  names(TData_list2) <- CTData_listnames
  
  return(list(dAD_ResDF, CData_list, TData_list, TData_list2))
  
}
