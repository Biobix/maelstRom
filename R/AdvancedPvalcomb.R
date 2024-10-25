
#' @export

AdvancedPvalcomb <- function(DATAS, SE, nperm = 10000){
  
  # A helper function for inverting a matrix, providing some error handling:
  InvHelper <- function(MatOI){
    Out <- tryCatch({matlib::inv(MatOI)}, error = function(e) NULL)  
    if(is.null(Out)){
      return(NA)
    }else{
      return(Out)
    }
  }
  
  # Vectors for storing the per-gene combined p-values: Fisher, Lancaster, correlation-corrected Lancaster and Caucy.
  # CLPV_OffDiagSum_vec is a measure of the correlations between test statistics as used in correlation-corrected Lancaster,
  # but is scales with the number of tests being combined
  FisherPVAL_vec <- LancasterPVAL_vec <-  CorLancasterPVAL_vec <-  CauchyPVAL_vec <- c()
  CLPV_OffDiagSum_vec <- c()
  
  goi <- DATAS[[1]]
  ResDF_OIs <- DATAS[[2]]
  ControlOIs <- DATAS[[3]]
  CaseOIs <- DATAS[[4]]
  
  for(GT in 1:length(goi)){

    GENE <- goi[GT]
    
    ResDF_OI <- ResDF_OIs[[GT]]
    LocsOI <- ResDF_OI$Locus
    ControlOI <- ControlOIs[[GT]]
    CaseOI <- CaseOIs[[GT]]
    
    if(nrow(ResDF_OI) > 1){ # We need to actually have more than one p-value in order to combine anything.
      
      # I COULD use regular old (unweighted) FISHER
      
      FisherPVAL <- pchisq(-2*sum(log(ResDF_OI$dAD_pval)), df = 2*length(ResDF_OI$dAD_pval), lower.tail = FALSE) 
      
      # GOOD EXPLANATION of the Lancaster procedure: https://www.frontiersin.org/articles/10.3389/fgene.2014.00032/full
      # SOME DISCUSSION ABOUT THE IDEAL WEIGHTING SCHEME; this publication suggests one that keeps the degrees of freedom
      # equal to those in the (unweighted) Fisher procedure, that being 2n (2 times the number of hypotheses):
      # https://www.nature.com/articles/s41598-021-86465-y#Sec16
      
      # So, yeah... First up: Lancaster. Papers suggest weights ~ sample size. Here, of course, we also have coverage to take into account.
      # So I'm going for sample size * coverage. You could argue for or against putting a square root in there,
      # but for now I'm going to leave it at just the product of these terms
      # ... I'll take the minimum of {controls, cases} as a weight...
      
      LancasterWeights <- c()
      for(i in 1:nrow(ResDF_OI)){
        LancasterWeights <- c(LancasterWeights, min(c(ResDF_OI[i,"NumHetC"]*ResDF_OI[i,"CovC_Med"], ResDF_OI[i,"NumHetT"]*ResDF_OI[i,"CovT_Med"])))
      }
      
      LancasterWeights <- (LancasterWeights/sum(LancasterWeights))*length(LancasterWeights) # Need to sum to number of hypotheses
      
      LancasterTS <- c() # The test statistic
      for(i in 1:nrow(ResDF_OI)){
        LancasterTS <- c(LancasterTS, qgamma(ResDF_OI[i,"dAD_pval"], shape = LancasterWeights[i], scale = 2, lower.tail = FALSE))
      }
      
      LancasterPVAL <- pchisq(sum(LancasterTS), df = 2*length(LancasterTS), lower.tail = FALSE) 
    
      # First one down, two more to go.
      # Correlation-corrected Lancaster... This one's the most time-consuming, as it involves permutation-skedaddle.
      # IMPROVEMENT (or at least speedup): we'll use score tests!
      
      CurDF_List <- vector(mode = "list", length=length(LocsOI))
      PermCheckDF_List <- vector(mode = "list", length=length(LocsOI))
      
      for(i in 1:length(LocsOI)){
        
        CurDF <- rbind((ControlOI[[LocsOI[i]]])[,c("Sample", "Outlier", "ref_count", "var_count")],
                       (CaseOI[[LocsOI[i]]])[,c("Sample", "Outlier", "ref_count", "var_count")])
        # Make sure that "Sample" is actually a unique column!
        # If you share sample-ids between controls and cases, e.g. add something to these ids to make them distinct.
        
        # We can use the H0 fit for the score test:
        PiH0 <- ResDF_OI[i, "PiH0"]
        rho_rr <- ResDF_OI[i, "phi_rr_H0"]
        rho_rv <- ResDF_OI[i, "phi_rv_H0"]
        rho_vv <- ResDF_OI[i, "phi_vv_H0"]
        ThetaHomH0 <- ResDF_OI[i, "ThetaHomH0"]
        ThetaHetH0 <- ResDF_OI[i, "ThetaHetH0"]
        
        LikTot <- maelstRom::pmf_betabinomMix(CurDF[CurDF$Outlier==0,]$ref_count, CurDF[CurDF$Outlier==0,]$var_count, probshift = PiH0, SE, 
                                         rho_rr, rho_vv, rho_rv, theta_hom = ThetaHomH0, theta_het = ThetaHetH0)
        
        # Partial derivatives to all parameters;
        # remember that one of {phi_rr, phi_rv, phi_vv} is not actually a parameter as the three are related:
        # phi_rr + phi_rv + phi_vv
        # so one of these can be written as a function of the other two.
        Varis <- c("probshift", "theta_hom", "theta_het", "prr", "prv")
        NewColNames <- c()
        for(V in Varis){
          NewColNames <- c(NewColNames, paste0("D", V))
        }
        for(V in Varis){
          for(V2 in Varis){
            NewColNames <- c(NewColNames, paste0("D", V, "D", V2))
          }
        }
        
        for(V in Varis){
          
          CurDF <- cbind(CurDF, rep(NA, nrow(CurDF)))
          TssVec <- maelstRom::BetaBinomMix_LLDeriv(CurDF[CurDF$Outlier==0,]$ref_count, CurDF[CurDF$Outlier==0,]$var_count, probshift = PiH0, 
                                               SE = SE, prr = rho_rr, pvv = rho_vv, prv = rho_rv, 
                                               theta_hom = ThetaHomH0, theta_het = ThetaHetH0, Der1 = V, Der2 = NULL)
          CurDF[CurDF$Outlier==0, ncol(CurDF)] <- TssVec
          
        }
        for(V in Varis){
          for(V2 in Varis){
            
            CurDF <- cbind(CurDF, rep(NA, nrow(CurDF)))
            TssVec <- maelstRom::BetaBinomMix_LLDeriv(CurDF[CurDF$Outlier==0,]$ref_count, CurDF[CurDF$Outlier==0,]$var_count, probshift = PiH0, 
                                                 SE = SE, prr = rho_rr, pvv = rho_vv, prv = rho_rv, 
                                                 theta_hom = ThetaHomH0, theta_het = ThetaHetH0, Der1 = V, Der2 = V2)
            CurDF[CurDF$Outlier==0, ncol(CurDF)] <- TssVec
            
          }
        }
        
        colnames(CurDF)[5:ncol(CurDF)] <- NewColNames
        
        CurDF_List[[i]] <- CurDF
        
      }
      
      ControlSamps <- c(); CaseSamps <- c()
      for(i in 1:length(LocsOI)){
        ControlSamps <- c(ControlSamps, ControlOI[[LocsOI[i]]]$Sample)
        CaseSamps <- c(CaseSamps, CaseOI[[LocsOI[i]]]$Sample)
      }
      ControlSamps <- unique(ControlSamps); CaseSamps <- unique(CaseSamps)
      AllSamps <- c(ControlSamps, CaseSamps)
      SampsLabels <- c(rep("Control", length(ControlSamps)), rep("Case", length(CaseSamps)))
      PermCheckDF <- data.frame("Sample" = AllSamps, "Label" = SampsLabels)
      
      for(i in 1:length(LocsOI)){
        prr_control <- maelstRom::dBetaBinom(ControlOI[[LocsOI[i]]]$ref_count, ControlOI[[LocsOI[i]]]$ref_count + ControlOI[[LocsOI[i]]]$var_count,
                                        pi = 1 - SE, theta = ResDF_OI[i, "ThetaHomOnlyC"], LOG = FALSE)
        pvv_control <- maelstRom::dBetaBinom(ControlOI[[LocsOI[i]]]$var_count, ControlOI[[LocsOI[i]]]$ref_count + ControlOI[[LocsOI[i]]]$var_count,
                                        pi = 1 - SE, theta = ResDF_OI[i, "ThetaHomOnlyC"], LOG = FALSE)
        prv_control <- maelstRom::dBetaBinom(ControlOI[[LocsOI[i]]]$ref_count, ControlOI[[LocsOI[i]]]$ref_count + ControlOI[[LocsOI[i]]]$var_count,
                                        pi = ResDF_OI[i, "PiOnlyC"], theta = ResDF_OI[i, "ThetaHetOnlyC"], LOG = FALSE)
        pdata_control <- rowSums(cbind(prr_control, prv_control, pvv_control)) 
        if (any(pdata_control==0)){ 
          ProblemCases <- which(pdata_control==0)
          for(case in ProblemCases){
            var_part<-ControlOI[[LocsOI[i]]]$var_count[case]
            ref_part<-ControlOI[[LocsOI[i]]]$ref_count[case]
            spr_part<-maelstRom::dBetaBinom(ref_part, ref_part+var_part, pi = 1-SE, theta = ResDF_OI[i, "ThetaHomOnlyC"], LOG = TRUE)
            spv_part<-maelstRom::dBetaBinom(var_part, ref_part+var_part, pi = 1-SE, theta = ResDF_OI[i, "ThetaHomOnlyC"], LOG = TRUE)
            sprv_part<-maelstRom::dBetaBinom(ref_part, ref_part+var_part, pi = ResDF_OI[i, "PiOnlyC"], theta = ResDF_OI[i, "ThetaHetOnlyC"], LOG = TRUE)
            spvec<-c(spr_part, sprv_part, spv_part)
            if(spr_part==max(spvec)){
              prr_control[case]<-1
              pdata_control[case]<-1
            } else if(sprv_part==max(spvec)){
              prv_control[case]<-1
              pdata_control[case]<-1
            } else{
              pvv_control[case]<-1
              pdata_control[case]<-1
            }
          }
        }
        prr_case <- maelstRom::dBetaBinom(CaseOI[[LocsOI[i]]]$ref_count, CaseOI[[LocsOI[i]]]$ref_count + CaseOI[[LocsOI[i]]]$var_count,
                                     pi = 1 - SE, theta = ResDF_OI[i, "ThetaHomOnlyT"], LOG = FALSE)
        pvv_case <- maelstRom::dBetaBinom(CaseOI[[LocsOI[i]]]$var_count, CaseOI[[LocsOI[i]]]$ref_count + CaseOI[[LocsOI[i]]]$var_count,
                                     pi = 1 - SE, theta = ResDF_OI[i, "ThetaHomOnlyT"], LOG = FALSE)
        prv_case <- maelstRom::dBetaBinom(CaseOI[[LocsOI[i]]]$ref_count, CaseOI[[LocsOI[i]]]$ref_count + CaseOI[[LocsOI[i]]]$var_count,
                                     pi = ResDF_OI[i, "PiOnlyT"], theta = ResDF_OI[i, "ThetaHetOnlyT"], LOG = FALSE)
        pdata_case <- rowSums(cbind(prr_case, prv_case, pvv_case)) 
        if (any(pdata_case==0)){ 
          ProblemCases <- which(pdata_case==0)
          for(case in ProblemCases){
            var_part<-CaseOI[[LocsOI[i]]]$var_count[case]
            ref_part<-CaseOI[[LocsOI[i]]]$ref_count[case]
            spr_part<-maelstRom::dBetaBinom(ref_part, ref_part+var_part, pi = 1-SE, theta = ResDF_OI[i, "ThetaHomOnlyT"], LOG = TRUE)
            spv_part<-maelstRom::dBetaBinom(var_part, ref_part+var_part, pi = 1-SE, theta = ResDF_OI[i, "ThetaHomOnlyT"], LOG = TRUE)
            sprv_part<-maelstRom::dBetaBinom(ref_part, ref_part+var_part, pi = ResDF_OI[i, "PiOnlyT"], theta = ResDF_OI[i, "ThetaHetOnlyT"], LOG = TRUE)
            spvec<-c(spr_part, sprv_part, spv_part)
            if(spr_part==max(spvec)){
              prr_case[case]<-1
              pdata_case[case]<-1
            } else if(sprv_part==max(spvec)){
              prv_case[case]<-1
              pdata_case[case]<-1
            } else{
              pvv_case[case]<-1
              pdata_case[case]<-1
            }
          }
        }
        ControlOI[[LocsOI[i]]]$prv_OnlyTC_ROB <- prv_control/pdata_control
        CaseOI[[LocsOI[i]]]$prv_OnlyTC_ROB <- prv_case/pdata_case
        ControlOI[[LocsOI[i]]]$prv_OnlyTC_ROB[ControlOI[[LocsOI[i]]]$Outlier == 1] <- 0
        CaseOI[[LocsOI[i]]]$prv_OnlyTC_ROB[CaseOI[[LocsOI[i]]]$Outlier == 1] <- 0
        PermCheckDF <- rbind( merge(PermCheckDF[PermCheckDF$Label == "Control",], (ControlOI[[LocsOI[i]]])[,c("Sample", "prv_OnlyTC_ROB")],
                                    by = "Sample", all = TRUE),
                              merge(PermCheckDF[PermCheckDF$Label == "Case",], (CaseOI[[LocsOI[i]]])[,c("Sample", "prv_OnlyTC_ROB")],
                                    by = "Sample", all = TRUE))
        colnames(PermCheckDF)[colnames(PermCheckDF) == "prv_OnlyTC_ROB"] <- LocsOI[i]
        
      }
      
      PermCheckDF[is.na(PermCheckDF)] <- 0
      
      
      
      for(m in 1:nperm){
        
        GoOn <- TRUE
        while(GoOn){
          
          dAD_pvalVec <- c()
          
          PermTry <- TRUE
          while(PermTry){
            PermCheckDF$Perm <- sample(PermCheckDF$Label)
            
            MinNumHetControl <-  min(colSums(PermCheckDF[PermCheckDF$Perm == "Control", colnames(PermCheckDF)%in%LocsOI, drop = FALSE]))
            MinNumHetCase <- min(colSums(PermCheckDF[PermCheckDF$Perm == "Case", colnames(PermCheckDF)%in%LocsOI, drop = FALSE]))
            
            if(MinNumHetControl > 6 & MinNumHetCase > 6){
              PermTry <- FALSE
            }
            # Hopefully we don't get stuck in this dumb while-loop forever... but we shouldn't, really,
            # 8 heterozygotes is quite a bit less than the 12 per control and case we filtered on earlier...
          }
          
          for(i in 1:length(LocsOI)){
            
            CurDF <- CurDF_List[[i]]
            
            CurDF <- merge(CurDF, PermCheckDF[, c("Sample", "Perm")], by = "Sample")
            CurDF$isCase <- as.numeric(CurDF$Perm == "Case")
            
            VarisX <- c("probshift", "theta_hom", "theta_het_control", "theta_het_case", "prr", "prv")
            
            CurDF$Dtheta_het_control <- CurDF$Dtheta_het
            CurDF$Dtheta_het_control[CurDF$isCase == 0] <- 0
            CurDF$Dtheta_het_case <- CurDF$Dtheta_het
            CurDF$Dtheta_het_case[CurDF$isCase == 1] <- 0
            
            EfficientScoresVec <- c()
            for(V in VarisX){
              EfficientScoresVec <- c(EfficientScoresVec, sum(CurDF[CurDF$Outlier == 0, paste0("D", V)]))
            }
            
            ObservedFisher <- c()
            for(V in VarisX){
              if(V == "theta_het_control" | V == "theta_het_case"){
                VT <- "theta_het"
              }else{
                VT <- V
              }
              for(V2 in VarisX){
                if(V2 == "theta_het_control" | V2 == "theta_het_case"){
                  V2T <- "theta_het"
                }else{
                  V2T <- V2
                }
                
                if( (V=="theta_het_control" & V2=="theta_het_case")|(V=="theta_het_case" & V2=="theta_het_control") ){
                  ObservedFisher <- c(ObservedFisher, 0)
                } else if (V=="theta_het_control" | V2=="theta_het_control"){
                  ObservedFisher <- c(ObservedFisher, -sum(CurDF[CurDF$Outlier == 0 & CurDF$isCase == 0, paste0("D", VT, "D", V2T)]))
                } else if (V=="theta_het_case" | V2=="theta_het_case"){
                  ObservedFisher <- c(ObservedFisher, -sum(CurDF[CurDF$Outlier == 0 & CurDF$isCase == 1, paste0("D", VT, "D", V2T)]))
                } else {
                  ObservedFisher <- c(ObservedFisher, -sum(CurDF[CurDF$Outlier == 0, paste0("D", VT, "D", V2T)]))
                }
              }
            }
            ObservedFisher <- matrix(ObservedFisher, ncol = length(VarisX))
            
            Finv <- InvHelper(ObservedFisher)
            
            if(length(Finv)==1 && is.na(Finv)){
              break
            }else{
              ScoreStat <- matrix(EfficientScoresVec, nrow = 1) %*% Finv %*% matrix(EfficientScoresVec, ncol = 1)
              
              dAD_pval <- pchisq((ScoreStat), df = 1, lower.tail = F)
              
              dAD_pvalVec <- c(dAD_pvalVec, dAD_pval)
              
              GoOn <- FALSE
            }
            
          }
          
        }
        
        if(!GoOn){
          if(m==1){
            dAD_pvalMat <- dAD_pvalVec
          }else{
            dAD_pvalMat <- rbind(dAD_pvalMat, dAD_pvalVec)
          }
        }
        
      }
      
      
      # VarCovMat <- cor(dAD_pvalMat, method = "pearson")
      # I'll need to calculate the variance covariance matrix on transformed p-values...
      
      # The LancasterWeights remain the same...
      E_T <- 2*length(LocsOI) # Which equals two times the sum of the Lancaster weights
      # BE CAREFUL WITH THIS; a lot of publications use weight/2 in the degrees of freedom, I don't,
      # SO I have to multiply by two whenever they use the regular "w" on itself
      
      dAD_pvalMat_TR <- dAD_pvalMat
      dAD_pvalMat_TR[dAD_pvalMat_TR < 1e-4] <- NA
      for(Li in 1:length(LocsOI)){
        dAD_pvalMat_TR[,Li] <- qgamma(sapply(dAD_pvalMat_TR[,Li], FUN = function(x){max(x, 1e-300)}, simplify = TRUE), shape = LancasterWeights[Li], scale = 2, lower.tail = FALSE)
      }
      
      NAfracs <- c()
      for(MI in 1:ncol(dAD_pvalMat_TR)){
        NAfracs <- c(NAfracs,sum(is.na(dAD_pvalMat_TR[,MI]))/nrow(dAD_pvalMat_TR))
      }
      
      if(!any(NAfracs > 0.8)){
        # Not CORrelation but COVariance, here... 
        VarCovMat <- cov(dAD_pvalMat_TR, method = "pearson", use = "pairwise.complete.obs")
        VarCovMat[VarCovMat < 0] <- 0   # Negatives can occur, but... just aren't logical?
        
        CorMat <- VarCovMat
        for(iii in 1:ncol(dAD_pvalMat_TR)){
          for(jjj in 1:ncol(dAD_pvalMat_TR)){
            CorMat[iii,jjj] <- sqrt(LancasterWeights[iii]*4)*sqrt(LancasterWeights[jjj]*4)
          }
        }
        
        VarCovMat <- pmin(VarCovMat, CorMat)
      }else{
        MU <- which(NAfracs > 0.8)
        dAD_pvalMat_TR[,MU] <- NA
        
        VarCovMat <- cov(dAD_pvalMat_TR, method = "pearson", use = "pairwise.complete.obs")
        VarCovMat[VarCovMat < 0] <- 0   # Negatives can occur, but... just aren't logical?
        CorMat <- VarCovMat
        for(iii in 1:ncol(dAD_pvalMat_TR)){
          for(jjj in 1:ncol(dAD_pvalMat_TR)){
            CorMat[iii,jjj] <- sqrt(LancasterWeights[iii]*4)*sqrt(LancasterWeights[jjj]*4)
          }
        }
        
        for(iii in 1:nrow(VarCovMat)){
          for(jjj in MU){
            VarCovMat[iii,jjj] <- CorMat[iii,jjj]
            VarCovMat[jjj,iii] <- CorMat[jjj,iii]
          }
        }
        
        VarCovMat <- pmin(VarCovMat, CorMat)
      }
      
      Var_T <- 2*2*length(LocsOI) + 2*sum(VarCovMat[upper.tri(VarCovMat)])
      
      # The test statistic, LancasterTS, also remains the same...
      # We just need to multiply it by a constant to then compare it to a chi square distribution
      
      CL_v <- (2*E_T^2)/Var_T
      CL_c <- CL_v/E_T
      CorLancasterPVAL <- pchisq(CL_c*sum(LancasterTS), df = CL_v, lower.tail = FALSE) 
      
      CLPV_DiagSum <- sum(VarCovMat[upper.tri(VarCovMat)])
      
      # Only Cauchy remains... this one is relatively straightforward
      
      # For Cauchy, the SUM OF THE WEIGHTS MUST BE EQUAL TO 1
      # Here, I take them proportional to the Lancaster weights, which may or may not be a good weighting scheme. Further theory and/or testing might be needed...
      
      CauchyWeights <- LancasterWeights/sum(LancasterWeights)
      CauchyPVAL <- pcauchy(sum(CauchyWeights*tan((0.5-ResDF_OI$dAD_pval)*pi)), lower.tail = FALSE)
      
    }else{
      
      # Only one p-value; no combination necessary.
      FisherPVAL <- LancasterPVAL <-  CorLancasterPVAL <-  CauchyPVAL <- ResDF_OI$dAD_pval
      CLPV_DiagSum <- 0
      
    }
    
    FisherPVAL_vec <- c(FisherPVAL_vec, FisherPVAL)
    LancasterPVAL_vec <- c(LancasterPVAL_vec, LancasterPVAL) 
    CorLancasterPVAL_vec <- c(CorLancasterPVAL_vec, CorLancasterPVAL)
    CauchyPVAL_vec <- c(CauchyPVAL_vec, CauchyPVAL)
    CLPV_OffDiagSum_vec <- c(CLPV_OffDiagSum_vec, CLPV_DiagSum)
    
  }
  
  return(list(FisherPVAL_vec, LancasterPVAL_vec, CorLancasterPVAL_vec, CauchyPVAL_vec, CLPV_OffDiagSum_vec))
  
}
