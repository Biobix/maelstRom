#' AnnotateLater 
#' @export

AD_correlator <- function(CountDF, dAD_fitres, SE, CustomCor = NULL, MinCount = 0, method = "spearman"){
  
  CountDF$GT<- "none"
  
  CountDF$EX <- NA
  
  for(AAi in 1:nrow(CountDF)){
    
    EXTR <- min(maelstRom::pBetaBinom(CountDF$ref_count[AAi],
                                 CountDF$ref_count[AAi] + CountDF$var_count[AAi],
                                 dAD_fitres$Pi, dAD_fitres$ThetaHetT),
                maelstRom::pBetaBinom(CountDF$ref_count[AAi],
                                 CountDF$ref_count[AAi] + CountDF$var_count[AAi],
                                 dAD_fitres$Pi, dAD_fitres$ThetaHetT, lower.tail = FALSE) + 
                  maelstRom::dBetaBinom(CountDF$ref_count[AAi],
                                   CountDF$ref_count[AAi] + CountDF$var_count[AAi],
                                   dAD_fitres$Pi, dAD_fitres$ThetaHetT))
    
    CountDF$EX[AAi] <- EXTR
    
    GT <- which(c(dAD_fitres$phi_rr*maelstRom::dBetaBinom(CountDF$ref_count[AAi],
                                                     CountDF$ref_count[AAi] + CountDF$var_count[AAi],
                                                     1-SE, dAD_fitres$ThetaHom),
                  dAD_fitres$phi_rv*maelstRom::dBetaBinom(CountDF$ref_count[AAi],
                                                     CountDF$ref_count[AAi] + CountDF$var_count[AAi],
                                                     dAD_fitres$Pi, dAD_fitres$ThetaHetT),
                  dAD_fitres$phi_vv*maelstRom::dBetaBinom(CountDF$ref_count[AAi],
                                                     CountDF$ref_count[AAi] + CountDF$var_count[AAi],
                                                     SE, dAD_fitres$ThetaHom)) == 
                  max(c(dAD_fitres$phi_rr*maelstRom::dBetaBinom(CountDF$ref_count[AAi],
                                                           CountDF$ref_count[AAi] + CountDF$var_count[AAi],
                                                           1-SE, dAD_fitres$ThetaHom),
                        dAD_fitres$phi_rv*maelstRom::dBetaBinom(CountDF$ref_count[AAi],
                                                           CountDF$ref_count[AAi] + CountDF$var_count[AAi],
                                                           dAD_fitres$Pi, dAD_fitres$ThetaHetT),
                        dAD_fitres$phi_vv*maelstRom::dBetaBinom(CountDF$ref_count[AAi],
                                                           CountDF$ref_count[AAi] + CountDF$var_count[AAi],
                                                           SE, dAD_fitres$ThetaHom))))
    
    if(any(GT == 2)){
      CountDF$GT[AAi]<- "rv"
    }else if(any(GT == 1)){
      CountDF$GT[AAi]<- "rr"
    }else{
      CountDF$GT[AAi]<- "vv"
    }
    
  }
  
  CountDF$EX <- -log(CountDF$EX)
  
  CountDF <- CountDF[!is.na(CountDF$EX),]
  
  CountDF <- CountDF[CountDF$GT == "rv" &
                       (CountDF$ref_count + CountDF$var_count > MinCount) &
                       CountDF$Outlier != 1,]
  
  if(is.null(CustomCor)){
    CountDF$ToCorr <- CountDF$ref_count + CountDF$var_count
  }else{
    CountDF$ToCorr <- CustomCor
  }
  
  CountDF <- CountDF[!is.na(CountDF$ToCorr),]
  
  if(nrow(CountDF) <= 1){
    ReturnObj <- c(NA, NA)
    names(ReturnObj) <- c("correlation", "p-value")
    return(ReturnObj)
  }else{
    ReturnObj <- c()
    ReturnObj <- c(ReturnObj, cor(CountDF$EX, CountDF$ToCorr, method = method))
    ReturnObj <- c(ReturnObj, cor.test(CountDF$EX, CountDF$ToCorr, method = method)$p.value)
    names(ReturnObj) <- c("correlation", "p-value")
    return(ReturnObj)
  }
  
}
