
#' @export

AdvancedFitter <- function(ThetaForm = Theta ~ 1, PiForm = NULL, SNPdata, RefCol = "ref_count", VarCol = "var_count",
                           WeightCol = NULL, Pi_start, Theta_start, PiLink = "identity", ThetaLink = "identity", center = TRUE,
                           ResetThetaMin = 10^-10, ResetThetaMax = 10^-1){
  
  # I'm writing this function very prototype-like for now. In the future, it may be adviseable to handle input (the formulas) better,
  # and maybe put some restrictions on the parameter(fit)s. 
  
  if(Theta_start > ResetThetaMax){
    Theta_start <- ResetThetaMax
  } else if (Theta_start < ResetThetaMin){
    Theta_start <- ResetThetaMin
  }
  
  if(PiLink == "identity"){
    PiTransFun <- I
    PiTransInv <- I
  }else if(PiLink == "log"){
    PiTransFun <- log
    PiTransInv <- exp
  }else if(PiLink == "sqrt"){
    PiTransFun <- sqrt
    PiTransInv <- function(x){return(x^2)}
  }
  
  if(ThetaLink == "identity"){
    ThetaTransFun <- I
    ThetaTransInv <- I
  }else if(ThetaLink == "log"){
    ThetaTransFun <- log
    ThetaTransInv <- exp
  }else if(ThetaLink == "sqrt"){
    ThetaTransFun <- sqrt
    ThetaTransInv <- function(x){return(x^2)}
  }
  
  ThetaTerms <- terms(ThetaForm)
  if(is.null(PiForm)){
    PiFix <- Pi_start
    Pi_start <- NULL
    PiTerms <- NULL
  }else{
    PiFix <- NULL
    PiTerms <- terms(PiForm)
  }
  
  
  ThetaLabs <- attr(ThetaTerms, "term.labels")
  if(is.null(PiForm)){
    PiLabs <- NULL
  }else{
    PiLabs <- attr(PiTerms, "term.labels")
  }
  
  if(!is.null(Pi_start)){
    Pi_start <- gtools::logit(Pi_start, max = 1-10^-16, min = 10^-16)
  }

  paramvec <- c(Pi_start, log(Theta_start), rep(0, length(PiLabs)+length(ThetaLabs)))
  
  AdvancedLik <- function(paramvec, SNPdata, PiLabs, ThetaLabs, RefCol, VarCol, WeightCol, PiTransFun, PiTransInv, ThetaTransFun, ThetaTransInv, center, PiFix = NULL){
    
    Fixer <- 0
    if(!is.null(PiFix)){
      Fixer <- 1
    }else{
      pivec <- rep(PiTransFun(gtools::inv.logit(paramvec[1-Fixer], max = 1-10^-16, min = 10^-16)), nrow(SNPdata))
    }
    
    thetavec <- rep(ThetaTransFun(exp(paramvec[2-Fixer])), nrow(SNPdata))
    
    teller <- 3-Fixer
    for(PL in PiLabs){
      if(center){
        pivec <- (pivec + paramvec[teller]*(SNPdata[,PL]-mean(SNPdata[,PL])))
        teller <- teller+1
      }else{
        pivec <- (pivec + paramvec[teller]*SNPdata[,PL])
        teller <- teller+1
      }
    }
    if(is.null(PiFix)){
      pivec <- PiTransInv(pivec)
    }
    
    for(TL in ThetaLabs){
      if(center){
        thetavec <- (thetavec + paramvec[teller]*(SNPdata[,TL]-mean(SNPdata[,TL])))
        teller <- teller+1
      }else{
        thetavec <- (thetavec + paramvec[teller]*SNPdata[,TL])
        teller <- teller+1
      }
    }
    
    if(identical(sqrt, ThetaTransFun)){
      if(any(thetavec < 0)){
        return(10000000*sum(thetavec < 0))
      }
      thetavec[thetavec >= 0] <- ThetaTransInv(thetavec)
    }else if(identical(I, ThetaTransFun)){
      if(any(thetavec < 0)){
        return(10000000*sum(thetavec < 0))
      }
      thetavec <- ThetaTransInv(thetavec)
    }else{
      thetavec <- ThetaTransInv(thetavec)
    }
    
    if(!is.null(PiFix)){
      LikVec <- c()
      for(i in 1:length(thetavec)){
        if(is.null(WeightCol)){
          LikVec <- c(LikVec, maelstRom::dBetaBinom(ms = SNPdata[i, RefCol], ns = SNPdata[i, RefCol]+SNPdata[i, VarCol],pi = PiFix, theta = thetavec[i], LOG = TRUE))
        }else{
          LikVec <- c(LikVec, SNPdata[i, WeightCol]*maelstRom::dBetaBinom(ms = SNPdata[i, RefCol], ns = SNPdata[i, RefCol]+SNPdata[i, VarCol],pi = PiFix, theta = thetavec[i], LOG = TRUE))
        }
      }
    }else{
      LikVec <- c()
      for(i in 1:length(thetavec)){
        if(is.null(WeightCol)){
          LikVec <- c(LikVec, maelstRom::dBetaBinom(ms = SNPdata[i, RefCol], ns = SNPdata[i, RefCol]+SNPdata[i, VarCol],pi = pivec[i], theta = thetavec[i], LOG = TRUE))
        }else{
          LikVec <- c(LikVec, SNPdata[i, WeightCol]*maelstRom::dBetaBinom(ms = SNPdata[i, RefCol], ns = SNPdata[i, RefCol]+SNPdata[i, VarCol],pi = pivec[i], theta = thetavec[i], LOG = TRUE))
        }
      }
    }
    
    return(-sum(LikVec))
    
  }
  
  
  OptObj <- tryCatch( {optim(par = paramvec, fn = AdvancedLik, method = "Nelder-Mead",
                             SNPdata=SNPdata, PiLabs=PiLabs, ThetaLabs=ThetaLabs, RefCol=RefCol, VarCol=VarCol, WeightCol=WeightCol,
                             PiTransFun=PiTransFun, PiTransInv=PiTransInv, ThetaTransFun=ThetaTransFun, ThetaTransInv=ThetaTransInv, 
                             center=center, PiFix=PiFix)},
                      error = function(e) NULL)
  
  if(is.null(OptObj)){
    return(NULL)
  }else if(length(ThetaLabs) == 1){
    return(list(OptObj$value, OptObj$par, mean(SNPdata[,ThetaLabs[1]])))
  }else{
    return(list(OptObj$value, OptObj$par))
  }
  
}
