#' Detects relationships between beta-binomial parameters (pi, theta) and independent variables
#'
#' \code{AdvancedFitter} attempts to detect relationships between beta-binomial parameters (pi, theta) and any independent variable(s) suspected to play a role,
#' e.g. increased theta (overdispersion) in function of age.
#' It is currently, however, more of an exploratory function, as it does so by fitting these relationships assuming the supplied "WeightCol"
#' reflects the "true" probability that a certain observation belongs to the genotype under study (heterozygotes) as this is kept invariant;
#' there is of course no such thing (in reality, an observation has a one true underlying genotype, not "a probability) and, ideally,
#' the parameter-variable(s) relationships here studied are taken into account during the complete beta-binomial mixture model EM fit.
#' 
#' @param ThetaForm An object describing the (linear) relationship between theta and any independent variables present in `SNPdata`, like `Theta ~ Var1 + Var2`;
#' though the "ThetaLink" input (see later) allows for general linear relationships. 
#' @param PiForm An object describing the (linear) relationship between pi and any independent variables present in `SNPdata`, like `Pi ~ Var1 + Var2`;
#' though the "PiLink" input (see later) allows for general linear relationships. 
#' @param SNPdata Dataframe. A dataframe containing reference- and variant allele counts, in columns with names as given by the "RefCol" and "VarCol" inputs.
#' This dataframe also has to contain any independent variables passed to the "ThetaForm" and "PiForm" input arguments.
#' @param RefCol String. Name of the column in SNPdata containing reference allele counts.
#' @param VarCol String. Name of the column in SNPdata containing variant allele counts.
#' @param WeightCol String. Optional name of a column of SNPdata containinig per-sample weights, that are - if specified - are used in a weighted maximum likelihood fit
#' (maximizing sum(sample-weights * sample-log-likelihoods))
#' @param Pi_start Number. Starting pi value for numerical optimization (when fitting PiForm, this starting value is used as the intercept and all independent variables
#' start as having a regression coefficient of zero)
#' @param Theta_start Number. Starting theta value for numerical optimization (when fitting ThetaForm, this starting value is used as the intercept and all independent variables
#' start as having a regression coefficient of zero)
#' @param PiLink String. One of "identity", "log" or "sqrt"; the linear relationship specified by PiForm is fit as PiLink(pi) ~ PiForm
#' @param ThetaLink String. One of "identity", "log" or "sqrt"; the linear relationship specified by ThetaForm is fit as ThetaLink(theta) ~ ThetaForm
#' @param center Logical. If TRUE, centers all exploratory variables (relative to their mean) before performing the maximum likelihood fit.
#' This is recommended when the supplied Pi_start is an expected/mean pi-value across all samples (e.g. the result of a single pi-fit on these samples)
#' @param ResetThetaMin Number. When the supplied Theta_start value is lower than this input, it is reset to this input. Default 10^-10;
#' it is not recommended to change this value.
#' @param ResetThetaMax Number. When the supplied Theta_start value is higher than this input, it is reset to this input. Default 10^-1;
#' it is not recommended to change this value.
#' @export
#' @return A list containing the following components:
#' \item{first list-item}{The negative of the maximized log-likelihood.}
#' \item{second list-item}{Optimized parameter values, given in the order (1) pi-intercept, (2) theta-intercept, (3) pi regression coefficients, (4) theta regression coefficients.}
#' \item{optional third list-item}{If theta depends on only one independent variable, the mean of that variable across samples.}

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
