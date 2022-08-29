#' Returns parameter- and likelihood-distances based on leave-on-out (case-deletion) MLE's of \code{EMfit_betabinom_robust}'s fitted model
#'
#' \code{LikelyDistsHet} is mainly for internal use by \code{EMfit_betabinom_robust}. It assists in robustifying the function's EM-fit by iteratively re-fitting
#' the model on the entire input dataset except for one point, after which the difference in heterozygous pi- and theta-estimates and likelihoods is logged.
#' The difference with their respective full-data fit counterparts is a measure for the left-out data point's influence on the model fit;
#' if either one is sufficiently high the data point could be considered an outlier.
#'
#' @param ref_counts Numeric vector. reference counts.
#' @param var_counts Numeric vector. variant counts.
#' @param sprv Numeric vector. Each sample's EM-weight reflecting its likelihood to be part of the heterozygous population.
#' @param parvec_cur Numeric vector. Pi and theta (in that order) of the heterozyous peak of the full-data fit.
#' @param NoSplitHet Logical. If TRUE, don't allow the beta-binomial fit for heterozygotes to be bimodal
#' @param ResetThetaMin Number. Initial theta values in numeric optimization get capped at this minimum (e.g. in case the moment estimate is even lower)
#' @param ResetThetaMax Number. Initial theta values in numeric optimization get capped at this maximum (e.g. in case the moment estimate is even higher)
#' @param SE Number. Sequencing error rate.
#' @param ReEstPars Logical. If TRUE, re-estimates \code{parvec_cur} given \code{ref_counts} and \code{var_counts}. This is useless if these are the actual
#' counts of the full dataset, but are useful for an emperical approach in which "expected" parameter- and likelihood-distances if the assumed model is 100%
#' correct are simulated by drawing \code{ref_counts} and \code{var_counts} from this assumed model (see \code{EMfit_betabinom_robust})
#' @export
#' @return A list containing the following components:
#' \item{LikDists}{A vector containing likelihood distances per sample (2 times full-data log-likelihood minus re-fitted log-likelihood leaving out the sample).}
#' \item{PiDists}{A vector containing pi distances per sample (leave-sample-out refitted pi minus full-data pi).}
#' \item{ThetaDists}{A vector containing theta distances per sample (leave-sample-out refitted theta minus full-data theta).}


LikelyDistsHet <- function(ref_counts, var_counts, sprv, parvec_cur, NoSplitHet, ResetThetaMin, ResetThetaMax, SE, ReEstPars = FALSE){
  
  LogLikComp_het <- function(parvec, ref_counts, var_counts, sprv){
    return(-sum(sprv*(dBetaBinom(ref_counts, ref_counts+var_counts, pi=gtools::inv.logit(parvec[1], max = 1-10^-16, min = 10^-16), theta=exp(parvec[2]), LOG = TRUE))))
  }
  GradComp_het <- function(parvec, ref_counts, var_counts, sprv){
    Grad1 <- sum(sprv*grad_pi(ref_counts, ref_counts+var_counts, pi=gtools::inv.logit(parvec[1], max = 1-10^-16, min = 10^-16), theta=exp(parvec[2])) * (exp(parvec[1])/(1+exp(parvec[1]))^2) )
    Grad2 <- sum(sprv*grad_theta(ref_counts, ref_counts+var_counts, pi=gtools::inv.logit(parvec[1], max = 1-10^-16, min = 10^-16), theta=exp(parvec[2])) * exp(parvec[2]) )
    return(c(-Grad1, -Grad2))
  }
  eval_g_f_2 <- function(parvec, ref_counts, var_counts, sprv) { # Specifies the inequality constraint between parameters in case NoSplitHet == TRUE
    return(( max(1 - gtools::inv.logit(parvec[1], max = 1-10^-16, min = 10^-16), gtools::inv.logit(parvec[1], max = 1-10^-16, min = 10^-16)) ) - exp(parvec[2]))
  }
  
  # Current log-likelihood
  if(ReEstPars){
    
    probshift <- parvec_cur[1]
    theta_het <- parvec_cur[2]
    if(NoSplitHet){
      theta_het_clone <- theta_het
      OptObj <- tryCatch( {optim(par = c(gtools::logit(probshift), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, gr = GradComp_het, method = "BFGS",
                                 ref_counts = ref_counts, var_counts=var_counts, sprv = sprv)},
                          error = function(e) NULL)
      if(is.null(OptObj)){
        OptObj <- optim(par = c(gtools::logit(probshift), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, method = "BFGS",
                        ref_counts = ref_counts, var_counts=var_counts, sprv = sprv)
      }
      
      probshift_i <- gtools::inv.logit(OptObj$par[1], max = 1-10^-16, min = 10^-16)
      theta_het_i <- exp(OptObj$par[2])
      
      if(theta_het_i  > max(probshift_i, 1-probshift_i)){
        
        theta_het_i <- min(c(theta_het_clone, max((1-probshift_i)/10, probshift_i/10)) ) # Take some distance from the boundary at which bimodality occurs
        
        AObj <- tryCatch( {alabama::auglag(c(gtools::logit(probshift_i), log(min(max(theta_het_i, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, gr = GradComp_het, hin = eval_g_f_2,
                                           ref_counts = ref_counts, var_counts=var_counts, sprv = sprv, control.outer = list("trace" = FALSE))},
                          error = function(e) NULL)
        if(is.null(AObj)){
          AObj <- alabama::auglag(c(gtools::logit(probshift_i), log(min(max(theta_het_i, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, hin = eval_g_f_2,
                                  ref_counts = ref_counts, var_counts=var_counts, sprv = sprv, control.outer = list("trace" = FALSE))
        }
        
        probshift_i <- gtools::inv.logit(AObj$par[1], max = 1-10^-16, min = 10^-16)
        theta_het_i <- exp(AObj$par[2])
        
      }
    } else{
      OptObj <- tryCatch( {optim(par = c(gtools::logit(probshift), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, gr = GradComp_het, method = "BFGS",
                                 ref_counts = ref_counts, var_counts=var_counts, sprv = sprv)},
                          error = function(e) NULL)
      if(is.null(OptObj)){
        OptObj <- optim(par = c(gtools::logit(probshift), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, method = "BFGS",
                        ref_counts = ref_counts, var_counts=var_counts, sprv = sprv)
      }
      
      probshift_i <- gtools::inv.logit(OptObj$par[1], max = 1-10^-16, min = 10^-16)
      theta_het_i <- exp(OptObj$par[2])
    }
    
    if(probshift_i > (1-10*SE)){
      probshift_i <- (1-10*SE)
    } else if(probshift_i < 10*SE){
      probshift_i <- 10*SE
    }
    
    probshift <- probshift_i
    theta_het <- theta_het_i
    CurLL <- LogLikComp_het(c(probshift, theta_het), ref_counts, var_counts, sprv)
    
  }else{
    CurLL <- LogLikComp_het(parvec_cur, ref_counts, var_counts, sprv)
    probshift <- parvec_cur[1]
    theta_het <- parvec_cur[2]
  }
  
  
  # Now we need to, one-by-one, drop each datapoint and redo the parameter calculation
  NewLLs <- c()
  PiDists <- c()
  ThetaDists <- c()
  for(i in 1:length(ref_counts)){
    ref_counts_i <- ref_counts[-i]
    var_counts_i <- var_counts[-i]
    sprv_i <- sprv[-i]
    
    
    
    if(NoSplitHet){
      theta_het_clone <- theta_het
      OptObj <- tryCatch( {optim(par = c(gtools::logit(probshift), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, gr = GradComp_het, method = "BFGS",
                                 ref_counts = ref_counts_i, var_counts=var_counts_i, sprv = sprv_i)},
                          error = function(e) NULL)
      if(is.null(OptObj)){
        OptObj <- optim(par = c(gtools::logit(probshift), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, method = "BFGS",
                        ref_counts = ref_counts_i, var_counts=var_counts_i, sprv = sprv_i)
      }
      
      probshift_i <- gtools::inv.logit(OptObj$par[1], max = 1-10^-16, min = 10^-16)
      theta_het_i <- exp(OptObj$par[2])
      
      if(theta_het_i  > max(probshift_i, 1-probshift_i)){
        
        theta_het_i <- min(c(theta_het_clone, max((1-probshift_i)/10, probshift_i/10)) ) # Take some distance from the boundary at which bimodality occurs
        
        AObj <- tryCatch( {alabama::auglag(c(gtools::logit(probshift_i), log(min(max(theta_het_i, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, gr = GradComp_het, hin = eval_g_f_2,
                                           ref_counts = ref_counts_i, var_counts=var_counts_i, sprv = sprv_i, control.outer = list("trace" = FALSE))},
                          error = function(e) NULL)
        if(is.null(AObj)){
          AObj <- alabama::auglag(c(gtools::logit(probshift_i), log(min(max(theta_het_i, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, hin = eval_g_f_2,
                                  ref_counts = ref_counts_i, var_counts=var_counts_i, sprv = sprv_i, control.outer = list("trace" = FALSE))
        }
        
        probshift_i <- gtools::inv.logit(AObj$par[1], max = 1-10^-16, min = 10^-16)
        theta_het_i <- exp(AObj$par[2])
        
      }
    } else{
      OptObj <- tryCatch( {optim(par = c(gtools::logit(probshift), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, gr = GradComp_het, method = "BFGS",
                                 ref_counts = ref_counts_i, var_counts=var_counts_i, sprv = sprv_i)},
                          error = function(e) NULL)
      if(is.null(OptObj)){
        OptObj <- optim(par = c(gtools::logit(probshift), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, method = "BFGS",
                        ref_counts = ref_counts_i, var_counts=var_counts_i, sprv = sprv_i)
      }
      
      probshift_i <- gtools::inv.logit(OptObj$par[1], max = 1-10^-16, min = 10^-16)
      theta_het_i <- exp(OptObj$par[2])
    }
    
    if(probshift_i > (1-10*SE)){
      probshift_i <- (1-10*SE)
    } else if(probshift_i < 10*SE){
      probshift_i <- 10*SE
    }
    
    
    LLi <- LogLikComp_het(c(probshift_i, theta_het_i), ref_counts, var_counts, sprv)
    NewLLs <- c(NewLLs, LLi)
    PiDists <- c(PiDists, probshift_i - probshift)
    ThetaDists <- c(ThetaDists, theta_het_i - theta_het)
    
  }
  
  out <- 2*(CurLL - NewLLs)
  return(list("LikDists" = out, "PiDists" = PiDists, "ThetaDists" = ThetaDists))
  
}
