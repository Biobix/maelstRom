#' Models and estimates the allelic shift for cis-eQTL data, using beta-binomial models.
#'
#' \code{EMfit_betabinom_popcomb}, similar to \code{EMfit_betabinom}, fits a beta-binomial mixture model to input reference- and allele RNAseq counts,
#' but with those counts originating from two different populations (a control- and case-population). All distribution models are jointly fitted except
#' for the heterozygous overdispersion parameters, which is allowed to differ between the populations. The main purpose of this fit is getting this model's
#' likelihood to be used in a likelihood ratio test, checking for statistical evidence that this heterozygous overdispersion parameter 
#' truly differs between the populations, which would indicate differential allelic divergence.
#'
#' @param data_counts Data frame. Data frame of a SNP with reference and variant counts ("ref_count" and
#'     "var_count", respectively) for each sample ("sample"). Also required is the "isCase" column, equaling 1 for case-samples and 0 for control-samples
#' @param allelefreq Number. Allele frequency. Only used when pInitEst is TRUE (default = 0.5)
#' @param SE Number. Sequencing error rate.
#' @param inbr Number. Degree of inbreeding (default = 0).
#' @param dltaco Number. Minimal difference between 2 iterations (default = 0.001).
#' @param HWE Logical. Should HWE be used for allele frequency estimation, not recommended (default = FALSE).
#' @param p_InitEst Logical. Calculate initial estimates of pr, pv and prv from allelefreq, not recommended (default = FALSE).
#' @param ThetaInits Numeric vector, "moment" or "TryThree". A vector of length two containing the initial dispersion estimate for the homozygous peaks, followed by the one for the heterozygous peak,
#'    or "moment", in which case moment initial estimates are generated. In case of "TryThree", three initial values given by the \code{thetaTRY} input argument are tried out and the one yielding the best fit is retained.
#' @param ReEstThetas String. Accepts the methods "moment" or "simple" as arguments, in which case dispersions are re-estimated during every EM-step using a moment estimate or a simple custom procedure respectively, to be used as a starting
#   point for numerical optimization. If set to anything else, the dispersion estimates resulting from the previous EM iterations are used.
#' @param NoSplitHom Logical. If TRUE, don't allow the beta-binomial fits for homozygotes to be bimodal
#' @param NoSplitHet Logical. If TRUE, don't allow the beta-binomial fit for heterozygotes to be bimodal
#' @param ResetThetaMin Number. Initial theta values in numeric optimization get capped at this minimum (e.g. in case the moment estimate is even lower)
#' @param ResetThetaMax Number. Initial theta values in numeric optimization get capped at this maximum (e.g. in case the moment estimate is even higher)
#' @param thetaTRY Numeric vector. See \code{ReEstThetas}, in case this argument equals "TryThree".
#' @param probshift_init Number. Initial estimate position heterozygous peak (allelic bias; expected fraction of reference reads)
#' @export
#' @return A list containing the following components:
#' \item{ParamVec}{A vector containing fitted model parameters: genotype frequencies, heterozygous peak position (allelic bias), and overdispersion of the homozygous peaks and of the separate heterozygous control- and case-peaks.}
#' \item{GenoDF}{A dataframe containing per-samples genotypes likelihoods according to the fitted model.}
#' \item{genotypeN}{A character vector containing the most likely per-sample genotypes according to the fitted model.}
#' \item{nrep}{The number of iterations.}
#' \item{quality}{An "!" indicates a bad quality locus or fit (no apparent heterozygotes)}

EMfit_betabinom_popcomb <- function(data_counts, allelefreq=0.5, SE, inbr = 0, dltaco = 10^-6, HWE = FALSE, p_InitEst = FALSE, 
                               ThetaInits = "moment", ReEstThetas = "moment", NoSplitHom = TRUE, NoSplitHet = TRUE,
                               ResetThetaMin = 10^-10, ResetThetaMax = 10^-1, thetaTRY = c(10^-1, 10^-3, 10^-7),
                               probshift_init = 0.5) {
  
  LogLikComp_hom <- function(theta, SE, ref_counts, var_counts, spr, spv){
    return(-sum(spr*(dBetaBinom(ref_counts, ref_counts+var_counts, pi=1-SE, theta=exp(theta), LOG = TRUE)) +
                  spv*(dBetaBinom(var_counts,ref_counts+var_counts, pi=1-SE, theta=exp(theta), LOG = TRUE))))
  }
  GradComp_hom <- function(theta, SE, ref_counts, var_counts, spr, spv){
    Grad <- sum(spr*grad_theta(ref_counts, ref_counts+var_counts, pi=1-SE, theta=exp(theta)) * exp(theta) +
                  spv*grad_theta(var_counts, ref_counts+var_counts, pi=1-SE, theta=exp(theta)) * exp(theta))
    return(-Grad)
  }
  LogLikComp_het <- function(parvec, ref_counts, var_counts, sprv){
    return(-sum(sprv*(dBetaBinom(ref_counts, ref_counts+var_counts, pi=gtools::inv.logit(parvec[1], max = 1-10^-16, min = 10^-16), theta=exp(parvec[2]), LOG = TRUE))))
  }
  GradComp_het <- function(parvec, ref_counts, var_counts, sprv){
    Grad1 <- sum(sprv*grad_pi(ref_counts, ref_counts+var_counts, pi=gtools::inv.logit(parvec[1], max = 1-10^-16, min = 10^-16), theta=exp(parvec[2])) * (exp(parvec[1])/(1+exp(parvec[1]))^2) )
    Grad2 <- sum(sprv*grad_theta(ref_counts, ref_counts+var_counts, pi=gtools::inv.logit(parvec[1], max = 1-10^-16, min = 10^-16), theta=exp(parvec[2])) * exp(parvec[2]) )
    return(c(-Grad1, -Grad2))
  }
  
  LogLikComp_het_H0 <- function(theta, probshift, ref_counts, var_counts, sprv){
    return(-sum(sprv*(dBetaBinom(ref_counts, ref_counts+var_counts, pi=probshift, theta=exp(theta), LOG = TRUE))))
  }
  GradComp_het_H0 <- function(theta, probshift, ref_counts, var_counts, sprv){
    return(-sum(sprv*grad_theta(ref_counts, ref_counts+var_counts, pi=probshift, theta=exp(theta)) * exp(theta) ))
  }
  
  eval_g_f <- function(theta, SE, ref_counts, var_counts, spr, spv) { # Specifies the inequality constraint between parameters in case NoSplitHom == TRUE
    return(( max(1 - SE, SE) ) - exp(theta))
  }
  
  eval_g_f_2 <- function(parvec, ref_counts, var_counts, sprv) { # Specifies the inequality constraint between parameters in case NoSplitHet == TRUE
    return(( max(1 - gtools::inv.logit(parvec[1], max = 1-10^-16, min = 10^-16), gtools::inv.logit(parvec[1], max = 1-10^-16, min = 10^-16)) ) - exp(parvec[2]))
  }
  
  eval_g_f_2_H0 <- function(theta, probshift, ref_counts, var_counts, sprv) { # Specifies the inequality constraint between parameters in case NoSplitHet == TRUE
    return(( max(1 - probshift, probshift) ) - exp(theta))
  }
  
  LogLikComp_het_CnT <- function(parvec, ref_counts, var_counts, isCase, sprv){
    return(-sum(sprv*(  c(dBetaBinom(ref_counts[isCase==0], (ref_counts+var_counts)[isCase==0], pi=gtools::inv.logit(parvec[1], max = 1-10^-16, min = 10^-16), theta=exp(parvec[2]), LOG = TRUE),
                          dBetaBinom(ref_counts[isCase==1], (ref_counts+var_counts)[isCase==1], pi=gtools::inv.logit(parvec[1], max = 1-10^-16, min = 10^-16), theta=exp(parvec[3]), LOG = TRUE))  )))
  }
  GradComp_het_CnT <- function(parvec, ref_counts, var_counts, isCase, sprv){
    Grad1 <- sum(sprv*c(grad_pi(ref_counts[isCase==0], (ref_counts+var_counts)[isCase==0], pi=gtools::inv.logit(parvec[1], max = 1-10^-16, min = 10^-16), theta=exp(parvec[2])),
                        grad_pi(ref_counts[isCase==1], (ref_counts+var_counts)[isCase==1], pi=gtools::inv.logit(parvec[1], max = 1-10^-16, min = 10^-16), theta=exp(parvec[3])))* (exp(parvec[1])/(1+exp(parvec[1]))^2) )
    Grad2 <- sum(sprv[isCase==0]*grad_theta(ref_counts[isCase==0], (ref_counts+var_counts)[isCase==0], pi=gtools::inv.logit(parvec[1], max = 1-10^-16, min = 10^-16), theta=exp(parvec[2])) * exp(parvec[2]) )
    Grad3 <- sum(sprv[isCase==1]*grad_theta(ref_counts[isCase==1], (ref_counts+var_counts)[isCase==1], pi=gtools::inv.logit(parvec[1], max = 1-10^-16, min = 10^-16), theta=exp(parvec[3])) * exp(parvec[3]) )
    return(c(-Grad1, -Grad2, -Grad3))
  }
  eval_g_f_2_CnT <- function(parvec, ref_counts, var_counts, isCase, sprv) { # Specifies the inequality constraint between parameters in case NoSplitHet == TRUE
    return(  c(((max(1 - gtools::inv.logit(parvec[1], max = 1-10^-16, min = 10^-16), gtools::inv.logit(parvec[1], max = 1-10^-16, min = 10^-16)) ) - exp(parvec[2])),
             ((max(1 - gtools::inv.logit(parvec[1], max = 1-10^-16, min = 10^-16), gtools::inv.logit(parvec[1], max = 1-10^-16, min = 10^-16)) ) - exp(parvec[2]))) )
  }
  
  data_counts <- rbind(data_counts[data_counts$isCase==0,], data_counts[data_counts$isCase==1,])
  
  probshift <- probshift_init # Initial estimate position heterozygous peak
  dlta <- 1
  nrep <- 0 
  if(p_InitEst){
    prv <- 2 * allelefreq * (1 - allelefreq) * (1 - inbr) 
    pr <- allelefreq^2 + inbr * allelefreq * (1 - allelefreq)
    pv <- (1 - allelefreq)^2 + inbr * allelefreq * (1 - allelefreq) 
  } else{
    pr <- 1/3
    prv <- 1/3
    pv <- 1/3
  }
  
  pr_H0 <- pr
  pv_H0 <- pv
  prv_H0 <- prv
  
  if(is.null(ThetaInits)){
    theta_hom <- ResetThetaMin
    theta_het_control <- ResetThetaMin
    theta_het_tumor <- ResetThetaMin
  }else if(ThetaInits == "moment"){
    
    # First fit some regular binomial models to get an initial categorisation of the data points
    spr_est <- pr * dbinom(data_counts$ref_count, data_counts$ref_count + data_counts$var_count, prob = 1 - SE)
    spv_est <- pv * dbinom(data_counts$var_count,data_counts$ref_count + data_counts$var_count, prob = 1 - SE)
    sprv_est <- prv * dbinom(data_counts$ref_count, data_counts$ref_count + data_counts$var_count, prob = 0.5)
    
    pdata <- rowSums(cbind(spv_est, sprv_est, spr_est)) 
    if (any(pdata==0)){
      ProblemCases <- which(pdata==0)
      for(case in ProblemCases){
        var_part<-data_counts$var_count[case]
        ref_part<-data_counts$ref_count[case]
        spv_part<-dbinom(var_part, ref_part+var_part, 1-SE, log = TRUE)
        sprv_part<-dbinom(ref_part, ref_part+var_part, 0.5, log = TRUE)
        spr_part<-dbinom(ref_part, ref_part+var_part, 1-SE, log = TRUE)
        spvec<-c(spr_part, sprv_part, spv_part)
        if(spr_part==max(spvec)){
          spr_est[case]<-1
          pdata[case]<-1
        } else if(sprv_part==max(spvec)){
          sprv_est[case]<-1
          pdata[case]<-1
        } else{
          spv_est[case]<-1
          pdata[case]<-1
        }
      }
    }
    spr_est <- spr_est/pdata
    spv_est <- spv_est/pdata
    sprv_est <- sprv_est/pdata
    
    MomentEst_complete <- MomentEst_MixedBetaBinom(ref_counts = data_counts$ref_count, var_counts = data_counts$var_count, spr = spr_est, 
                                                   spv = spv_est, sprv = sprv_est, pi_hom_fix = 1-SE, pi_het_fix = 0.5)
    MomentEst_control <- MomentEst_MixedBetaBinom(ref_counts = data_counts$ref_count[data_counts$isCase==0], 
                                                  var_counts = data_counts$var_count[data_counts$isCase==0], spr = spr_est[data_counts$isCase==0], 
                                                  spv = spv_est[data_counts$isCase==0], sprv = sprv_est[data_counts$isCase==0], pi_hom_fix = 1-SE, pi_het_fix = 0.5)
    MomentEst_tumor <- MomentEst_MixedBetaBinom(ref_counts = data_counts$ref_count[data_counts$isCase==1], 
                                                var_counts = data_counts$var_count[data_counts$isCase==1], spr = spr_est[data_counts$isCase==1], 
                                                spv = spv_est[data_counts$isCase==1], sprv = sprv_est[data_counts$isCase==1], pi_hom_fix = 1-SE, pi_het_fix = 0.5)
    
    theta_hom <- min(max(MomentEst_complete["theta_hom"], ResetThetaMin), ResetThetaMax)
    theta_het_control <- min(max(MomentEst_control["theta_het"], ResetThetaMin), ResetThetaMax)
    theta_het_tumor <- min(max(MomentEst_tumor["theta_het"], ResetThetaMin), ResetThetaMax)
    if(is.nan(theta_het_control) | is.infinite(theta_het_control)){
      theta_het_control <- ResetThetaMin
    }
    if(is.nan(theta_het_tumor) | is.infinite(theta_het_tumor)){
      theta_het_tumor <- ResetThetaMin
    }
    if(is.nan(theta_hom) | is.infinite(theta_hom)){
      theta_hom <- ResetThetaMin
    }
    theta_hom_init <- theta_hom
    theta_het_control_init <- theta_het_control
    theta_het_tumor_init <- theta_het_tumor
  } else{
    theta_hom <- min(max(ThetaInits[1], ResetThetaMin), ResetThetaMax)
    theta_het_control <- min(max(ThetaInits[2], ResetThetaMin), ResetThetaMax)
    theta_het_tumor <- min(max(ThetaInits[3], ResetThetaMin), ResetThetaMax)
  }
  
  theta_hom <- min(c(theta_hom, max((1-SE)/10, SE/10)) ) # Take some distance from the boundary at which bimodality occurs
  theta_het_control <- min(c(theta_het_control, max((1-probshift)/10, probshift/10)) ) # Take some distance from the boundary at which bimodality occurs
  theta_het_tumor <- min(c(theta_het_tumor, max((1-probshift)/10, probshift/10)) ) 
  
  
  quality <- ""
  
  ########################################################################################################################################################################
  ########################################################################################################################################################################
  
  Q <- 1000
  
  while (dlta > dltaco & nrep < 100) {
    Qold <- Q
    nrep <- nrep + 1
    
    spr <- pr * dBetaBinom(data_counts$ref_count, data_counts$ref_count + data_counts$var_count, pi = 1 - SE, theta = theta_hom, LOG = FALSE)
    spv <- pv * dBetaBinom(data_counts$var_count, data_counts$ref_count + data_counts$var_count, pi = 1 - SE, theta = theta_hom, LOG = FALSE)
    sprv <- prv * c(dBetaBinom(data_counts$ref_count[data_counts$isCase==0], (data_counts$ref_count + data_counts$var_count)[data_counts$isCase==0], pi = probshift, theta = theta_het_control, LOG = FALSE), 
                    dBetaBinom(data_counts$ref_count[data_counts$isCase==1], (data_counts$ref_count + data_counts$var_count)[data_counts$isCase==1], pi = probshift, theta = theta_het_tumor, LOG = FALSE))
    
    pdata <- rowSums(cbind(spv, sprv, spr)) 
    
    if (any(pdata==0)){ 
      ProblemCases <- which(pdata==0)
      for(case in ProblemCases){
        var_part<-data_counts$var_count[case]
        ref_part<-data_counts$ref_count[case]
        tumor_part <- data_counts$isCase[case]
        spr_part<-dBetaBinom(ref_part, ref_part+var_part, pi = 1-SE, theta = theta_hom, LOG = TRUE)
        spv_part<-dBetaBinom(var_part, ref_part+var_part, pi = 1-SE, theta = theta_hom, LOG = TRUE)
        
        if(tumor_part==0){
          sprv_part<-dBetaBinom(ref_part, ref_part+var_part, pi = probshift, theta = theta_het_control, LOG = TRUE)
        }else{
          sprv_part<-dBetaBinom(ref_part, ref_part+var_part, pi = probshift, theta = theta_het_tumor, LOG = TRUE)
        }
        
        spvec<-c(spr_part, sprv_part, spv_part)
        if(spr_part==max(spvec)){
          spr[case]<-1
          pdata[case]<-1
        } else if(sprv_part==max(spvec)){
          sprv[case]<-1
          pdata[case]<-1
        } else{
          spv[case]<-1
          pdata[case]<-1
        }
      }
    }
    
    spv <- spv/pdata
    sprv <- sprv/pdata
    spr <- spr/pdata
    
    if (sum(sprv)==0) {
      quality <- "!"
      break
    } 
    
    allelefreq <- mean(spr) + mean(sprv)/2 
    
    if (HWE) {
      prv <- 2 * allelefreq * (1 - allelefreq) * (1 - inbr)
      pr <- allelefreq^2 + inbr * allelefreq * (1 - allelefreq)
      pv <- (1 - allelefreq)^2 + inbr * allelefreq * (1 - allelefreq)
    }else {
      pv <- mean(spv) 
      prv <- mean(sprv)
      pr <- mean(spr)
    }
    
    if(ReEstThetas == "moment"){
      MomentEst_complete <- MomentEst_MixedBetaBinom(ref_counts = data_counts$ref_count, var_counts = data_counts$var_count, spr = spr, spv = spv, sprv = sprv, 
                                            pi_hom_fix = 1-SE, pi_het_fix = probshift)
      MomentEst_control <- MomentEst_MixedBetaBinom(ref_counts = data_counts$ref_count[data_counts$isCase==0], var_counts = data_counts$var_count[data_counts$isCase==0],
                                                    spr = spr[data_counts$isCase==0], spv = spv[data_counts$isCase==0], sprv = sprv[data_counts$isCase==0], 
                                                    pi_hom_fix = 1-SE, pi_het_fix = probshift)
      MomentEst_tumor <- MomentEst_MixedBetaBinom(ref_counts = data_counts$ref_count[data_counts$isCase==1], var_counts = data_counts$var_count[data_counts$isCase==1],
                                                  spr = spr[data_counts$isCase==1], spv = spv[data_counts$isCase==1], sprv = sprv[data_counts$isCase==1], 
                                                  pi_hom_fix = 1-SE, pi_het_fix = probshift)
      if(sum(spr+spv)!=0){
        theta_hom <- MomentEst_complete["theta_hom"]
      }
      theta_het_control <- MomentEst_control["theta_het"]
      theta_het_tumor <- MomentEst_tumor["theta_het"]
      if(is.nan(theta_het_control) | is.infinite(theta_het_control) | theta_het_control < ResetThetaMin){
        theta_het_control <- ResetThetaMin
      }
      if(is.nan(theta_het_tumor) | is.infinite(theta_het_tumor) | theta_het_tumor < ResetThetaMin){
        theta_het_tumor <- ResetThetaMin
      }
      if(is.nan(theta_hom) | is.infinite(theta_hom) | theta_hom < ResetThetaMin){
        theta_hom <- ResetThetaMin
      }
      
      if(sum(spr+spv)!=0){
        theta_hom_vec <- c()
        for(ii in c(10^(-5:-1), theta_hom) ){
          Likely <- -LogLikComp_hom(log(ii), SE = SE, ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr, spv = spv)
          theta_hom_vec <- c(theta_hom_vec, Likely)
        }
        WhichOne <- which(theta_hom_vec==max(theta_hom_vec))
        theta_hom <- c(10^(-5:-1), theta_hom)[WhichOne]
      }
      
      theta_het_control_vec <- c()
      for(ii in c(10^(-5:-1), theta_het_control)){
        Likely <- -LogLikComp_het(c(gtools::logit(probshift), log(ii)), ref_counts=data_counts$ref_count[data_counts$isCase==0], var_counts=data_counts$var_count[data_counts$isCase==0], sprv = sprv[data_counts$isCase==0])
        theta_het_control_vec <- c(theta_het_control_vec, Likely)
      }
      WhichOne <- which(theta_het_control_vec==max(theta_het_control_vec))
      theta_het_control <- c(10^(-5:-1), theta_het_control)[WhichOne]
      
      theta_het_tumor_vec <- c()
      for(ii in c(10^(-5:-1), theta_het_tumor)){
        Likely <- -LogLikComp_het(c(gtools::logit(probshift), log(ii)), ref_counts=data_counts$ref_count[data_counts$isCase==1], var_counts=data_counts$var_count[data_counts$isCase==1], sprv = sprv[data_counts$isCase==1])
        theta_het_tumor_vec <- c(theta_het_tumor_vec, Likely)
      }
      WhichOne <- which(theta_het_tumor_vec==max(theta_het_tumor_vec))
      theta_het_tumor <- c(10^(-5:-1), theta_het_tumor)[WhichOne]
      
    } else if(ReEstThetas == "simple" | nrep == 1){
      if(sum(spr+spv)!=0){
        theta_hom_vec <- c()
        for(ii in 10^(-5:-1)){
          Likely <- -LogLikComp_hom(log(ii), SE = SE, ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr, spv = spv)
          theta_hom_vec <- c(theta_hom_vec, Likely)
        }
        WhichOne <- which(theta_hom_vec==max(theta_hom_vec))
        theta_hom <- 10^(-5:-1)[WhichOne]
      }
      
      theta_het_control_vec <- c()
      for(ii in 10^(-5:-1)){
        Likely <- -LogLikComp_het(c(gtools::logit(probshift), log(ii)), ref_counts=data_counts$ref_count[data_counts$isCase==0], var_counts=data_counts$var_count[data_counts$isCase==0], sprv = sprv[data_counts$isCase==0])
        theta_het_control_vec <- c(theta_het_control_vec, Likely)
      }
      WhichOne <- which(theta_het_control_vec==max(theta_het_control_vec))
      theta_het_control <- 10^(-5:-1)[WhichOne]
      
      theta_het_tumor_vec <- c()
      for(ii in 10^(-5:-1)){
        Likely <- -LogLikComp_het(c(gtools::logit(probshift), log(ii)), ref_counts=data_counts$ref_count[data_counts$isCase==1], var_counts=data_counts$var_count[data_counts$isCase==1], sprv = sprv[data_counts$isCase==1])
        theta_het_tumor_vec <- c(theta_het_tumor_vec, Likely)
      }
      WhichOne <- which(theta_het_tumor_vec==max(theta_het_tumor_vec))
      theta_het_tumor <- 10^(-5:-1)[WhichOne]
    }
    
    #compareDerivatives(LogLikComp_hom, GradComp_hom, t0 = log(0.01), SE = 0.002, ref_counts=ref_counts, var_counts=var_counts, spr = spr, spv = spv)
    
    if(ReEstThetas != "TryThree"){
      
      if(sum(spr+spv)!=0){
        if(NoSplitHom){
          theta_hom_clone <- theta_hom
          OptObj <- optim(par = log(min(max(theta_hom, ResetThetaMin), ResetThetaMax)), fn = LogLikComp_hom, gr = GradComp_hom, method = "BFGS", SE = SE, 
                          ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr, spv = spv)
          theta_hom <- exp(OptObj$par)
          
          if(theta_hom > max(SE, 1-SE)){
            
            theta_hom <- min(c(theta_hom_clone, max((1-SE)/10, SE/10)) ) # Take some distance from the boundary at which bimodality occurs
            
            AObj <- alabama::auglag(log(min(max(theta_hom, ResetThetaMin), ResetThetaMax)), fn = LogLikComp_hom, gr = GradComp_hom, hin = eval_g_f, SE = SE,
                                    ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr, spv = spv, control.outer = list("trace" = FALSE))
            
            theta_hom <- exp(AObj$par[1])
            
          }
          
          
        } else{
          OptObj <- optim(par = log(min(max(theta_hom, ResetThetaMin), ResetThetaMax)), fn = LogLikComp_hom, gr = GradComp_hom, method = "BFGS", SE = SE, 
                          ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr, spv = spv)
          theta_hom <- exp(OptObj$par)
        }
      }
      
      #compareDerivatives(LogLikComp_hom, GradComp_hom, HessComp_hom, t0 = 0.3, SE=SE, ref_counts=ref_counts, var_counts=var_counts, spr=spr, spv=spv)
      #compareDerivatives(LogLikComp_het, GradComp_het, t0 = c(gtools::logit(1-0.002), log(0.5)), ref_counts=ref_counts, var_counts=var_counts, sprv = sprv)
      
      if(NoSplitHet){
        theta_het_control_clone <- theta_het_control
        theta_het_tumor_clone <- theta_het_tumor
        OptObj <- tryCatch( {optim(par = c(gtools::logit(probshift), log(min(max(theta_het_control, ResetThetaMin), ResetThetaMax)), log(min(max(theta_het_tumor, ResetThetaMin), ResetThetaMax))), 
                                   fn = LogLikComp_het_CnT, gr = GradComp_het_CnT, method = "BFGS",
                                   ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, isCase = data_counts$isCase)},
                            error = function(e) NULL)
        if(is.null(OptObj)){
          OptObj <- optim(par = c(gtools::logit(probshift), log(min(max(theta_het_control, ResetThetaMin), ResetThetaMax)), log(min(max(theta_het_tumor, ResetThetaMin), ResetThetaMax))), 
                          fn = LogLikComp_het_CnT, method = "BFGS",
                          ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, isCase = data_counts$isCase)
        }
        
        probshift <- gtools::inv.logit(OptObj$par[1], max = 1-10^-16, min = 10^-16)
        theta_het_control <- exp(OptObj$par[2])
        theta_het_tumor <- exp(OptObj$par[3])
        
        if(theta_het_control  > max(probshift, 1-probshift) | theta_het_tumor  > max(probshift, 1-probshift)){
          
          theta_het_control <- min(c(theta_het_control_clone, max((1-probshift)/10, probshift/10)) ) # Take some distance from the boundary at which bimodality occurs
          theta_het_tumor <- min(c(theta_het_tumor_clone, max((1-probshift)/10, probshift/10)) )
          
          AObj <- tryCatch( {alabama::auglag(c(gtools::logit(probshift), log(min(max(theta_het_control, ResetThetaMin), ResetThetaMax)), log(min(max(theta_het_tumor, ResetThetaMin), ResetThetaMax))),
                                             fn = LogLikComp_het_CnT, gr = GradComp_het_CnT, hin = eval_g_f_2_CnT,
                                             ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, isCase = data_counts$isCase, control.outer = list("trace" = FALSE))},
                            error = function(e) NULL)
          if(is.null(AObj)){
            AObj <- alabama::auglag(c(gtools::logit(probshift), log(min(max(theta_het_control, ResetThetaMin), ResetThetaMax)), log(min(max(theta_het_tumor, ResetThetaMin), ResetThetaMax))), 
                                    fn = LogLikComp_het_CnT, hin = eval_g_f_2_CnT,
                                    ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, isCase = data_counts$isCase, control.outer = list("trace" = FALSE))
          }
          
          probshift <- gtools::inv.logit(AObj$par[1], max = 1-10^-16, min = 10^-16)
          theta_het_control <- exp(AObj$par[2])
          theta_het_tumor <- exp(AObj$par[3])
          
        }
      } else{
        OptObj <- tryCatch( {optim(par = c(gtools::logit(probshift), log(min(max(theta_het_control, ResetThetaMin), ResetThetaMax)), log(min(max(theta_het_tumor, ResetThetaMin), ResetThetaMax))), 
                                   fn = LogLikComp_het_CnT, gr = GradComp_het_CnT, method = "BFGS",
                                   ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, isCase = data_counts$isCase)},
                            error = function(e) NULL)
        if(is.null(OptObj)){
          OptObj <- optim(par = c(gtools::logit(probshift), log(min(max(theta_het_control, ResetThetaMin), ResetThetaMax)), log(min(max(theta_het_tumor, ResetThetaMin), ResetThetaMax))), 
                          fn = LogLikComp_het_CnT, method = "BFGS",
                          ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, isCase = data_counts$isCase)
        }
        
        probshift <- gtools::inv.logit(OptObj$par[1], max = 1-10^-16, min = 10^-16)
        theta_het_control <- exp(OptObj$par[2])
        theta_het_tumor <- exp(OptObj$par[3])
      }
      
    } else{
      
      if(sum(spr+spv)!=0){
        if(NoSplitHom){
          
          thetavec <- c()
          likvec <- c()
          for(TH in thetaTRY){
            OptObj <- optim(par = log(TH), fn = LogLikComp_hom, gr = GradComp_hom, method = "BFGS", SE = SE, 
                            ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr, spv = spv)
            thetavec <- c(thetavec, exp(OptObj$par))
            likvec <- c(likvec, LogLikComp_hom(OptObj$par, SE=SE, ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, spr=spr, spv=spv))
          }
          theta_hom_loc <- which(likvec == min(likvec))[1]
          theta_hom <- thetavec[theta_hom_loc]
          
          if(theta_hom > max(SE, 1-SE)){
            
            thetavec <- c()
            likvec <- c()
            for(TH in thetaTRY){
              AObj <- alabama::auglag(log(TH), fn = LogLikComp_hom, gr = GradComp_hom, hin = eval_g_f, SE = SE,
                                      ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr, spv = spv, control.outer = list("trace" = FALSE))
              thetavec <- c(thetavec, exp(AObj$par[1]))
              likvec <- c(likvec, LogLikComp_hom(AObj$par[1], SE=SE, ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, spr=spr, spv=spv))
            }
            theta_hom_loc <- which(likvec == min(likvec))[1]
            theta_hom <- thetavec[theta_hom_loc]
          }
          
          
        } else{
          thetavec <- c()
          likvec <- c()
          for(TH in thetaTRY){
            OptObj <- optim(par = log(TH), fn = LogLikComp_hom, gr = GradComp_hom, method = "BFGS", SE = SE, 
                            ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr, spv = spv)
            thetavec <- c(thetavec, exp(OptObj$par))
            likvec <- c(likvec, LogLikComp_hom(OptObj$par, SE=SE, ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, spr=spr, spv=spv))
          }
          theta_hom_loc <- which(likvec == min(likvec))[1]
          theta_hom <- thetavec[theta_hom_loc]
        }
      }
      
      #compareDerivatives(LogLikComp_hom, GradComp_hom, HessComp_hom, t0 = 0.3, SE=SE, ref_counts=ref_counts, var_counts=var_counts, spr=spr, spv=spv)
      #compareDerivatives(LogLikComp_het, GradComp_het, t0 = c(gtools::logit(1-0.002), log(0.5)), ref_counts=ref_counts, var_counts=var_counts, sprv = sprv)
      
      if(NoSplitHet){
        
        thetavec_control <- c()
        thetavec_tumor <- c()
        pivec <- c()
        likvec <- c()
        
        
        
        for(TH1 in thetaTRY){
          for(TH2 in thetaTRY){
          
            OptObj <- tryCatch( {optim(par = c(gtools::logit(probshift), log(TH1), log(TH2)), fn = LogLikComp_het_CnT, gr = GradComp_het_CnT, method = "BFGS",
                                       ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, isCase = data_counts$isCase)},
                                error = function(e) NULL)
            if(is.null(OptObj)){
              OptObj <- optim(par = c(gtools::logit(probshift), log(TH1), log(TH2)), fn = LogLikComp_het_CnT, method = "BFGS",
                              ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, isCase = data_counts$isCase)
            }
            thetavec_control <- c(thetavec_control, exp(OptObj$par[2]))
            thetavec_tumor <- c(thetavec_tumor, exp(OptObj$par[3]))
            pivec <- c(pivec, gtools::inv.logit(OptObj$par[1], max = 1-10^-16, min = 10^-16))
            likvec <- c(likvec, LogLikComp_het_CnT(c(OptObj$par[1], OptObj$par[2], OptObj$par[3]), ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, sprv=sprv, isCase = data_counts$isCase))
          
          }
        }
        
        
        
        theta_het_loc <- which(likvec == min(likvec))[1]
        theta_het_control <- thetavec_control[theta_het_loc]
        theta_het_tumor <- thetavec_tumor[theta_het_loc]
        probshift <- pivec[theta_het_loc]
        
        if((theta_het_control  > max(probshift, 1-probshift)) | (theta_het_tumor  > max(probshift, 1-probshift))){
          
          thetavec_control <- c()
          thetavec_tumor <- c()
          pivec <- c()
          likvec <- c()
          for(TH1 in thetaTRY){
            for(TH2 in thetaTRY){
            
              AObj <- tryCatch( {alabama::auglag(c(gtools::logit(probshift), log(TH1), log(TH2)), fn = LogLikComp_het_CnT, gr = GradComp_het_CnT, hin = eval_g_f_2_CnT,
                                                 ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, isCase = data_counts$isCase, control.outer = list("trace" = FALSE))},
                                error = function(e) NULL)
              if(is.null(AObj)){
                AObj <- alabama::auglag(c(gtools::logit(probshift), log(TH1), log(TH2)), fn = LogLikComp_het_CnT, hin = eval_g_f_2_CnT,
                                        ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, isCase = data_counts$isCase, control.outer = list("trace" = FALSE))
              }
              thetavec_control <- c(thetavec_control, exp(AObj$par[2]))
              thetavec_tumor <- c(thetavec_tumor, exp(AObj$par[3]))
              pivec <- c(pivec, gtools::inv.logit(AObj$par[1], max = 1-10^-16, min = 10^-16))
              likvec <- c(likvec, LogLikComp_het_CnT(c(AObj$par[1], AObj$par[2], AObj$par[3]), ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, sprv=sprv, isCase = data_counts$isCase))
          
            }
          }
          
          theta_het_loc <- which(likvec == min(likvec))[1]
          theta_het_control <- thetavec_control[theta_het_loc]
          theta_het_tumor <- thetavec_tumor[theta_het_loc]
          probshift <- pivec[theta_het_loc]
          
        }
        
      } else{
        
        thetavec_control <- c()
        thetavec_tumor <- c()
        pivec <- c()
        likvec <- c()
        for(TH1 in thetaTRY){
          for(TH2 in thetaTRY){
          
            OptObj <- tryCatch( {optim(par = c(gtools::logit(probshift), log(TH1), log(TH2)), fn = LogLikComp_het_CnT, gr = GradComp_het_CnT, method = "BFGS",
                                       ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, isCase = data_counts$isCase)},
                                error = function(e) NULL)
            if(is.null(OptObj)){
              OptObj <- optim(par = c(gtools::logit(probshift), log(TH1), log(TH2)), fn = LogLikComp_het_CnT, method = "BFGS",
                              ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, isCase = data_counts$isCase)
            }
            thetavec_control <- c(thetavec_control, exp(OptObj$par[2]))
            thetavec_tumor <- c(thetavec_tumor, exp(OptObj$par[3]))
            pivec <- c(pivec, gtools::inv.logit(OptObj$par[1], max = 1-10^-16, min = 10^-16))
            likvec <- c(likvec, LogLikComp_het_CnT(c(OptObj$par[1], OptObj$par[2], OptObj$par[3]), ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, sprv=sprv, isCase = data_counts$isCase))
        
          }
        }
        
        theta_het_loc <- which(likvec == min(likvec))[1]
        theta_het_control <- thetavec_control[theta_het_loc]
        theta_het_tumor <- thetavec_tumor[theta_het_loc]
        probshift <- pivec[theta_het_loc]
        
      }
      
      
    }
    
    if(probshift > (1-10*SE)){
      probshift <- (1-10*SE)
    } else if(probshift < 10*SE){
      probshift <- 10*SE
    }
    
    Q <- sum(ifelse(pr>0, spr*log(pr), 0) + spr*dBetaBinom(data_counts$ref_count, data_counts$ref_count+data_counts$var_count, pi = 1-SE, theta = theta_hom, LOG = TRUE) +
               ifelse(prv>0, sprv*log(prv), 0) + sprv*c(dBetaBinom(data_counts$ref_count[data_counts$isCase==0], (data_counts$ref_count+data_counts$var_count)[data_counts$isCase==0], pi = probshift, theta = theta_het_control, LOG = TRUE),
                                                        dBetaBinom(data_counts$ref_count[data_counts$isCase==1], (data_counts$ref_count+data_counts$var_count)[data_counts$isCase==1], pi = probshift, theta = theta_het_tumor, LOG = TRUE))+
               ifelse(pv>0, spv*log(pv), 0) + spv*dBetaBinom(data_counts$var_count, data_counts$ref_count+data_counts$var_count, pi = 1-SE, theta = theta_hom, LOG = TRUE))
    
    dlta <- abs(Qold-Q)
    
    allelefreq <- mean(spr) + mean(sprv)/2
    data_counts$allelefreq <- allelefreq
  }
  
  #############################################################################################################################################################
  #############################################################################################################################################################
  
  GenoDF <- data.frame("pvv" = spv, "prr" = spr, "prv" = sprv)
  genotypeN <- c("rr", "rv", "vv")[apply(GenoDF[, c("prr", "prv", "pvv")], 1, function(x) which(x==max(x))[1])]
  
  ParamVec <- c(pr, prv, pv, probshift, theta_hom, theta_het_control, theta_het_tumor)
  names(ParamVec) <- c("pr", "prv", "pv", "probshift", "theta_hom", "theta_het_control", "theta_het_case")
  
  return(list("ParamVec" = ParamVec, "GenoDF" = GenoDF, "genotypeN" = genotypeN, "nrep" = nrep, "quality" = quality))
}




