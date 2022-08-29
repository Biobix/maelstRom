#' Models allelic RNAseq counts using an expectation-maximization fit of a beta-binomial mixture distribution 
#'
#' \code{EMfit_betabinom} estimates, per locus, parameters of an assumed beta-binomial mixture model via expectation maximization (genotype
#' frequencies, allelic bias i.e. heterozygous pi parameter, homozygous and heterozygous allelic divergence i.e. their respective theta- (overdispersion) parameters.
#' Through this fit, per-samples genotype probabilities are obtained as well. Optionally, a fit assuming no allelic bias (heterozygous pi-parameter = 0.5) is
#' performed as well, for the purpose of significant allelic bias detection via likelihood ratio test.
#' 
#' @param data_counts Data frame. Data frame of a SNP with reference and variant counts ("ref_count" and
#'     "var_count", respectively) for each sample ("sample").
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
#' @param fitH0 Logical. If TRUE, performs a fit assuming no allelic bias (heterozygous pi-parameter = 0.5) and a subsequent likihood ratio test for significant allelic bias.
#' @export
#' @return A list containing the following components:
#' \item{AB}{The estimated Allelic Bias (estimated reference read fraction for heterozygous samples).}
#' \item{AB_lrt}{The test statistic of the likelihood ratio test against perfectly balanced Allelic Bias (0.5). Not included if \code{fitH0==FALSE}.}
#' \item{AB_p}{The p-value of the likelihood ratio test against perfectly balanced Allelic Bias (0.5). Not included if \code{fitH0==FALSE}.}
#' \item{nrep}{The number of iterations.}
#' \item{quality}{Indicates the quality of a locus. An "!" indicates bad data or a bad fit due to no apparent heterozygosity.}
#' \item{rho_vv}{The variant homozygote genotype frequency for the locus.}
#' \item{rho_rr}{The reference homozygote genotype frequency for the locus.}
#' \item{rho_rv}{The heterozyous genotype frequency for the locus.}
#' \item{rho_vv_H0}{The variant homozygote genotype frequency with Allelic Bias = 0.5. Not included if \code{fitH0==FALSE}.}
#' \item{rho_rr_H0}{The reference homozygote genotype frequency with Allelic Bias = 0.5. Not included if \code{fitH0==FALSE}.}
#' \item{rho_rv_H0}{The heterozyous genotype frequency with Allelic Bias = 0.5. Not included if \code{fitH0==FALSE}.}
#' \item{data_hash}{Data frame. Input data frame with extra columns: allelefreq, (most likely) genotype, genotype probabilities (prr, prv, pvv), and whether the data point is an outlier ("Outlier", equals 1 if so) per sample.}
#' \item{theta_hom}{Final overdispersion estimate for the homozygous peaks.}
#' \item{theta_het}{Final overdispersion estimate for the heterozygous peak.}
#' \item{theta_hom_H0}{Final overdispersion estimate for the homozygous peaks with Allelic Bias = 0.5. Not included if \code{fitH0==FALSE}.}
#' \item{theta_het_H0}{Final overdispersion estimate for the heterozygous peak with Allelic Bias = 0.5. Not included if \code{fitH0==FALSE}.}
#' \item{GOF, GOFaltMEAN, GOFaltMEDIAN, GOFaltPERDIST, GOFaltPERDIST_S, GOFaltONLYHET, GOFexactMEAN, GOFexactMEANLOG}{Various Goodness-Of-Fit heuristics.}

EMfit_betabinom <- function(data_counts, allelefreq=0.5, SE, inbr = 0, dltaco = 10^-6, HWE = FALSE, p_InitEst = FALSE, 
                               ThetaInits = "moment", ReEstThetas = "moment", NoSplitHom = TRUE, NoSplitHet = TRUE,
                               ResetThetaMin = 10^-10, ResetThetaMax = 10^-1, thetaTRY = c(10^-1, 10^-3, 10^-7),
                               fitH0 = TRUE) {

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
  
  probshift <- 0.5 # Initial estimate position heterozygous peak
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
    theta_het <- ResetThetaMin
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
    MomentEst <- MomentEst_MixedBetaBinom(ref_counts = data_counts$ref_count, var_counts = data_counts$var_count, spr = spr_est, spv = spv_est, sprv = sprv_est, 
                                          pi_hom_fix = 1-SE, pi_het_fix = 0.5)
    theta_hom <- min(max(MomentEst["theta_hom"], ResetThetaMin), ResetThetaMax)
    theta_het <- min(max(MomentEst["theta_het"], ResetThetaMin), ResetThetaMax)
    if(is.nan(theta_het) | is.infinite(theta_het)){
      theta_het <- ResetThetaMin
    }
    if(is.nan(theta_hom) | is.infinite(theta_hom)){
      theta_hom <- ResetThetaMin
    }
    theta_hom_init <- theta_hom
    theta_het_init <- theta_het
  } else{
    theta_hom <- min(max(ThetaInits[1], ResetThetaMin), ResetThetaMax)
    theta_het <- min(max(ThetaInits[2], ResetThetaMin), ResetThetaMax)
  }
  
  theta_hom <- min(c(theta_hom, max((1-SE)/10, SE/10)) ) # Take some distance from the boundary at which bimodality occurs
  theta_het <- min(c(theta_het, max((1-probshift)/10, probshift/10)) ) # Take some distance from the boundary at which bimodality occurs
  
  theta_hom_H0 <- theta_hom
  theta_het_H0 <- theta_het
  
  quality <- ""
  
  ########################################################################################################################################################################
  ########################################################################################################################################################################
  
  Q <- 1000
    
  while (dlta > dltaco & nrep < 100) {
    Qold <- Q
    nrep <- nrep + 1
    
    spr <- pr * dBetaBinom(data_counts$ref_count, data_counts$ref_count + data_counts$var_count, pi = 1 - SE, theta = theta_hom, LOG = FALSE)
    spv <- pv * dBetaBinom(data_counts$var_count, data_counts$ref_count + data_counts$var_count, pi = 1 - SE, theta = theta_hom, LOG = FALSE)
    sprv <- prv * dBetaBinom(data_counts$ref_count, data_counts$ref_count + data_counts$var_count, pi = probshift, theta = theta_het, LOG = FALSE)
    
    pdata <- rowSums(cbind(spv, sprv, spr)) 
    
    if (any(pdata==0)){ 
      ProblemCases <- which(pdata==0)
      for(case in ProblemCases){
        var_part<-data_counts$var_count[case]
        ref_part<-data_counts$ref_count[case]
        spr_part<-dBetaBinom(ref_part, ref_part+var_part, pi = 1-SE, theta = theta_hom, LOG = TRUE)
        spv_part<-dBetaBinom(var_part, ref_part+var_part, pi = 1-SE, theta = theta_hom, LOG = TRUE)
        sprv_part<-dBetaBinom(ref_part, ref_part+var_part, pi = probshift, theta = theta_het, LOG = TRUE)
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
      MomentEst <- MomentEst_MixedBetaBinom(ref_counts = data_counts$ref_count, var_counts = data_counts$var_count, spr = spr, spv = spv, sprv = sprv, 
                                            pi_hom_fix = 1-SE, pi_het_fix = probshift)
      if(sum(spr+spv)!=0){
        theta_hom <- MomentEst["theta_hom"]
      }
      theta_het <- MomentEst["theta_het"]
      if(is.nan(theta_het) | is.infinite(theta_het) | theta_het < ResetThetaMin){
        theta_het <- ResetThetaMin
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
        theta_hom <- mean(c(10^(-5:-1), theta_hom)[WhichOne])
      }
      theta_het_vec <- c()
      for(ii in c(10^(-5:-1), theta_het)){
        Likely <- -LogLikComp_het(c(gtools::logit(probshift), log(ii)), ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv)
        theta_het_vec <- c(theta_het_vec, Likely)
      }
      WhichOne <- which(theta_het_vec==max(theta_het_vec))
      theta_het <- mean(c(10^(-5:-1), theta_het)[WhichOne])
      
    } else if(ReEstThetas == "simple" | nrep == 1){
      if(sum(spr+spv)!=0){
        theta_hom_vec <- c()
        for(ii in 10^(-5:-1)){
          Likely <- -LogLikComp_hom(log(ii), SE = SE, ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr, spv = spv)
          theta_hom_vec <- c(theta_hom_vec, Likely)
        }
        WhichOne <- which(theta_hom_vec==max(theta_hom_vec))
        theta_hom <- mean(10^(-5:-1)[WhichOne])
      }
      theta_het_vec <- c()
      for(ii in 10^(-5:-1)){
        Likely <- -LogLikComp_het(c(gtools::logit(probshift), log(ii)), ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv)
        theta_het_vec <- c(theta_het_vec, Likely)
      }
      WhichOne <- which(theta_het_vec==max(theta_het_vec))
      theta_het <- mean(10^(-5:-1)[WhichOne])
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
      theta_het_clone <- theta_het
      OptObj <- tryCatch( {optim(par = c(gtools::logit(probshift), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, gr = GradComp_het, method = "BFGS",
                                 ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv)},
                          error = function(e) NULL)
      if(is.null(OptObj)){
        OptObj <- optim(par = c(gtools::logit(probshift), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, method = "BFGS",
                        ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv)
      }
      
      probshift <- gtools::inv.logit(OptObj$par[1], max = 1-10^-16, min = 10^-16)
      theta_het <- exp(OptObj$par[2])
      
      if(theta_het  > max(probshift, 1-probshift)){
        
        theta_het <- min(c(theta_het_clone, max((1-probshift)/10, probshift/10)) ) # Take some distance from the boundary at which bimodality occurs
        
        AObj <- tryCatch( {alabama::auglag(c(gtools::logit(probshift), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, gr = GradComp_het, hin = eval_g_f_2,
                                ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, control.outer = list("trace" = FALSE))},
                          error = function(e) NULL)
        if(is.null(AObj)){
          AObj <- alabama::auglag(c(gtools::logit(probshift), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, hin = eval_g_f_2,
                                  ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, control.outer = list("trace" = FALSE))
        }
        
        probshift <- gtools::inv.logit(AObj$par[1], max = 1-10^-16, min = 10^-16)
        theta_het <- exp(AObj$par[2])
        
      }
    } else{
      OptObj <- tryCatch( {optim(par = c(gtools::logit(probshift), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, gr = GradComp_het, method = "BFGS",
                                 ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv)},
                          error = function(e) NULL)
      if(is.null(OptObj)){
        OptObj <- optim(par = c(gtools::logit(probshift), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, method = "BFGS",
                        ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv)
      }
      
      probshift <- gtools::inv.logit(OptObj$par[1], max = 1-10^-16, min = 10^-16)
      theta_het <- exp(OptObj$par[2])
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
        
        thetavec <- c()
        pivec <- c()
        likvec <- c()
        for(TH in thetaTRY){
          
          OptObj <- tryCatch( {optim(par = c(gtools::logit(probshift), log(TH)), fn = LogLikComp_het, gr = GradComp_het, method = "BFGS",
                                     ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv)},
                              error = function(e) NULL)
          if(is.null(OptObj)){
            OptObj <- optim(par = c(gtools::logit(probshift), log(TH)), fn = LogLikComp_het, method = "BFGS",
                            ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv)
          }
          thetavec <- c(thetavec, exp(OptObj$par[2]))
          pivec <- c(pivec, gtools::inv.logit(OptObj$par[1], max = 1-10^-16, min = 10^-16))
          likvec <- c(likvec, LogLikComp_het(c(OptObj$par[1], OptObj$par[2]), ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, sprv=sprv))
        }
        theta_het_loc <- which(likvec == min(likvec))[1]
        theta_het <- thetavec[theta_het_loc]
        probshift <- pivec[theta_het_loc]
        
        if(theta_het  > max(probshift, 1-probshift)){
          
          thetavec <- c()
          pivec <- c()
          likvec <- c()
          for(TH in thetaTRY){
            
            AObj <- tryCatch( {alabama::auglag(c(gtools::logit(probshift), log(TH)), fn = LogLikComp_het, gr = GradComp_het, hin = eval_g_f_2,
                                               ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, control.outer = list("trace" = FALSE))},
                              error = function(e) NULL)
            if(is.null(AObj)){
              AObj <- alabama::auglag(c(gtools::logit(probshift), log(TH)), fn = LogLikComp_het, hin = eval_g_f_2,
                                      ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, control.outer = list("trace" = FALSE))
            }
            thetavec <- c(thetavec, exp(AObj$par[2]))
            pivec <- c(pivec, gtools::inv.logit(AObj$par[1], max = 1-10^-16, min = 10^-16))
            likvec <- c(likvec, LogLikComp_het(c(AObj$par[1], AObj$par[2]), ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, sprv=sprv))
          }
            
          theta_het_loc <- which(likvec == min(likvec))[1]
          theta_het <- thetavec[theta_het_loc]
          probshift <- pivec[theta_het_loc]
          
        }
      } else{
        
        thetavec <- c()
        pivec <- c()
        likvec <- c()
        for(TH in thetaTRY){
          OptObj <- tryCatch( {optim(par = c(gtools::logit(probshift), log(TH)), fn = LogLikComp_het, gr = GradComp_het, method = "BFGS",
                                     ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv)},
                              error = function(e) NULL)
          if(is.null(OptObj)){
            OptObj <- optim(par = c(gtools::logit(probshift), log(TH)), fn = LogLikComp_het, method = "BFGS",
                            ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv)
          }
          thetavec <- c(thetavec, exp(OptObj$par[2]))
          pivec <- c(pivec, gtools::inv.logit(OptObj$par[1], max = 1-10^-16, min = 10^-16))
          likvec <- c(likvec, LogLikComp_het(c(OptObj$par[1], OptObj$par[2]), ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, sprv=sprv))
        }
        
        theta_het_loc <- which(likvec == min(likvec))[1]
        theta_het <- thetavec[theta_het_loc]
        probshift <- pivec[theta_het_loc]
        
      }
      
      
    }
    
    if(probshift > (1-10*SE)){
      probshift <- (1-10*SE)
    } else if(probshift < 10*SE){
      probshift <- 10*SE
    }

    Q <- sum(ifelse(pr>0, spr*log(pr), 0) + spr*dBetaBinom(data_counts$ref_count, data_counts$ref_count+data_counts$var_count, pi = 1-SE, theta = theta_hom, LOG = TRUE) +
               ifelse(prv>0, sprv*log(prv), 0) + sprv*dBetaBinom(data_counts$ref_count, data_counts$ref_count+data_counts$var_count, pi = probshift, theta = theta_het, LOG = TRUE) +
               ifelse(pv>0, spv*log(pv), 0) + spv*dBetaBinom(data_counts$var_count, data_counts$ref_count+data_counts$var_count, pi = 1-SE, theta = theta_hom, LOG = TRUE))
    
    dlta <- abs(Qold-Q)
    
    allelefreq <- mean(spr) + mean(sprv)/2
    data_counts$allelefreq <- allelefreq
  }
  
  #############################################################################################################################################################
  #############################################################################################################################################################

  data_counts$pvv <- spv; data_counts$prr <- spr; data_counts$prv <- sprv
  data_counts$genotypeN <- c("rr", "rv", "vv")[apply(data_counts[, c("prr", "prv", "pvv")], 1, function(x) which(x==max(x))[1])]
  
  if(!fitH0){
    results <- list(AB = probshift, nrep = nrep, quality = quality, rho_rr = pr, rho_vv = pv, rho_rv = prv, data_hash = data_counts, 
                    theta_hom = theta_hom, theta_het = theta_het)
    return(results)
  }
  
  # A whole bunch of potential GOF-measures
  data_counts$FitRR<-0
  data_counts$FitRV<-0
  data_counts$FitVV<-0
  data_counts$ExactPVAL<-0
  data_counts$TotPVAL<-0
  for(DCR in 1:nrow(data_counts)){
    data_counts$FitRR[DCR] <- (BetaBinom_test_pvalTS(m = data_counts$ref_count[DCR], n = data_counts$ref_count[DCR] + data_counts$var_count[DCR], pi = 1-SE, theta = theta_hom))
    data_counts$FitRV[DCR] <- (BetaBinom_test_pvalTS(m = data_counts$ref_count[DCR], n = data_counts$ref_count[DCR] + data_counts$var_count[DCR], pi = probshift, theta = theta_het))
    data_counts$FitVV[DCR] <- (BetaBinom_test_pvalTS(m = data_counts$var_count[DCR], n = data_counts$ref_count[DCR] + data_counts$var_count[DCR], pi = 1-SE, theta = theta_hom))
    data_counts$ExactPVAL[DCR] <- data_counts$prr[DCR]*data_counts$FitRR[DCR] +
      data_counts$prv[DCR]*data_counts$FitRV[DCR] +
      data_counts$pvv[DCR]*data_counts$FitVV[DCR]
    # A TOTAL exact test
    TotC <- data_counts$ref_count[DCR] + data_counts$var_count[DCR]
    CountCombos <- cbind(0:TotC, TotC:0)
    
    ComboProbs <- pr * dBetaBinom(CountCombos[,1], rep(TotC, TotC+1), pi = 1 - SE, theta = theta_hom, LOG = FALSE) + pv * dBetaBinom(CountCombos[,2], rep(TotC, TotC+1), pi = 1 - SE, theta = theta_hom, LOG = FALSE) + prv * dBetaBinom(CountCombos[,1], rep(TotC, TotC+1), pi = probshift, theta = theta_het, LOG = FALSE)
    
    PROBOI <- ComboProbs[(data_counts$ref_count[DCR]+1)]
    TOTPVAL <- sum(ComboProbs[ComboProbs<=PROBOI])
    data_counts$TotPVAL[DCR] <- TOTPVAL
  }
  
  GOFaltMEAN_S <- mean(data_counts$ExactPVAL)
  GOFaltMEDIAN_S <- median(data_counts$ExactPVAL)
  GOFaltPERDIST_S <- (sum(data_counts$prr*data_counts$FitRR)/sum(data_counts$prr) + sum(data_counts$prv*data_counts$FitRV)/sum(data_counts$prv) + sum(data_counts$pvv*data_counts$FitVV)/sum(data_counts$pvv))/3
  GOFaltONLYHET_S <- sum(data_counts$prv*data_counts$FitRV)/sum(data_counts$prv)
  
  GOFexactMEAN_S <- mean(data_counts$TotPVAL)
  GOFexactMEANLOG_S <- mean(log(data_counts$TotPVAL))
  
  # The old GOF-heuristic
  dmixase_corrected <- pmf_betabinomMix(data_counts$ref_count, data_counts$var_count, probshift, SE, mean(spr), mean(spv), mean(sprv), theta_hom = theta_hom, theta_het = theta_het) * (data_counts$ref_count + data_counts$var_count + 1)
  logLikelihood <- mean(log(dmixase_corrected))
  #probshift_adj <- -(as.numeric(probshift) - 0.5)/0.5
  

  # For the LRT, we should actually re-fit the unshifted data while keeping the sequencing error constant...
  # So this means an additional EM run!
  
  dlta <- 1
  nrep2 <- 0
  MyFlag <- ""
  Q <- 1000
  # We already have initial theta estimates under H0
  
  while (dlta > dltaco & nrep2 < 100) {
    Qold <- Q
    nrep2 <- nrep2 + 1
    # Beta binom:
    spr_H0 <- pr_H0 * dBetaBinom(data_counts$ref_count, data_counts$ref_count + data_counts$var_count, pi = 1 - SE, theta = theta_hom_H0, LOG = FALSE)
    spv_H0 <- pv_H0 * dBetaBinom(data_counts$var_count, data_counts$ref_count + data_counts$var_count, pi = 1 - SE, theta = theta_hom_H0, LOG = FALSE)
    sprv_H0 <- prv_H0 * dBetaBinom(data_counts$ref_count, data_counts$ref_count + data_counts$var_count, pi = 0.5, theta = theta_het_H0, LOG = FALSE)
    
    pdata <- rowSums(cbind(spr_H0, spv_H0, sprv_H0)) 
    
    if (any(pdata==0)){ 
      ProblemCases <- which(pdata==0)
      for(case in ProblemCases){
        var_part<-data_counts$var_count[case]
        ref_part<-data_counts$ref_count[case]
        spr_part<-dBetaBinom(ref_part, ref_part+var_part, pi = 1-SE, theta = theta_hom_H0, LOG = TRUE)
        spv_part<-dBetaBinom(var_part, ref_part+var_part, pi = 1-SE, theta = theta_hom_H0, LOG = TRUE)
        sprv_part<-dBetaBinom(ref_part, ref_part+var_part, pi = 0.5, theta = theta_het_H0, LOG = TRUE)
        spvec<-c(spr_part, sprv_part, spv_part)
        if(spr_part==max(spvec)){
          spr_H0[case]<-1
          pdata[case]<-1
        } else if(sprv_part==max(spvec)){
          sprv_H0[case]<-1
          pdata[case]<-1
        } else{
          spv_H0[case]<-1
          pdata[case]<-1
        }
      }
    }
    
    #if(sum(spr_H0+spv_H0)==0){
    #  break
    #}
    
    spv_H0 <- spv_H0/pdata
    sprv_H0 <- sprv_H0/pdata
    spr_H0 <- spr_H0/pdata
    
    if (sum(sprv_H0)==0) {
      MyFlag <- "!"
      break
    } 

    allelefreq_H0 <- mean(spr_H0) + mean(sprv_H0)/2 
    if (HWE) {
      prv_H0 <- 2 * allelefreq_H0 * (1 - allelefreq_H0) * (1 - inbr)
      pr_H0 <- allelefreq_H0^2 + inbr * allelefreq_H0 * (1 - allelefreq_H0)
      pv_H0 <- (1 - allelefreq_H0)^2 + inbr * allelefreq_H0 * (1 - allelefreq_H0)
    }else {
      pv_H0 <- mean(spv_H0) 
      prv_H0 <- mean(sprv_H0)
      pr_H0 <- mean(spr_H0)
    }
    
    if(ReEstThetas == "moment"){
      MomentEst <- MomentEst_MixedBetaBinom(ref_counts = data_counts$ref_count, var_counts = data_counts$var_count, spr = spr_H0, spv = spv_H0, sprv = sprv_H0, 
                                            pi_hom_fix = 1-SE, pi_het_fix = 0.5)
      if(sum(spr_H0+spv_H0)!=0){
        theta_hom_H0 <- MomentEst["theta_hom"]
      }
      theta_het_H0 <- MomentEst["theta_het"]
      if(is.nan(theta_het_H0) | is.infinite(theta_het_H0) | theta_het_H0 < ResetThetaMin){
        theta_het_H0 <- ResetThetaMin
      }
      if(is.nan(theta_hom_H0) | is.infinite(theta_hom_H0) | theta_hom_H0 < ResetThetaMin){
        theta_hom_H0 <- ResetThetaMin
      }
      
      if(sum(spr_H0+spv_H0)!=0){
        theta_hom_vec <- c()
        for(ii in c(10^(-5:-1), theta_hom_H0) ){
          Likely <- -LogLikComp_hom(log(ii), SE = SE, ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, spv = spv_H0)
          theta_hom_vec <- c(theta_hom_vec, Likely)
        }
        WhichOne <- which(theta_hom_vec==max(theta_hom_vec))
        theta_hom_H0 <- mean(c(10^(-5:-1), theta_hom_H0)[WhichOne])
      }
      theta_het_vec <- c()
      for(ii in c(10^(-5:-1), theta_het_H0)){
        Likely <- -LogLikComp_het_H0(log(ii), probshift = 0.5, ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0)
        theta_het_vec <- c(theta_het_vec, Likely)
      }
      WhichOne <- which(theta_het_vec==max(theta_het_vec))
      theta_het_H0 <- mean(c(10^(-5:-1), theta_het_H0)[WhichOne])
      
    } else if(ReEstThetas == "simple" | nrep2 == 1){
      if(sum(spr_H0+spv_H0)!=0){
        theta_hom_vec <- c()
        for(ii in 10^(-5:-1)){
          Likely <- -LogLikComp_hom(log(ii), SE = SE, ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, spv = spv_H0)
          theta_hom_vec <- c(theta_hom_vec, Likely)
        }
        WhichOne <- which(theta_hom_vec==max(theta_hom_vec))
        theta_hom_H0 <- mean(10^(-5:-1)[WhichOne])
      }
      theta_het_vec <- c()
      for(ii in 10^(-5:-1)){
        Likely <- -LogLikComp_het_H0(log(ii), probshift = 0.5, ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0)
        theta_het_vec <- c(theta_het_vec, Likely)
      }
      WhichOne <- which(theta_het_vec==max(theta_het_vec))
      theta_het_H0 <- mean(10^(-5:-1)[WhichOne])
    }
    
    if(ReEstThetas != "TryThree"){
    if(sum(spr_H0+spv_H0)!=0){
      if(NoSplitHom){
        theta_hom_H0_clone <- theta_hom_H0
        OptObj <- optim(par = log(min(max(theta_hom_H0, ResetThetaMin), ResetThetaMax)), fn = LogLikComp_hom, gr = GradComp_hom, method = "BFGS", SE = SE, 
                        ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, spv = spv_H0)
        theta_hom_H0 <- exp(OptObj$par)
        
        if(theta_hom_H0 > max(SE, 1-SE)){
          theta_hom_H0 <- min(c(theta_hom_H0_clone, max((1-SE)/10, SE/10)) ) # Take some distance from the boundary at which bimodality occurs
          
          AObj <- alabama::auglag(log(min(max(theta_hom_H0, ResetThetaMin), ResetThetaMax)), fn = LogLikComp_hom, gr = GradComp_hom, hin = eval_g_f, SE = SE,
                                  ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, spv = spv_H0, control.outer = list("trace" = FALSE))
          
          theta_hom_H0 <- exp(AObj$par[1])
        }
        
      } else{
        OptObj <- optim(par = log(min(max(theta_hom_H0, ResetThetaMin), ResetThetaMax)), fn = LogLikComp_hom, gr = GradComp_hom, method = "BFGS", SE = SE, 
                        ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, spv = spv_H0)
        theta_hom_H0 <- exp(OptObj$par)
      }
    }
    
    
    if(NoSplitHet){
      theta_het_H0_clone <- theta_het_H0
      OptObj <- tryCatch( {optim(par = log(min(max(theta_het_H0, ResetThetaMin), ResetThetaMax)), probshift = 0.5, fn = LogLikComp_het_H0, gr = GradComp_het_H0, method = "BFGS",
                                 ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0)},
                          error = function(e) NULL)
      if(is.null(OptObj)){
        OptObj <- optim(par = log(min(max(theta_het_H0, ResetThetaMin), ResetThetaMax)), probshift = 0.5, fn = LogLikComp_het_H0, method = "BFGS",
                        ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0)
      }
      
      theta_het_H0 <- exp(OptObj$par[1])
      
      if(theta_het_H0  > 0.5){
        
        theta_het_H0 <- min(c(theta_het_H0_clone, 0.5/10) ) # Take some distance from the boundary at which bimodality occurs
        
        AObj <- tryCatch( {alabama::auglag(log(min(max(theta_het_H0, ResetThetaMin), ResetThetaMax)), probshift = 0.5, fn = LogLikComp_het_H0, gr = GradComp_het_H0, hin = eval_g_f_2_H0,
                                           ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0, control.outer = list("trace" = FALSE))},
                          error = function(e) NULL)
        if(is.null(AObj)){
          AObj <- alabama::auglag(log(min(max(theta_het_H0, ResetThetaMin), ResetThetaMax)), probshift = 0.5, fn = LogLikComp_het_H0, hin = eval_g_f_2_H0,
                                  ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0, control.outer = list("trace" = FALSE))
        }

        theta_het_H0 <- exp(AObj$par[1])
        
      }
    } else{
      OptObj <- tryCatch( {optim(par = log(min(max(theta_het_H0, ResetThetaMin), ResetThetaMax)), probshift = 0.5, fn = LogLikComp_het_H0, gr = GradComp_het_H0, method = "BFGS",
                                 ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0)},
                          error = function(e) NULL)
      if(is.null(OptObj)){
        OptObj <- optim(par = log(min(max(theta_het_H0, ResetThetaMin), ResetThetaMax)), probshift = 0.5, fn = LogLikComp_het_H0, method = "BFGS",
                        ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0)
      }
      
      theta_het_H0 <- exp(OptObj$par[1])
    }
    } else{
      
      if(sum(spr_H0+spv_H0)!=0){
        if(NoSplitHom){
          thetavec <- c()
          likvec <- c()
          for(TH in thetaTRY){
            OptObj <- optim(par = log(TH), fn = LogLikComp_hom, gr = GradComp_hom, method = "BFGS", SE = SE, 
                            ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, spv = spv_H0)
            thetavec <- c(thetavec, exp(OptObj$par))
            likvec <- c(likvec, LogLikComp_hom(OptObj$par, SE=SE, ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, spr=spr_H0, spv=spv_H0))
          }
          theta_hom_loc <- which(likvec == min(likvec))[1]
          theta_hom_H0 <- thetavec[theta_hom_loc]

          if(theta_hom_H0 > max(SE, 1-SE)){
            thetavec <- c()
            likvec <- c()
            for(TH in thetaTRY){
              AObj <- alabama::auglag(log(TH), fn = LogLikComp_hom, gr = GradComp_hom, hin = eval_g_f, SE = SE,
                                      ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, spv = spv_H0, control.outer = list("trace" = FALSE))
              thetavec <- c(thetavec, exp(AObj$par[1]))
              likvec <- c(likvec, LogLikComp_hom(AObj$par[1], SE=SE, ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, spr=spr_H0, spv=spv_H0))
            }

            theta_hom_loc <- which(likvec == min(likvec))[1]
            theta_hom_H0 <- thetavec[theta_hom_loc]
          }
          
        } else{
          thetavec <- c()
          likvec <- c()
          for(TH in thetaTRY){
            OptObj <- optim(par = log(TH), fn = LogLikComp_hom, gr = GradComp_hom, method = "BFGS", SE = SE, 
                            ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, spv = spv_H0)
            thetavec <- c(thetavec, exp(OptObj$par))
            likvec <- c(likvec, LogLikComp_hom(OptObj$par, SE=SE, ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, spr=spr_H0, spv=spv_H0))
          }
          theta_hom_loc <- which(likvec == min(likvec))[1]
          theta_hom_H0 <- thetavec[theta_hom_loc]
        }
      }
      
      
      if(NoSplitHet){
        thetavec <- c()
        likvec <- c()
        for(TH in thetaTRY){
          OptObj <- tryCatch( {optim(par = log(TH), probshift = 0.5, fn = LogLikComp_het_H0, gr = GradComp_het_H0, method = "BFGS",
                                     ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0)},
                              error = function(e) NULL)
          if(is.null(OptObj)){
            OptObj <- optim(par = log(TH), probshift = 0.5, fn = LogLikComp_het_H0, method = "BFGS",
                            ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0)
          }
          thetavec <- c(thetavec, exp(OptObj$par[1]))
          likvec <- c(likvec, LogLikComp_het_H0(OptObj$par[1], probshift=0.5, ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, sprv=sprv_H0))
        }
        theta_het_loc <- which(likvec == min(likvec))[1]
        theta_het_H0 <- thetavec[theta_het_loc]

        if(theta_het_H0  > 0.5){
          thetavec <- c()
          likvec <- c()
          for(TH in thetaTRY){
            AObj <- tryCatch( {alabama::auglag(log(TH), probshift = 0.5, fn = LogLikComp_het_H0, gr = GradComp_het_H0, hin = eval_g_f_2_H0,
                                               ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0, control.outer = list("trace" = FALSE))},
                              error = function(e) NULL)
            if(is.null(AObj)){
              AObj <- alabama::auglag(log(TH), probshift = 0.5, fn = LogLikComp_het_H0, hin = eval_g_f_2_H0,
                                      ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0, control.outer = list("trace" = FALSE))
            }
            thetavec <- c(thetavec, exp(AObj$par[1]))
            likvec <- c(likvec, LogLikComp_het_H0(AObj$par[1], probshift=0.5, ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, sprv=sprv_H0))
            
          }
          theta_het_loc <- which(likvec == min(likvec))[1]
          theta_het_H0 <- thetavec[theta_het_loc]
        }
      } else{
        thetavec <- c()
        likvec <- c()
        for(TH in thetaTRY){
          OptObj <- tryCatch( {optim(par = log(TH), probshift = 0.5, fn = LogLikComp_het_H0, gr = GradComp_het_H0, method = "BFGS",
                                     ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0)},
                              error = function(e) NULL)
          if(is.null(OptObj)){
            OptObj <- optim(par = log(TH), probshift = 0.5, fn = LogLikComp_het_H0, method = "BFGS",
                            ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0)
          }
          thetavec <- c(thetavec, exp(OptObj$par[1]))
          likvec <- c(likvec, LogLikComp_het_H0(OptObj$par[1], probshift=0.5, ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, sprv=sprv_H0))
        }
        
        theta_het_loc <- which(likvec == min(likvec))[1]
        theta_het_H0 <- thetavec[theta_het_loc]
      }
      
    }
    
    

    Q <- sum(ifelse(pr_H0>0, spr_H0*log(pr_H0), 0) + spr_H0*dBetaBinom(data_counts$ref_count, data_counts$ref_count+data_counts$var_count, pi = 1-SE, theta = theta_hom_H0, LOG = TRUE) +
               ifelse(prv_H0>0, sprv_H0*log(prv_H0), 0) + sprv_H0*dBetaBinom(data_counts$ref_count, data_counts$ref_count+data_counts$var_count, pi = 0.5, theta = theta_het_H0, LOG = TRUE) +
               ifelse(pv_H0>0, spv_H0*log(pv_H0), 0) + spv_H0*dBetaBinom(data_counts$var_count, data_counts$ref_count+data_counts$var_count, pi = 1-SE, theta = theta_hom_H0, LOG = TRUE))
    
    dlta <- abs(Qold-Q)
    
    allelefreq_H0 <- mean(spr_H0) + mean(sprv_H0)/2
  }
  
  
  if(MyFlag == "!"){
    lrtstat <- NA
    pval <- NA
  } else{
    dmixh0 <- pmf_betabinomMix(data_counts$ref_count, data_counts$var_count, probshift = 0.5, SE,  mean(spr_H0), mean(spv_H0), mean(sprv_H0) , theta_hom = theta_hom_H0, theta_het = theta_het_H0)
    dmixase <- pmf_betabinomMix(data_counts$ref_count, data_counts$var_count, probshift = probshift, SE, mean(spr), mean(spv), mean(sprv), theta_hom = theta_hom, theta_het = theta_het)
    lrtstat <- -2 * (sum(log(dmixh0)) - sum(log(dmixase)))
    pval <- pchisq(lrtstat, df = 1, lower.tail = F)
  }
  
  data_counts$prr_H0 <- spr_H0
  data_counts$pvv_H0 <- spv_H0
  data_counts$prv_H0 <- sprv_H0
  
  results <- list(AB = probshift, AB_lrt = lrtstat, AB_p = pval, GOF = logLikelihood, GOFaltMEAN = GOFaltMEAN_S, GOFaltMEDIAN = GOFaltMEDIAN_S, 
                  GOFaltPERDIST = GOFaltPERDIST_S, GOFaltONLYHET = GOFaltONLYHET_S, GOFexactMEAN = GOFexactMEAN_S, GOFexactMEANLOG = GOFexactMEANLOG_S, 
                  nrep = nrep, quality = quality, rho_rr = pr, rho_vv = pv, rho_rv = prv, rho_rr_H0 = pr_H0, rho_vv_H0 = pv_H0, rho_rv_H0 = prv_H0, 
                  data_hash = data_counts, theta_hom = theta_hom, theta_het = theta_het, theta_hom_NoShift = theta_hom_H0, theta_het_NoShift = theta_het_H0)
  return(results)
}




