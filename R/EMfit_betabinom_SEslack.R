#' This function is a work in progress. As such, it is not exported yet; but it allows for some leniency on the sequencing error parameter, hence "SESlack"

EMfit_betabinom_SEslack <- function(data_counts, allelefreq=0.5, SE, inbr = 0, dltaco = 10^-6, HWE = FALSE, p_InitEst = FALSE, 
                               ThetaInits = "moment", ReEstThetas = "moment", NoSplitHom = TRUE, NoSplitHet = TRUE,
                               ResetThetaMin = 10^-10, ResetThetaMax = 10^-1, thetaTRY = c(10^-1, 10^-3, 10^-7),
                               fitH0 = TRUE, SESlack = 10, probshift_InitEst = TRUE) {
  
  LogLikComp_RR <- function(parvec, ref_counts, var_counts, spr, pim, M){
    return(-sum(spr*(dBetaBinom(ref_counts, ref_counts+var_counts, pi=gtools::inv.logit(parvec[1], max = M, min =pim), theta=exp(parvec[2]), LOG = TRUE))))
  }
  
  LogLikComp_VV <- function(parvec, ref_counts, var_counts, spv, pim, M){
    return(-sum(spv*(dBetaBinom(var_counts, ref_counts+var_counts, pi=gtools::inv.logit(parvec[1], max = M, min =pim), theta=exp(parvec[2]), LOG = TRUE))))
  }
  
  #LogLikComp_hom <- function(theta, SE, ref_counts, var_counts, spr, spv){
  #  return(-sum(spr*(dBetaBinom(ref_counts, ref_counts+var_counts, pi=1-SE, theta=exp(theta), LOG = TRUE)) +
  #                spv*(dBetaBinom(var_counts,ref_counts+var_counts, pi=1-SE, theta=exp(theta), LOG = TRUE))))
  #}
  
  GradComp_RR <- function(parvec, ref_counts, var_counts, spr,pim, M){
    Grad1 <- sum(spr*grad_pi(ref_counts, ref_counts+var_counts, pi=gtools::inv.logit(parvec[1], max = M, min =pim), theta=exp(parvec[2])) * (exp(parvec[1])/(1+exp(parvec[1]))^2) * (M-m) )
    Grad2 <- sum(spr*grad_theta(ref_counts, ref_counts+var_counts, pi=gtools::inv.logit(parvec[1], max = M, min =pim), theta=exp(parvec[2])) * exp(parvec[2]) )
    return(c(-Grad1, -Grad2))
  }
  
  GradComp_VV <- function(parvec, ref_counts, var_counts, spv,pim, M){
    Grad1 <- sum(spv*grad_pi(var_counts, ref_counts+var_counts, pi=gtools::inv.logit(parvec[1], max = M, min =pim), theta=exp(parvec[2])) * (exp(parvec[1])/(1+exp(parvec[1]))^2) * (M-m) )
    Grad2 <- sum(spv*grad_theta(var_counts, ref_counts+var_counts, pi=gtools::inv.logit(parvec[1], max = M, min =pim), theta=exp(parvec[2])) * exp(parvec[2]) )
    return(c(-Grad1, -Grad2))
  }
  
  #GradComp_hom <- function(theta, SE, ref_counts, var_counts, spr, spv){
  #  Grad <- sum(spr*grad_theta(ref_counts, ref_counts+var_counts, pi=1-SE, theta=exp(theta)) * exp(theta) +
  #                spv*grad_theta(var_counts, ref_counts+var_counts, pi=1-SE, theta=exp(theta)) * exp(theta))
  #  return(-Grad)
  #}
  
  LogLikComp_het <- function(parvec, ref_counts, var_counts, sprv, piRR, piVV){
    return(-sum(sprv*(dBetaBinom(ref_counts, ref_counts+var_counts, pi=gtools::inv.logit(parvec[1], max = piRR, min = 1-piVV), theta=exp(parvec[2]), LOG = TRUE))))
  }
  GradComp_het <- function(parvec, ref_counts, var_counts, sprv, piRR, piVV){
    Grad1 <- sum(sprv*grad_pi(ref_counts, ref_counts+var_counts, pi=gtools::inv.logit(parvec[1], max = piRR, min = 1-piVV), theta=exp(parvec[2])) * (exp(parvec[1])/(1+exp(parvec[1]))^2) * (piRR+piVV-1) )
    Grad2 <- sum(sprv*grad_theta(ref_counts, ref_counts+var_counts, pi=gtools::inv.logit(parvec[1], max = piRR, min = 1-piVV), theta=exp(parvec[2])) * exp(parvec[2]) )
    return(c(-Grad1, -Grad2))
  }
  
  LogLikComp_het_H0 <- function(theta, probshift, ref_counts, var_counts, sprv){
    return(-sum(sprv*(dBetaBinom(ref_counts, ref_counts+var_counts, pi=probshift, theta=exp(theta), LOG = TRUE))))
  }
  GradComp_het_H0 <- function(theta, probshift, ref_counts, var_counts, sprv){
    return(-sum(sprv*grad_theta(ref_counts, ref_counts+var_counts, pi=probshift, theta=exp(theta)) * exp(theta) ))
  }
  
  eval_g_f_RR <- function(parvec, ref_counts, var_counts, spr,pim, M) { # Specifies the inequality constraint between parameters in case NoSplitHet == TRUE
    return(( max(1 - gtools::inv.logit(parvec[1], max = M, min =pim), gtools::inv.logit(parvec[1], max = M, min =pim)) ) - exp(parvec[2]))
  }
  
  eval_g_f_VV <- function(parvec, ref_counts, var_counts, spv,pim, M) { # Specifies the inequality constraint between parameters in case NoSplitHet == TRUE
    return(( max(1 - gtools::inv.logit(parvec[1], max = M, min =pim), gtools::inv.logit(parvec[1], max = M, min =pim)) ) - exp(parvec[2]))
  }
  
  #eval_g_f <- function(theta, SE, ref_counts, var_counts, spr, spv) { # Specifies the inequality constraint between parameters in case NoSplitHom == TRUE
  #  return(( max(1 - SE, SE) ) - exp(theta))
  #}
  
  eval_g_f_2 <- function(parvec, ref_counts, var_counts, sprv, piRR, piVV) { # Specifies the inequality constraint between parameters in case NoSplitHet == TRUE
    return(( max(1 - gtools::inv.logit(parvec[1], max = piRR, min = 1-piVV), gtools::inv.logit(parvec[1], max = piRR, min = 1-piVV)) ) - exp(parvec[2]))
  }
  
  eval_g_f_2_H0 <- function(theta, probshift, ref_counts, var_counts, sprv) { # Specifies the inequality constraint between parameters in case NoSplitHet == TRUE
    return(( max(1 - probshift, probshift) ) - exp(theta))
  }
  
  
  piRR <- 1-SESlack/2*SE # Initial estimate position homozygous reference peak
  piVV <- 1-SESlack/2*SE # Initial estimate position homozygous variant peak
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
  
  if(probshift_InitEst){
    TotLikVec <- c()
    for(PS in seq(0.1,0.9,0.1)){
      TotLik <- sum(log(pmf_betabinomMix(data_counts$ref_count, data_counts$var_count, probshift = PS, SE = SE, pr = pr, pv = pv, prv = prv, theta_hom = 0, theta_het = 0)))
      TotLikVec <- c(TotLikVec, TotLik)
    }
    probshift <- mean(seq(0.1,0.9,0.1)[which(TotLikVec==max(TotLikVec))])
  }else{
    probshift <- 0.5 # Initial estimate position heterozygous peak
  }
  
  pr_H0 <- pr
  pv_H0 <- pv
  prv_H0 <- prv
  
  if(is.null(ThetaInits)){
    theta_RR <- ResetThetaMin
    theta_VV <- ResetThetaMin
    theta_het <- ResetThetaMin
  }else if(ThetaInits == "moment"){
    
    # First fit some regular binomial models to get an initial categorisation of the data points
    spr_est <- pr * dbinom(data_counts$ref_count, data_counts$ref_count + data_counts$var_count, prob = piRR)
    spv_est <- pv * dbinom(data_counts$var_count,data_counts$ref_count + data_counts$var_count, prob = piVV)
    sprv_est <- prv * dbinom(data_counts$ref_count, data_counts$ref_count + data_counts$var_count, prob = 0.5)
    
    pdata <- rowSums(cbind(spv_est, sprv_est, spr_est)) 
    if (any(pdata==0)){
      ProblemCases <- which(pdata==0)
      for(case in ProblemCases){
        var_part<-data_counts$var_count[case]
        ref_part<-data_counts$ref_count[case]
        spv_part<-dbinom(var_part, ref_part+var_part, piRR, log = TRUE)
        sprv_part<-dbinom(ref_part, ref_part+var_part, 0.5, log = TRUE)
        spr_part<-dbinom(ref_part, ref_part+var_part, piVV, log = TRUE)
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
    MomentEst <- MomentEst_MixedBetaBinom_SESlack(ref_counts = data_counts$ref_count, var_counts = data_counts$var_count, spr = spr_est, spv = spv_est, sprv = sprv_est, 
                                          pi_RR_fix = piRR, pi_VV_fix = piVV, pi_het_fix = 0.5)
    theta_RR <- min(max(MomentEst["theta_RR"], ResetThetaMin), ResetThetaMax)
    theta_VV <- min(max(MomentEst["theta_VV"], ResetThetaMin), ResetThetaMax)
    theta_het <- min(max(MomentEst["theta_het"], ResetThetaMin), ResetThetaMax)
    if(is.nan(theta_het) | is.infinite(theta_het)){
      theta_het <- ResetThetaMin
    }
    if(is.nan(theta_RR) | is.infinite(theta_RR)){
      theta_RR <- ResetThetaMin
    }
    if(is.nan(theta_VV) | is.infinite(theta_VV)){
      theta_VV <- ResetThetaMin
    }
    theta_RR_init <- theta_RR
    theta_VV_init <- theta_VV
    theta_het_init <- theta_het
  } else{
    theta_RR <- min(max(ThetaInits[1], ResetThetaMin), ResetThetaMax)
    theta_VV <- min(max(ThetaInits[2], ResetThetaMin), ResetThetaMax)
    theta_het <- min(max(ThetaInits[3], ResetThetaMin), ResetThetaMax)
    theta_RR_init <- theta_RR
    theta_VV_init <- theta_VV
    theta_het_init <- theta_het
  }
  
  theta_RR <- min(c(theta_RR, max((1-piRR)/10, piRR/10)) ) # Take some distance from the boundary at which bimodality occurs
  theta_VV <- min(c(theta_VV, max((1-piVV)/10, piVV/10)) ) # Take some distance from the boundary at which bimodality occurs
  theta_het <- min(c(theta_het, max((1-probshift)/10, probshift/10)) ) # Take some distance from the boundary at which bimodality occurs
  
  theta_RR_H0 <- theta_RR
  theta_VV_H0 <- theta_VV
  theta_het_H0 <- theta_het
  piRR_H0 <- piRR
  piVV_H0 <- piVV
  
  quality <- ""
  
  ########################################################################################################################################################################
  ########################################################################################################################################################################
  
  Q <- 1000
  maxC <- max(data_counts$ref_count + data_counts$var_count)
  eps <- 10^-15
  
  while (dlta > dltaco & nrep < 100) {
    Qold <- Q
    nrep <- nrep + 1
    
    spr <- pr * dBetaBinom(data_counts$ref_count, data_counts$ref_count + data_counts$var_count, pi = piRR, theta = theta_RR, LOG = FALSE)
    spv <- pv * dBetaBinom(data_counts$var_count, data_counts$ref_count + data_counts$var_count, pi = piVV, theta = theta_VV, LOG = FALSE)
    sprv <- prv * dBetaBinom(data_counts$ref_count, data_counts$ref_count + data_counts$var_count, pi = probshift, theta = theta_het, LOG = FALSE)
    
    pdata <- rowSums(cbind(spv, sprv, spr)) 
    
    if (any(pdata==0)){ 
      ProblemCases <- which(pdata==0)
      for(case in ProblemCases){
        var_part<-data_counts$var_count[case]
        ref_part<-data_counts$ref_count[case]
        spr_part<-dBetaBinom(ref_part, ref_part+var_part, pi = piRR, theta = theta_RR, LOG = TRUE)
        spv_part<-dBetaBinom(var_part, ref_part+var_part, pi = piVV, theta = theta_VV, LOG = TRUE)
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
      MomentEst <- MomentEst_MixedBetaBinom_SESlack(ref_counts = data_counts$ref_count, var_counts = data_counts$var_count, spr = spr, spv = spv, sprv = sprv, 
                                            pi_RR_fix = piRR, pi_VV_fix = piVV, pi_het_fix = probshift)
      if(sum(spr+spv)!=0){
        theta_RR <- MomentEst["theta_RR"]
        theta_VV <- MomentEst["theta_VV"]
      }
      theta_het <- MomentEst["theta_het"]
      if(is.nan(theta_het) | is.infinite(theta_het) | theta_het < ResetThetaMin){
        theta_het <- ResetThetaMin
      }
      if(is.nan(theta_RR) | is.infinite(theta_RR) | theta_RR < ResetThetaMin){
        theta_RR <- ResetThetaMin
      }
      if(is.nan(theta_VV) | is.infinite(theta_VV) | theta_VV < ResetThetaMin){
        theta_VV <- ResetThetaMin
      }
      
      if(sum(spr+spv)!=0){
        theta_RR_vec <- c()
        for(ii in c(10^(-5:-1), theta_RR) ){
          Likely <- -LogLikComp_RR(parvec = c(gtools::logit(ifelse(piRR >= 1-SE, 1-SE-eps, ifelse(piRR <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR)), max = 1-SE, min=1-SESlack*SE), log(ii)), ref_counts = data_counts$ref_count, 
                                   var_counts=data_counts$var_count, spr = spr,pim = 1-SESlack*SE, M = 1-SE)
          theta_RR_vec <- c(theta_RR_vec, Likely)
        }
        WhichOne <- which(theta_RR_vec==max(theta_RR_vec))
        theta_RR <- mean(c(10^(-5:-1), theta_RR)[WhichOne])
        
        theta_VV_vec <- c()
        for(ii in c(10^(-5:-1), theta_VV) ){
          Likely <- -LogLikComp_VV(parvec = c(gtools::logit(ifelse(piVV >= 1-SE, 1-SE-eps, ifelse(piVV <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV)), max = 1-SE, min=1-SESlack*SE), log(ii)), ref_counts = data_counts$ref_count, 
                                   var_counts=data_counts$var_count, spv = spv,pim = 1-SESlack*SE, M = 1-SE)
          theta_VV_vec <- c(theta_VV_vec, Likely)
        }
        WhichOne <- which(theta_VV_vec==max(theta_VV_vec))
        theta_VV <- mean(c(10^(-5:-1), theta_VV)[WhichOne])
      }
      theta_het_vec <- c()
      for(ii in c(10^(-5:-1), theta_het)){
        Likely <- -LogLikComp_het(parvec = c(gtools::logit(probshift), log(ii)), ref_counts=data_counts$ref_count, 
                                  var_counts=data_counts$var_count, sprv = sprv, piRR = piRR, piVV =piVV)
        theta_het_vec <- c(theta_het_vec, Likely)
      }
      WhichOne <- which(theta_het_vec==max(theta_het_vec))
      theta_het <- mean(c(10^(-5:-1), theta_het)[WhichOne])
      
    } else if(ReEstThetas == "simple" | nrep == 1){
      if(sum(spr+spv)!=0){
        theta_RR_vec <- c()
        for(ii in c(10^(-5:-1)) ){
          Likely <- -LogLikComp_RR(parvec = c(gtools::logit(ifelse(piRR >= 1-SE, 1-SE-eps, ifelse(piRR <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR)), max = 1-SE, min=1-SESlack*SE), log(ii)), ref_counts = data_counts$ref_count, 
                                   var_counts=data_counts$var_count, spr = spr,pim = 1-SESlack*SE, M = 1-SE)
          theta_RR_vec <- c(theta_RR_vec, Likely)
        }
        WhichOne <- which(theta_RR_vec==max(theta_RR_vec))
        theta_RR <- mean(c(10^(-5:-1))[WhichOne])
        
        theta_VV_vec <- c()
        for(ii in c(10^(-5:-1)) ){
          Likely <- -LogLikComp_VV(parvec = c(gtools::logit(ifelse(piVV >= 1-SE, 1-SE-eps, ifelse(piVV <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV)), max = 1-SE, min=1-SESlack*SE), log(ii)), ref_counts = data_counts$ref_count, 
                                   var_counts=data_counts$var_count, spv = spv,pim = 1-SESlack*SE, M = 1-SE)
          theta_VV_vec <- c(theta_VV_vec, Likely)
        }
        WhichOne <- which(theta_VV_vec==max(theta_VV_vec))
        theta_VV <- mean(c(10^(-5:-1))[WhichOne])
      }
      theta_het_vec <- c()
      for(ii in c(10^(-5:-1))){
        Likely <- -LogLikComp_het(parvec = c(gtools::logit(probshift), log(ii)), ref_counts=data_counts$ref_count, 
                                  var_counts=data_counts$var_count, sprv = sprv, piRR = piRR, piVV =piVV)
        theta_het_vec <- c(theta_het_vec, Likely)
      }
      WhichOne <- which(theta_het_vec==max(theta_het_vec))
      theta_het <- mean(c(10^(-5:-1))[WhichOne])
      
    }
    
    #compareDerivatives(LogLikComp_hom, GradComp_hom, t0 = log(0.01), SE = 0.002, ref_counts=ref_counts, var_counts=var_counts, spr = spr, spv = spv)
    
    if(ReEstThetas != "TryThree"){
      
      if(sum(spr+spv)!=0){
        if(NoSplitHom){
          
          theta_RR_clone <- theta_RR
          OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piRR >= 1-SE, 1-SE-eps, ifelse(piRR <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_RR, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_RR, gr = GradComp_RR, 
                                     method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr,pim = 1-SESlack*SE, M = 1-SE)},
                              error = function(e) NULL)
          if(is.null(OptObj)){
            OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piRR >= 1-SE, 1-SE-eps, ifelse(piRR <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_RR, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_RR, 
                                       method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr,pim = 1-SESlack*SE, M = 1-SE)},
                            error = function(e) NULL)
            if(is.null(OptObj)){
              OptObj <- optim(par = c(gtools::logit(ifelse(piRR >= 1-SE, 1-SE-eps, ifelse(piRR <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_RR, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_RR, 
                                         method = "Nelder-Mead", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr,pim = 1-SESlack*SE, M = 1-SE)
            }
            
          }
          
          piRR <- gtools::inv.logit(OptObj$par[1], max = 1-SE, min = 1-SESlack*SE)
          theta_RR <- exp(OptObj$par[2])
          
          if(theta_RR  > max(piRR, 1-piRR)){
            
            theta_RR <- min(c(theta_RR_clone, max((1-piRR)/10, piRR/10)) ) # Take some distance from the boundary at which bimodality occurs
            
            AObj <- tryCatch( {alabama::auglag(c(gtools::logit(ifelse(piRR >= 1-SE, 1-SE-eps, ifelse(piRR <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_RR, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_RR, gr = GradComp_RR, 
                                               hin = eval_g_f_RR, ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr,pim = 1-SESlack*SE, M = 1-SE, 
                                               control.outer = list("trace" = FALSE))}, error = function(e) NULL)
            if(is.null(AObj)){
              AObj <- alabama::auglag(c(gtools::logit(ifelse(piRR >= 1-SE, 1-SE-eps, ifelse(piRR <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_RR, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_RR, 
                                      hin = eval_g_f_RR, ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr,pim = 1-SESlack*SE, M = 1-SE, 
                                      control.outer = list("trace" = FALSE))
            }
            
            piRR <- gtools::inv.logit(AObj$par[1], max = 1-SE, min = 1-SESlack*SE)
            theta_RR <- exp(AObj$par[2])
            
          }
          
          theta_VV_clone <- theta_VV
          OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piVV >= 1-SE, 1-SE-eps, ifelse(piVV <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_VV, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_VV, gr = GradComp_VV, 
                                     method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv,pim = 1-SESlack*SE, M = 1-SE)},
                              error = function(e) NULL)
          if(is.null(OptObj)){
            OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piVV >= 1-SE, 1-SE-eps, ifelse(piVV <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_VV, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_VV, 
                            method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv,pim = 1-SESlack*SE, M = 1-SE)},
                            error = function(e) NULL)
            if(is.null(OptObj)){
              OptObj <- optim(par = c(gtools::logit(ifelse(piVV >= 1-SE, 1-SE-eps, ifelse(piVV <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_VV, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_VV, 
                              method = "Nelder-Mead", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv,pim = 1-SESlack*SE, M = 1-SE)
            }
          }
          
          piVV <- gtools::inv.logit(OptObj$par[1], max = 1-SE, min = 1-SESlack*SE)
          theta_VV <- exp(OptObj$par[2])
          
          if(theta_VV  > max(piVV, 1-piVV)){
            
            theta_VV <- min(c(theta_VV_clone, max((1-piVV)/10, piVV/10)) ) # Take some distance from the boundary at which bimodality occurs
            
            AObj <- tryCatch( {alabama::auglag(c(gtools::logit(ifelse(piVV >= 1-SE, 1-SE-eps, ifelse(piVV <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_VV, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_VV, gr = GradComp_VV, 
                                               hin = eval_g_f_VV, ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv,pim = 1-SESlack*SE, M = 1-SE, 
                                               control.outer = list("trace" = FALSE))}, error = function(e) NULL)
            if(is.null(AObj)){
              AObj <- alabama::auglag(c(gtools::logit(ifelse(piVV >= 1-SE, 1-SE-eps, ifelse(piVV <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_VV, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_VV, 
                                      hin = eval_g_f_VV, ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv,pim = 1-SESlack*SE, M = 1-SE, 
                                      control.outer = list("trace" = FALSE))
            }
            
            piVV <- gtools::inv.logit(AObj$par[1], max = 1-SE, min = 1-SESlack*SE)
            theta_VV <- exp(AObj$par[2])
            
          }
          
        } else{
          OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piRR >= 1-SE, 1-SE-eps, ifelse(piRR <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_RR, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_RR, gr = GradComp_RR, 
                                     method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr,pim = 1-SESlack*SE, M = 1-SE)},
                              error = function(e) NULL)
          if(is.null(OptObj)){
            OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piRR >= 1-SE, 1-SE-eps, ifelse(piRR <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_RR, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_RR, 
                            method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr,pim = 1-SESlack*SE, M = 1-SE)},
                            error = function(e) NULL)
            if(is.null(OptObj)){
              OptObj <- optim(par = c(gtools::logit(ifelse(piRR >= 1-SE, 1-SE-eps, ifelse(piRR <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_RR, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_RR, 
                              method = "Nelder-Mead", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr,pim = 1-SESlack*SE, M = 1-SE)
            }
          }
          
          piRR <- gtools::inv.logit(OptObj$par[1], max = 1-SE, min = 1-SESlack*SE)
          theta_RR <- exp(OptObj$par[2])
          
          OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piVV >= 1-SE, 1-SE-eps, ifelse(piVV <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_VV, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_VV, gr = GradComp_VV, 
                                     method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv, pim = 1-SESlack*SE, M = 1-SE)},
                              error = function(e) NULL)
          if(is.null(OptObj)){
            OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piVV >= 1-SE, 1-SE-eps, ifelse(piVV <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_VV, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_VV, 
                            method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv, pim = 1-SESlack*SE, M = 1-SE)},
                            error = function(e) NULL)
            if(is.null(OptObj)){
              OptObj <- optim(par = c(gtools::logit(ifelse(piVV >= 1-SE, 1-SE-eps, ifelse(piVV <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_VV, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_VV, 
                    method = "Nelder-Mead", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv, pim = 1-SESlack*SE, M = 1-SE)
            }
          }
          
          piVV <- gtools::inv.logit(OptObj$par[1], max = 1-SE, min = 1-SESlack*SE)
          theta_VV <- exp(OptObj$par[2])
        }
        
      }
      
      #compareDerivatives(LogLikComp_hom, GradComp_hom, HessComp_hom, t0 = 0.3, SE=SE, ref_counts=ref_counts, var_counts=var_counts, spr=spr, spv=spv)
      #compareDerivatives(LogLikComp_het, GradComp_het, t0 = c(gtools::logit(1-0.002), log(0.5)), ref_counts=ref_counts, var_counts=var_counts, sprv = sprv)
      
      if(NoSplitHet){
        theta_het_clone <- theta_het
        OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(probshift >= piRR, piRR-eps, ifelse(probshift <= 1-piVV, (1-piVV)+eps, probshift)), max = piRR, min = 1-piVV), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, gr = GradComp_het, method = "BFGS",
                                   ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, piRR = piRR, piVV = piVV)},
                            error = function(e) NULL)
        if(is.null(OptObj)){
          OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(probshift >= piRR, piRR-eps, ifelse(probshift <= 1-piVV, (1-piVV)+eps, probshift)), max = piRR, min = 1-piVV), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, method = "BFGS",
                          ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, piRR = piRR, piVV = piVV)},
                          error = function(e) NULL)
          if(is.null(OptObj)){
            OptObj <- optim(par = c(gtools::logit(ifelse(probshift >= piRR, piRR-eps, ifelse(probshift <= 1-piVV, (1-piVV)+eps, probshift)), max = piRR, min = 1-piVV), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, method = "Nelder-Mead",
                  ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, piRR = piRR, piVV = piVV)
          }
        }
        
        probshift <- gtools::inv.logit(OptObj$par[1], max = piRR, min = 1-piVV)
        theta_het <- exp(OptObj$par[2])
        
        if(theta_het  > max(probshift, 1-probshift)){
          
          theta_het <- min(c(theta_het_clone, max((1-probshift)/10, probshift/10)) ) # Take some distance from the boundary at which bimodality occurs
          
          AObj <- tryCatch( {alabama::auglag(c(gtools::logit(ifelse(probshift >= piRR, piRR-eps, ifelse(probshift <= 1-piVV, (1-piVV)+eps, probshift)), max = piRR, min = 1-piVV), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, gr = GradComp_het, hin = eval_g_f_2,
                                             ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, piRR = piRR, piVV = piVV, control.outer = list("trace" = FALSE))},
                            error = function(e) NULL)
          if(is.null(AObj)){
            AObj <- alabama::auglag(c(gtools::logit(ifelse(probshift >= piRR, piRR-eps, ifelse(probshift <= 1-piVV, (1-piVV)+eps, probshift)), max = piRR, min = 1-piVV), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, hin = eval_g_f_2,
                                    ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, piRR = piRR, piVV = piVV, control.outer = list("trace" = FALSE))
          }
          
          probshift <- gtools::inv.logit(AObj$par[1], max = piRR, min = 1-piVV)
          theta_het <- exp(AObj$par[2])
          
        }
      } else{
        OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(probshift >= piRR, piRR-eps, ifelse(probshift <= 1-piVV, (1-piVV)+eps, probshift)), max = piRR, min = 1-piVV), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, gr = GradComp_het, method = "BFGS",
                                   ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, piRR = piRR, piVV = piVV)},
                            error = function(e) NULL)
        if(is.null(OptObj)){
          OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(probshift >= piRR, piRR-eps, ifelse(probshift <= 1-piVV, (1-piVV)+eps, probshift)), max = piRR, min = 1-piVV), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, method = "BFGS",
                          ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, piRR = piRR, piVV = piVV)},
                          error = function(e) NULL)
          if(is.null(OptObj)){
            OptObj <- optim(par = c(gtools::logit(ifelse(probshift >= piRR, piRR-eps, ifelse(probshift <= 1-piVV, (1-piVV)+eps, probshift)), max = piRR, min = 1-piVV), log(min(max(theta_het, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_het, method = "Nelder-Mead",
                            ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, piRR = piRR, piVV = piVV)
          }
        }
        
        probshift <- gtools::inv.logit(OptObj$par[1], max = piRR, min = 1-piVV)
        theta_het <- exp(OptObj$par[2])
      }
      
      f <- function(thetaX) maelstRom::dBetaBinom(maxC,maxC,piRR,thetaX, LOG=TRUE)-maelstRom::dBetaBinom(maxC,maxC,probshift,theta_het, LOG=TRUE)
      US <- tryCatch( {uniroot(f, interval = c(0,1), tol=10^-16)}, error = function(e) NULL)
      if(!is.null(US)){
        theta_RR <- US$root
      }
      f <- function(thetaX) maelstRom::dBetaBinom(0,maxC,1-piVV,thetaX, LOG=TRUE)-maelstRom::dBetaBinom(0,maxC,probshift,theta_het, LOG=TRUE)
      US <- tryCatch( {uniroot(f, interval = c(0,1), tol=10^-16)}, error = function(e) NULL)
      if(!is.null(US)){
        theta_VV <- US$root
      }
      
    } else{
      
      
      
      if(sum(spr+spv)!=0){
        if(NoSplitHom){
          
          thetavec <- c()
          pivec <- c()
          likvec <- c()
          for(TH in thetaTRY){
          OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piRR >= 1-SE, 1-SE-eps, ifelse(piRR <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_RR, gr = GradComp_RR, 
                                     method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr, pim = 1-SESlack*SE, M = 1-SE)},
                              error = function(e) NULL)
          if(is.null(OptObj)){
            OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piRR >= 1-SE, 1-SE-eps, ifelse(piRR <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_RR, 
                            method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr, pim = 1-SESlack*SE, M = 1-SE)},
                            error = function(e) NULL)
            if(is.null(OptObj)){
              OptObj <- optim(par = c(gtools::logit(ifelse(piRR >= 1-SE, 1-SE-eps, ifelse(piRR <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_RR, 
                    method = "Nelder-Mead", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr, pim = 1-SESlack*SE, M = 1-SE)
            }
          }
          
          pivec <- c(pivec, gtools::inv.logit(OptObj$par[1], max = 1-SE, min = 1-SESlack*SE))
          thetavec <- c(thetavec, exp(OptObj$par[2]))
          likvec <- c(likvec, LogLikComp_RR(parvec = c(OptObj$par[1], OptObj$par[2]), ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr, pim = 1-SESlack*SE, M = 1-SE))
          }
          theta_RR_loc <- which(likvec == min(likvec))[1]
          piRR <- pivec[theta_RR_loc]
          theta_RR <- thetavec[theta_RR_loc]
          
          
          if(theta_RR  > max(piRR, 1-piRR)){
            
            thetavec <- c()
            pivec <- c()
            likvec <- c()
            for(TH in thetaTRY){
            AObj <- tryCatch( {alabama::auglag(c(gtools::logit(ifelse(piRR >= 1-SE, 1-SE-eps, ifelse(piRR <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_RR, gr = GradComp_RR, 
                                               hin = eval_g_f_RR, ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr, pim = 1-SESlack*SE, M = 1-SE, 
                                               control.outer = list("trace" = FALSE))}, error = function(e) NULL)
            if(is.null(AObj)){
              AObj <- alabama::auglag(c(gtools::logit(ifelse(piRR >= 1-SE, 1-SE-eps, ifelse(piRR <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_RR, 
                                      hin = eval_g_f_RR, ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr, pim = 1-SESlack*SE, M = 1-SE, 
                                      control.outer = list("trace" = FALSE))
            }
            
            pivec <- c(pivec, gtools::inv.logit(AObj$par[1], max = 1-SE, min = 1-SESlack*SE))
            thetavec <- c(thetavec, exp(AObj$par[2]))
            likvec <- c(likvec, LogLikComp_RR(parvec = c(AObj$par[1], AObj$par[2]), ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr, pim = 1-SESlack*SE, M = 1-SE))
            }
            theta_RR_loc <- which(likvec == min(likvec))[1]
            piRR <- pivec[theta_RR_loc]
            theta_RR <- thetavec[theta_RR_loc]
          }
          
          
          
          
          thetavec <- c()
          pivec <- c()
          likvec <- c()
          for(TH in thetaTRY){
            OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piVV >= 1-SE, 1-SE-eps, ifelse(piVV <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_VV, gr = GradComp_VV, 
                                       method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv, pim = 1-SESlack*SE, M = 1-SE)},
                                error = function(e) NULL)
            if(is.null(OptObj)){
              OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piVV >= 1-SE, 1-SE-eps, ifelse(piVV <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_VV, 
                              method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv, pim = 1-SESlack*SE, M = 1-SE)},
                              error = function(e) NULL)
              if(is.null(OptObj)){
                OptObj <- optim(par = c(gtools::logit(ifelse(piVV >= 1-SE, 1-SE-eps, ifelse(piVV <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_VV, 
                      method = "Nelder-Mead", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv, pim = 1-SESlack*SE, M = 1-SE)
              }
            }
            
            pivec <- c(pivec, gtools::inv.logit(OptObj$par[1], max = 1-SE, min = 1-SESlack*SE))
            thetavec <- c(thetavec, exp(OptObj$par[2]))
            likvec <- c(likvec, LogLikComp_VV(parvec = c(OptObj$par[1], OptObj$par[2]), ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv, pim = 1-SESlack*SE, M = 1-SE))
          }
          theta_VV_loc <- which(likvec == min(likvec))[1]
          piVV <- pivec[theta_VV_loc]
          theta_VV <- thetavec[theta_VV_loc]
          
          
          if(theta_VV  > max(piVV, 1-piVV)){
            
            thetavec <- c()
            pivec <- c()
            likvec <- c()
            for(TH in thetaTRY){
              AObj <- tryCatch( {alabama::auglag(c(gtools::logit(ifelse(piVV >= 1-SE, 1-SE-eps, ifelse(piVV <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_VV, gr = GradComp_VV, 
                                                 hin = eval_g_f_VV, ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv, pim = 1-SESlack*SE, M = 1-SE, 
                                                 control.outer = list("trace" = FALSE))}, error = function(e) NULL)
              if(is.null(AObj)){
                AObj <- alabama::auglag(c(gtools::logit(ifelse(piVV >= 1-SE, 1-SE-eps, ifelse(piVV <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_VV, 
                                        hin = eval_g_f_VV, ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv, pim = 1-SESlack*SE, M = 1-SE, 
                                        control.outer = list("trace" = FALSE))
              }
              
              pivec <- c(pivec, gtools::inv.logit(AObj$par[1], max = 1-SE, min = 1-SESlack*SE))
              thetavec <- c(thetavec, exp(AObj$par[2]))
              likvec <- c(likvec, LogLikComp_VV(parvec = c(AObj$par[1], AObj$par[2]), ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv, pim = 1-SESlack*SE, M = 1-SE))
            }
            theta_VV_loc <- which(likvec == min(likvec))[1]
            piVV <- pivec[theta_VV_loc]
            theta_VV <- thetavec[theta_VV_loc]
          }
          
          
          
          
        } else{
          
          thetavec <- c()
          pivec <- c()
          likvec <- c()
          for(TH in thetaTRY){
            OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piRR >= 1-SE, 1-SE-eps, ifelse(piRR <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_RR, gr = GradComp_RR, 
                                       method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr, pim = 1-SESlack*SE, M = 1-SE)},
                                error = function(e) NULL)
            if(is.null(OptObj)){
              OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piRR >= 1-SE, 1-SE-eps, ifelse(piRR <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_RR, 
                              method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr, pim = 1-SESlack*SE, M = 1-SE)},
                              error = function(e) NULL)
              if(is.null(OptObj)){
                OptObj <- optim(par = c(gtools::logit(ifelse(piRR >= 1-SE, 1-SE-eps, ifelse(piRR <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_RR, 
                      method = "Nelder-Mead", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr, pim = 1-SESlack*SE, M = 1-SE)
              }
            }
            
            pivec <- c(pivec, gtools::inv.logit(OptObj$par[1], max = 1-SE, min = 1-SESlack*SE))
            thetavec <- c(thetavec, exp(OptObj$par[2]))
            likvec <- c(likvec, LogLikComp_RR(parvec = c(OptObj$par[1], OptObj$par[2]), ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr, pim = 1-SESlack*SE, M = 1-SE))
          }
          theta_RR_loc <- which(likvec == min(likvec))[1]
          piRR <- pivec[theta_RR_loc]
          theta_RR <- thetavec[theta_RR_loc]
          
          thetavec <- c()
          pivec <- c()
          likvec <- c()
          for(TH in thetaTRY){
            OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piVV >= 1-SE, 1-SE-eps, ifelse(piVV <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_VV, gr = GradComp_VV, 
                                       method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv, pim = 1-SESlack*SE, M = 1-SE)},
                                error = function(e) NULL)
            if(is.null(OptObj)){
              OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piVV >= 1-SE, 1-SE-eps, ifelse(piVV <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_VV, 
                              method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv, pim = 1-SESlack*SE, M = 1-SE)},
                              error = function(e) NULL)
              if(is.null(OptObj)){
                OptObj <- optim(par = c(gtools::logit(ifelse(piVV >= 1-SE, 1-SE-eps, ifelse(piVV <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_VV, 
                      method = "Nelder-Mead", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv, pim = 1-SESlack*SE, M = 1-SE)
              }
            }
            
            pivec <- c(pivec, gtools::inv.logit(OptObj$par[1], max = 1-SE, min = 1-SESlack*SE))
            thetavec <- c(thetavec, exp(OptObj$par[2]))
            likvec <- c(likvec, LogLikComp_VV(parvec = c(OptObj$par[1], OptObj$par[2]), ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv, pim = 1-SESlack*SE, M = 1-SE))
          }
          theta_VV_loc <- which(likvec == min(likvec))[1]
          piVV <- pivec[theta_VV_loc]
          theta_VV <- thetavec[theta_VV_loc]
          
          
        }
        
      }
      
      #compareDerivatives(LogLikComp_hom, GradComp_hom, HessComp_hom, t0 = 0.3, SE=SE, ref_counts=ref_counts, var_counts=var_counts, spr=spr, spv=spv)
      #compareDerivatives(LogLikComp_het, GradComp_het, t0 = c(gtools::logit(1-0.002), log(0.5)), ref_counts=ref_counts, var_counts=var_counts, sprv = sprv)
      
      if(NoSplitHet){
        thetavec <- c()
        pivec <- c()
        likvec <- c()
        for(TH in thetaTRY){
        OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(probshift >= piRR, piRR-eps, ifelse(probshift <= 1-piVV, (1-piVV)+eps, probshift)), max = piRR, min = 1-piVV), log(TH)), fn = LogLikComp_het, gr = GradComp_het, method = "BFGS",
                                   ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, piRR = piRR, piVV = piVV)},
                            error = function(e) NULL)
        if(is.null(OptObj)){
          OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(probshift >= piRR, piRR-eps, ifelse(probshift <= 1-piVV, (1-piVV)+eps, probshift)), max = piRR, min = 1-piVV), log(TH)), fn = LogLikComp_het, method = "BFGS",
                          ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, piRR = piRR, piVV = piVV)},
                          error = function(e) NULL)
          if(is.null(OptObj)){
            OptObj <- optim(par = c(gtools::logit(ifelse(probshift >= piRR, piRR-eps, ifelse(probshift <= 1-piVV, (1-piVV)+eps, probshift)), max = piRR, min = 1-piVV), log(TH)), fn = LogLikComp_het, 
                            method = "Nelder-Mead", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, piRR = piRR, piVV = piVV)
          }
        }
        
        pivec <- c(pivec, gtools::inv.logit(OptObj$par[1], max = piRR, min = 1-piVV))
        thetavec <- c(thetavec, exp(OptObj$par[2]))
        likvec <- c(likvec, LogLikComp_het(c(OptObj$par[1], OptObj$par[2]), ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, sprv=sprv, piRR = piRR, piVV = piVV))
        }
        
        theta_het_loc <- which(likvec == min(likvec))[1]
        theta_het <- thetavec[theta_het_loc]
        probshift <- pivec[theta_het_loc]
        
        if(theta_het  > max(probshift, 1-probshift)){
          thetavec <- c()
          pivec <- c()
          likvec <- c()
          for(TH in thetaTRY){
          AObj <- tryCatch( {alabama::auglag(c(gtools::logit(ifelse(probshift >= piRR, piRR-eps, ifelse(probshift <= 1-piVV, (1-piVV)+eps, probshift)), max = piRR, min = 1-piVV), log(TH)), fn = LogLikComp_het, gr = GradComp_het, hin = eval_g_f_2,
                                             ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, piRR = piRR, piVV = piVV, control.outer = list("trace" = FALSE))},
                            error = function(e) NULL)
          if(is.null(AObj)){
            AObj <- alabama::auglag(c(gtools::logit(ifelse(probshift >= piRR, piRR-eps, ifelse(probshift <= 1-piVV, (1-piVV)+eps, probshift)), max = piRR, min = 1-piVV), log(TH)), fn = LogLikComp_het, hin = eval_g_f_2,
                                    ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, piRR = piRR, piVV = piVV, control.outer = list("trace" = FALSE))
          }
          
          pivec <- c(pivec, gtools::inv.logit(AObj$par[1], max = piRR, min = 1-piVV))
          thetavec <- c(thetavec, exp(AObj$par[2]))
          likvec <- c(likvec, LogLikComp_het(c(AObj$par[1], AObj$par[2]), ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, sprv=sprv, piRR = piRR, piVV = piVV))
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
          OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(probshift >= piRR, piRR-eps, ifelse(probshift <= 1-piVV, (1-piVV)+eps, probshift)), max = piRR, min = 1-piVV), log(TH)), fn = LogLikComp_het, gr = GradComp_het, method = "BFGS",
                                     ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, piRR = piRR, piVV = piVV)},
                              error = function(e) NULL)
          if(is.null(OptObj)){
            OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(probshift >= piRR, piRR-eps, ifelse(probshift <= 1-piVV, (1-piVV)+eps, probshift)), max = piRR, min = 1-piVV), log(TH)), fn = LogLikComp_het, method = "BFGS",
                            ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, piRR = piRR, piVV = piVV)},
                            error = function(e) NULL)
            if(is.null(OptObj)){
              OptObj <- optim(par = c(gtools::logit(ifelse(probshift >= piRR, piRR-eps, ifelse(probshift <= 1-piVV, (1-piVV)+eps, probshift)), max = piRR, min = 1-piVV), log(TH)), fn = LogLikComp_het,
                              method = "Nelder-Mead", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, piRR = piRR, piVV = piVV)
            }
            
          }
          
          pivec <- c(pivec, gtools::inv.logit(OptObj$par[1], max = piRR, min = 1-piVV))
          thetavec <- c(thetavec, exp(OptObj$par[2]))
          likvec <- c(likvec, LogLikComp_het(c(OptObj$par[1], OptObj$par[2]), ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, sprv=sprv, piRR = piRR, piVV = piVV))
        }
        
        theta_het_loc <- which(likvec == min(likvec))[1]
        theta_het <- thetavec[theta_het_loc]
        probshift <- pivec[theta_het_loc]
        
      }
      
      f <- function(thetaX) maelstRom::dBetaBinom(maxC,maxC,piRR,thetaX, LOG=TRUE)-maelstRom::dBetaBinom(maxC,maxC,probshift,theta_het, LOG=TRUE)
      US <- tryCatch( {uniroot(f, interval = c(0,1), tol=10^-16)}, error = function(e) NULL)
      if(!is.null(US)){
        theta_RR <- US$root
      }
      f <- function(thetaX) maelstRom::dBetaBinom(0,maxC,1-piVV,thetaX, LOG=TRUE)-maelstRom::dBetaBinom(0,maxC,probshift,theta_het, LOG=TRUE)
      US <- tryCatch( {uniroot(f, interval = c(0,1), tol=10^-16)}, error = function(e) NULL)
      if(!is.null(US)){
        theta_VV <- US$root
      }
      
    }
    
    #if(probshift > (1-10*SE)){
    #  probshift <- (1-10*SE)
    #} else if(probshift < 10*SE){
    #  probshift <- 10*SE
    #}
    
    Q <- sum(ifelse(pr>0, spr*log(pr), 0) + spr*dBetaBinom(data_counts$ref_count, data_counts$ref_count+data_counts$var_count, pi = piRR, theta = theta_RR, LOG = TRUE) +
               ifelse(prv>0, sprv*log(prv), 0) + sprv*dBetaBinom(data_counts$ref_count, data_counts$ref_count+data_counts$var_count, pi = probshift, theta = theta_het, LOG = TRUE) +
               ifelse(pv>0, spv*log(pv), 0) + spv*dBetaBinom(data_counts$var_count, data_counts$ref_count+data_counts$var_count, pi = piVV, theta = theta_VV, LOG = TRUE))
    
    dlta <- abs(Qold-Q)
    
    allelefreq <- mean(spr) + mean(sprv)/2
    data_counts$allelefreq <- allelefreq
  }
  
  #############################################################################################################################################################
  #############################################################################################################################################################
  
  data_counts$pvv <- spv; data_counts$prr <- spr; data_counts$prv <- sprv
  
  spr_NS <- dBetaBinom(data_counts$ref_count, data_counts$ref_count + data_counts$var_count, pi = piRR, theta = theta_RR, LOG = FALSE)
  spv_NS <- dBetaBinom(data_counts$var_count, data_counts$ref_count + data_counts$var_count, pi = piVV, theta = theta_VV, LOG = FALSE)
  sprv_NS <- dBetaBinom(data_counts$ref_count, data_counts$ref_count + data_counts$var_count, pi = probshift, theta = theta_het, LOG = FALSE)
  pdata_NS <- rowSums(cbind(spv_NS, sprv_NS, spr_NS)) 
  if (any(pdata_NS==0)){ 
    ProblemCases <- which(pdata_NS==0)
    for(case in ProblemCases){
      var_part<-data_counts$var_count[case]
      ref_part<-data_counts$ref_count[case]
      spr_NS_part<-dBetaBinom(ref_part, ref_part+var_part, pi = piRR, theta = theta_RR, LOG = TRUE)
      spv_NS_part<-dBetaBinom(var_part, ref_part+var_part, pi = piVV, theta = theta_VV, LOG = TRUE)
      sprv_NS_part<-dBetaBinom(ref_part, ref_part+var_part, pi = probshift, theta = theta_het, LOG = TRUE)
      spvec<-c(spr_NS_part, sprv_NS_part, spv_NS_part)
      if(spr_NS_part==max(spvec)){
        spr_NS[case]<-1
        pdata_NS[case]<-1
      } else if(sprv_NS_part==max(spvec)){
        sprv_NS[case]<-1
        pdata_NS[case]<-1
      } else{
        spv_NS[case]<-1
        pdata_NS[case]<-1
      }
    }
  }
  spv_NS <- spv_NS/pdata_NS
  sprv_NS <- sprv_NS/pdata_NS
  spr_NS <- spr_NS/pdata_NS
  
  data_counts$genotypeN <- c("rr", "rv", "vv")[apply(data_counts[, c("prr", "prv", "pvv")], 1, function(x) which(x==max(x))[1])]
  data_counts$genotypeN_NS <- c("rr", "rv", "vv")[apply(cbind(spr_NS, sprv_NS, spv_NS), 1, function(x) which(x==max(x))[1])]
  
  if(!fitH0){
    results <- list(AB = probshift, piRR = piRR, piVV = piVV, nrep = nrep, quality = quality, rho_rr = pr, rho_vv = pv, rho_rv = prv, data_hash = data_counts, 
                    theta_RR = theta_RR, theta_VV = theta_VV, theta_het = theta_het)
    return(results)
  }
  
  # A whole bunch of potential GOF-measures
  data_counts$FitRR<-0
  data_counts$FitRV<-0
  data_counts$FitVV<-0
  data_counts$ExactPVAL<-0
  data_counts$TotPVAL<-0
  for(DCR in 1:nrow(data_counts)){
    data_counts$FitRR[DCR] <- (BetaBinom_test_pvalTS(m = data_counts$ref_count[DCR], n = data_counts$ref_count[DCR] + data_counts$var_count[DCR], pi = piRR, theta = theta_RR))
    data_counts$FitRV[DCR] <- (BetaBinom_test_pvalTS(m = data_counts$ref_count[DCR], n = data_counts$ref_count[DCR] + data_counts$var_count[DCR], pi = probshift, theta = theta_het))
    data_counts$FitVV[DCR] <- (BetaBinom_test_pvalTS(m = data_counts$var_count[DCR], n = data_counts$ref_count[DCR] + data_counts$var_count[DCR], pi = piVV, theta = theta_VV))
    data_counts$ExactPVAL[DCR] <- data_counts$prr[DCR]*data_counts$FitRR[DCR] +
      data_counts$prv[DCR]*data_counts$FitRV[DCR] +
      data_counts$pvv[DCR]*data_counts$FitVV[DCR]
    # A TOTAL exact test
    TotC <- data_counts$ref_count[DCR] + data_counts$var_count[DCR]
    CountCombos <- cbind(0:TotC, TotC:0)
    
    ComboProbs <- pr * dBetaBinom(CountCombos[,1], rep(TotC, TotC+1), pi = piRR, theta = theta_RR, LOG = FALSE) + pv * dBetaBinom(CountCombos[,2], rep(TotC, TotC+1), pi = piVV, theta = theta_VV, LOG = FALSE) + prv * dBetaBinom(CountCombos[,1], rep(TotC, TotC+1), pi = probshift, theta = theta_het, LOG = FALSE)
    
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
  dmixase_corrected <- pmf_betabinomMix_SESlack(data_counts$ref_count, data_counts$var_count, probshift, piRR, piVV, mean(spr), mean(spv), mean(sprv), 
                                                  theta_RR = theta_RR, theta_VV = theta_VV, theta_het = theta_het) * (data_counts$ref_count + data_counts$var_count + 1)
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
    spr_H0 <- pr_H0 * dBetaBinom(data_counts$ref_count, data_counts$ref_count + data_counts$var_count, pi = piRR_H0, theta = theta_RR_H0, LOG = FALSE)
    spv_H0 <- pv_H0 * dBetaBinom(data_counts$var_count, data_counts$ref_count + data_counts$var_count, pi = piVV_H0, theta = theta_VV_H0, LOG = FALSE)
    sprv_H0 <- prv_H0 * dBetaBinom(data_counts$ref_count, data_counts$ref_count + data_counts$var_count, pi = 0.5, theta = theta_het_H0, LOG = FALSE)
    
    pdata <- rowSums(cbind(spr_H0, spv_H0, sprv_H0)) 
    
    if (any(pdata==0)){ 
      ProblemCases <- which(pdata==0)
      for(case in ProblemCases){
        var_part<-data_counts$var_count[case]
        ref_part<-data_counts$ref_count[case]
        spr_part<-dBetaBinom(ref_part, ref_part+var_part, pi = piRR_H0, theta = theta_RR_H0, LOG = TRUE)
        spv_part<-dBetaBinom(var_part, ref_part+var_part, pi = piVV_H0, theta = theta_VV_H0, LOG = TRUE)
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
      MomentEst <- MomentEst_MixedBetaBinom_SESlack(ref_counts = data_counts$ref_count, var_counts = data_counts$var_count, spr = spr_H0, spv = spv_H0, sprv = sprv_H0, 
                                                    pi_RR_fix = piRR_H0, pi_VV_fix = piVV_H0, pi_het_fix = 0.5)
      if(sum(spr+spv)!=0){
        theta_RR_H0 <- MomentEst["theta_RR"]
        theta_VV_H0 <- MomentEst["theta_VV"]
      }
      theta_het_H0 <- MomentEst["theta_het"]
      if(is.nan(theta_het_H0) | is.infinite(theta_het_H0) | theta_het_H0 < ResetThetaMin){
        theta_het_H0 <- ResetThetaMin
      }
      if(is.nan(theta_RR_H0) | is.infinite(theta_RR_H0) | theta_RR_H0 < ResetThetaMin){
        theta_RR_H0 <- ResetThetaMin
      }
      if(is.nan(theta_VV_H0) | is.infinite(theta_VV_H0) | theta_VV_H0 < ResetThetaMin){
        theta_VV_H0 <- ResetThetaMin
      }
      
      if(sum(spr_H0+spv_H0)!=0){
        theta_RR_H0_vec <- c()
        for(ii in c(10^(-5:-1), theta_RR_H0) ){
          Likely <- -LogLikComp_RR(parvec = c(gtools::logit(ifelse(piRR_H0 >= 1-SE, 1-SE-eps, ifelse(piRR_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR_H0)), max = 1-SE, min=1-SESlack*SE), log(ii)), ref_counts = data_counts$ref_count, 
                                   var_counts=data_counts$var_count, spr = spr_H0, pim = 1-SESlack*SE, M = 1-SE)
          theta_RR_H0_vec <- c(theta_RR_H0_vec, Likely)
        }
        WhichOne <- which(theta_RR_H0_vec==max(theta_RR_H0_vec))
        theta_RR_H0 <- mean(c(10^(-5:-1), theta_RR_H0)[WhichOne])
        
        theta_VV_H0_vec <- c()
        for(ii in c(10^(-5:-1), theta_VV_H0) ){
          Likely <- -LogLikComp_VV(parvec = c(gtools::logit(ifelse(piVV_H0 >= 1-SE, 1-SE-eps, ifelse(piVV_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV_H0)), max = 1-SE, min=1-SESlack*SE), log(ii)), ref_counts = data_counts$ref_count, 
                                   var_counts=data_counts$var_count, spv = spv_H0, pim = 1-SESlack*SE, M = 1-SE)
          theta_VV_H0_vec <- c(theta_VV_H0_vec, Likely)
        }
        WhichOne <- which(theta_VV_H0_vec==max(theta_VV_H0_vec))
        theta_VV_H0 <- mean(c(10^(-5:-1), theta_VV_H0)[WhichOne])
      }
      theta_het_H0_vec <- c()
      for(ii in c(10^(-5:-1), theta_het_H0)){
        Likely <- -LogLikComp_het_H0(theta = log(ii), probshift = 0.5, ref_counts=data_counts$ref_count, 
                                  var_counts=data_counts$var_count, sprv = sprv_H0)
        theta_het_H0_vec <- c(theta_het_H0_vec, Likely)
      }
      WhichOne <- which(theta_het_H0_vec==max(theta_het_H0_vec))
      theta_het_H0 <- mean(c(10^(-5:-1), theta_het_H0)[WhichOne])
      
    } else if(ReEstThetas == "simple" | nrep == 1){
      if(sum(spr_H0+spv_H0)!=0){
        theta_RR_H0_vec <- c()
        for(ii in c(10^(-5:-1)) ){
          Likely <- -LogLikComp_RR(parvec = c(gtools::logit(ifelse(piRR_H0 >= 1-SE, 1-SE-eps, ifelse(piRR_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR_H0)), max = 1-SE, min=1-SESlack*SE), log(ii)), ref_counts = data_counts$ref_count, 
                                   var_counts=data_counts$var_count, spr = spr_H0, pim = 1-SESlack*SE, M = 1-SE)
          theta_RR_H0_vec <- c(theta_RR_H0_vec, Likely)
        }
        WhichOne <- which(theta_RR_H0_vec==max(theta_RR_H0_vec))
        theta_RR_H0 <- mean(c(10^(-5:-1))[WhichOne])
        
        theta_VV_H0_vec <- c()
        for(ii in c(10^(-5:-1)) ){
          Likely <- -LogLikComp_VV(parvec = c(gtools::logit(ifelse(piVV_H0 >= 1-SE, 1-SE-eps, ifelse(piVV_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV_H0)), max = 1-SE, min=1-SESlack*SE), log(ii)), ref_counts = data_counts$ref_count, 
                                   var_counts=data_counts$var_count, spv = spv_H0, pim = 1-SESlack*SE, M = 1-SE)
          theta_VV_H0_vec <- c(theta_VV_H0_vec, Likely)
        }
        WhichOne <- which(theta_VV_H0_vec==max(theta_VV_H0_vec))
        theta_VV_H0 <- mean(c(10^(-5:-1))[WhichOne])
      }
      theta_het_H0_vec <- c()
      for(ii in c(10^(-5:-1))){
        Likely <- -LogLikComp_het_H0(theta = log(ii), probshift = 0.5, ref_counts=data_counts$ref_count, 
                                  var_counts=data_counts$var_count, sprv = sprv_H0)
        theta_het_H0_vec <- c(theta_het_H0_vec, Likely)
      }
      WhichOne <- which(theta_het_H0_vec==max(theta_het_H0_vec))
      theta_het_H0 <- mean(c(10^(-5:-1))[WhichOne])
      
    }
    
    #compareDerivatives(LogLikComp_hom, GradComp_hom, t0 = log(0.01), SE = 0.002, ref_counts=ref_counts, var_counts=var_counts, spr = spr, spv = spv)
    
    if(ReEstThetas != "TryThree"){
      
      if(sum(spr+spv)!=0){
        if(NoSplitHom){
          
          theta_RR_H0_clone <- theta_RR_H0
          OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piRR_H0 >= 1-SE, 1-SE-eps, ifelse(piRR_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR_H0)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_RR_H0, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_RR, gr = GradComp_RR, 
                                     method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, pim = 1-SESlack*SE, M = 1-SE)},
                              error = function(e) NULL)
          if(is.null(OptObj)){
            OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piRR_H0 >= 1-SE, 1-SE-eps, ifelse(piRR_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR_H0)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_RR_H0, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_RR, 
                            method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, pim = 1-SESlack*SE, M = 1-SE)},
                            error = function(e) NULL)
            if(is.null(OptObj)){
              OptObj <- optim(par = c(gtools::logit(ifelse(piRR_H0 >= 1-SE, 1-SE-eps, ifelse(piRR_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR_H0)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_RR_H0, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_RR, 
                              method = "Nelder-Mead", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, pim = 1-SESlack*SE, M = 1-SE)
            }
          }
          
          piRR_H0 <- gtools::inv.logit(OptObj$par[1], max = 1-SE, min = 1-SESlack*SE)
          theta_RR_H0 <- exp(OptObj$par[2])
          
          if(theta_RR_H0  > max(piRR_H0, 1-piRR_H0)){
            
            theta_RR_H0 <- min(c(theta_RR_H0_clone, max((1-piRR_H0)/10, piRR_H0/10)) ) # Take some distance from the boundary at which bimodality occurs
            
            AObj <- tryCatch( {alabama::auglag(c(gtools::logit(ifelse(piRR_H0 >= 1-SE, 1-SE-eps, ifelse(piRR_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR_H0)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_RR_H0, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_RR, gr = GradComp_RR, 
                                               hin = eval_g_f_RR, ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, pim = 1-SESlack*SE, M = 1-SE, 
                                               control.outer = list("trace" = FALSE))}, error = function(e) NULL)
            if(is.null(AObj)){
              AObj <- alabama::auglag(c(gtools::logit(ifelse(piRR_H0 >= 1-SE, 1-SE-eps, ifelse(piRR_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR_H0)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_RR_H0, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_RR, 
                                      hin = eval_g_f_RR, ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, pim = 1-SESlack*SE, M = 1-SE, 
                                      control.outer = list("trace" = FALSE))
            }
            
            piRR_H0 <- gtools::inv.logit(AObj$par[1], max = 1-SE, min = 1-SESlack*SE)
            theta_RR_H0 <- exp(AObj$par[2])
            
          }
          
          theta_VV_H0_clone <- theta_VV_H0
          OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piVV_H0 >= 1-SE, 1-SE-eps, ifelse(piVV_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV_H0)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_VV_H0, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_VV, gr = GradComp_VV, 
                                     method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv_H0, pim = 1-SESlack*SE, M = 1-SE)},
                              error = function(e) NULL)
          if(is.null(OptObj)){
            OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piVV_H0 >= 1-SE, 1-SE-eps, ifelse(piVV_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV_H0)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_VV_H0, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_VV, 
                            method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv_H0, pim = 1-SESlack*SE, M = 1-SE)},
                            error = function(e) NULL)
            if(is.null(OptObj)){
              OptObj <- optim(par = c(gtools::logit(ifelse(piVV_H0 >= 1-SE, 1-SE-eps, ifelse(piVV_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV_H0)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_VV_H0, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_VV, 
                              method = "Nelder-Mead", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv_H0, pim = 1-SESlack*SE, M = 1-SE)
            }
          }
          
          piVV_H0 <- gtools::inv.logit(OptObj$par[1], max = 1-SE, min = 1-SESlack*SE)
          theta_VV_H0 <- exp(OptObj$par[2])
          
          if(theta_VV_H0  > max(piVV_H0, 1-piVV_H0)){
            
            theta_VV_H0 <- min(c(theta_VV_H0_clone, max((1-piVV_H0)/10, piVV_H0/10)) ) # Take some distance from the boundary at which bimodality occurs
            
            AObj <- tryCatch( {alabama::auglag(c(gtools::logit(ifelse(piVV_H0 >= 1-SE, 1-SE-eps, ifelse(piVV_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV_H0)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_VV_H0, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_VV, gr = GradComp_VV, 
                                               hin = eval_g_f_VV, ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv_H0, pim = 1-SESlack*SE, M = 1-SE, 
                                               control.outer = list("trace" = FALSE))}, error = function(e) NULL)
            if(is.null(AObj)){
              AObj <- alabama::auglag(c(gtools::logit(ifelse(piVV_H0 >= 1-SE, 1-SE-eps, ifelse(piVV_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV_H0)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_VV_H0, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_VV, 
                                      hin = eval_g_f_VV, ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv_H0, pim = 1-SESlack*SE, M = 1-SE, 
                                      control.outer = list("trace" = FALSE))
            }
            
            piVV_H0 <- gtools::inv.logit(AObj$par[1], max = 1-SE, min = 1-SESlack*SE)
            theta_VV_H0 <- exp(AObj$par[2])
            
          }
          
        } else{
          OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piRR_H0 >= 1-SE, 1-SE-eps, ifelse(piRR_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR_H0)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_RR_H0, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_RR, gr = GradComp_RR, 
                                     method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, pim = 1-SESlack*SE, M = 1-SE)},
                              error = function(e) NULL)
          if(is.null(OptObj)){
            OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piRR_H0 >= 1-SE, 1-SE-eps, ifelse(piRR_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR_H0)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_RR_H0, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_RR, 
                            method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, pim = 1-SESlack*SE, M = 1-SE)},
                            error = function(e) NULL)
            if(is.null(OptObj)){
              OptObj <- optim(par = c(gtools::logit(ifelse(piRR_H0 >= 1-SE, 1-SE-eps, ifelse(piRR_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR_H0)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_RR_H0, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_RR, 
                              method = "Nelder-Mead", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, pim = 1-SESlack*SE, M = 1-SE)
            }
          }
          
          piRR_H0 <- gtools::inv.logit(OptObj$par[1], max = 1-SE, min = 1-SESlack*SE)
          theta_RR_H0 <- exp(OptObj$par[2])
          
          OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piVV_H0 >= 1-SE, 1-SE-eps, ifelse(piVV_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV_H0)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_VV_H0, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_VV, gr = GradComp_VV, 
                                     method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv_H0, pim = 1-SESlack*SE, M = 1-SE)},
                              error = function(e) NULL)
          if(is.null(OptObj)){
            OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piVV_H0 >= 1-SE, 1-SE-eps, ifelse(piVV_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV_H0)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_VV_H0, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_VV, 
                            method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv_H0, pim = 1-SESlack*SE, M = 1-SE)},
                            error = function(e) NULL)
            if(is.null(OptObj)){
              OptObj <- optim(par = c(gtools::logit(ifelse(piVV_H0 >= 1-SE, 1-SE-eps, ifelse(piVV_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV_H0)), max = 1-SE, min = 1-SESlack*SE), log(min(max(theta_VV_H0, ResetThetaMin), ResetThetaMax))), fn = LogLikComp_VV, 
                              method = "Nelder-Mead", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv_H0, pim = 1-SESlack*SE, M = 1-SE)
            }
          }
          
          piVV_H0 <- gtools::inv.logit(OptObj$par[1], max = 1-SE, min = 1-SESlack*SE)
          theta_VV_H0 <- exp(OptObj$par[2])
        }
        
      }
      
      #compareDerivatives(LogLikComp_hom, GradComp_hom, HessComp_hom, t0 = 0.3, SE=SE, ref_counts=ref_counts, var_counts=var_counts, spr=spr, spv=spv)
      #compareDerivatives(LogLikComp_het, GradComp_het, t0 = c(gtools::logit(1-0.002), log(0.5)), ref_counts=ref_counts, var_counts=var_counts, sprv = sprv)
      
      if(NoSplitHet){
        theta_het_H0_clone <- theta_het_H0
        OptObj <- tryCatch( {optim(par = log(min(max(theta_het_H0, ResetThetaMin), ResetThetaMax)), probshift=0.5, fn = LogLikComp_het_H0, gr = GradComp_het_H0, method = "BFGS",
                                   ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0)},
                            error = function(e) NULL)
        if(is.null(OptObj)){
          OptObj <- optim(par = log(min(max(theta_het_H0, ResetThetaMin), ResetThetaMax)), probshift=0.5, fn = LogLikComp_het_H0, method = "BFGS",
                          ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0)
        }
        
        theta_het_H0 <- exp(OptObj$par[1])
        
        if(theta_het_H0  > max(0.5, 1-0.5)){
          
          theta_het_H0 <- min(c(theta_het_H0_clone, max((1-0.5)/10, 0.5/10)) ) # Take some distance from the boundary at which bimodality occurs
          
          AObj <- tryCatch( {alabama::auglag(log(min(max(theta_het_H0, ResetThetaMin), ResetThetaMax)), probshift=0.5, fn = LogLikComp_het_H0, gr = GradComp_het_H0, hin = eval_g_f_2_H0,
                                             ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0, control.outer = list("trace" = FALSE))},
                            error = function(e) NULL)
          if(is.null(AObj)){
            AObj <- alabama::auglag(log(min(max(theta_het_H0, ResetThetaMin), ResetThetaMax)), probshift=0.5, fn = LogLikComp_het_H0, hin = eval_g_f_2_H0,
                                    ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0, control.outer = list("trace" = FALSE))
          }
          
          theta_het_H0 <- exp(AObj$par[1])
          
        }
      } else{
        OptObj <- tryCatch( {optim(par = log(min(max(theta_het_H0, ResetThetaMin), ResetThetaMax)), probshift=0.5, fn = LogLikComp_het_H0, gr = GradComp_het_H0, method = "BFGS",
                                   ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0)},
                            error = function(e) NULL)
        if(is.null(OptObj)){
          OptObj <- optim(par = log(min(max(theta_het_H0, ResetThetaMin), ResetThetaMax)), probshift=0.5, fn = LogLikComp_het_H0, method = "BFGS",
                          ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0)
        }
        
        theta_het_H0 <- exp(OptObj$par[1])
      }
      
    } else{
      
      if(sum(spr_H0+spv_H0)!=0){
        if(NoSplitHom){
          
          thetavec <- c()
          pivec <- c()
          likvec <- c()
          for(TH in thetaTRY){
            OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piRR_H0 >= 1-SE, 1-SE-eps, ifelse(piRR_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR_H0)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_RR, gr = GradComp_RR, 
                                       method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, pim = 1-SESlack*SE, M = 1-SE)},
                                error = function(e) NULL)
            if(is.null(OptObj)){
              OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piRR_H0 >= 1-SE, 1-SE-eps, ifelse(piRR_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR_H0)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_RR, 
                              method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, pim = 1-SESlack*SE, M = 1-SE)},
                              error = function(e) NULL)
              if(is.null(OptObj)){
                OptObj <- optim(par = c(gtools::logit(ifelse(piRR_H0 >= 1-SE, 1-SE-eps, ifelse(piRR_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR_H0)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_RR, 
                                method = "Nelder-Mead", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, pim = 1-SESlack*SE, M = 1-SE)
              }
            }
            
            pivec <- c(pivec, gtools::inv.logit(OptObj$par[1], max = 1-SE, min = 1-SESlack*SE))
            thetavec <- c(thetavec, exp(OptObj$par[2]))
            likvec <- c(likvec, LogLikComp_RR(parvec = c(OptObj$par[1], OptObj$par[2]), ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, pim = 1-SESlack*SE, M = 1-SE))
          }
          theta_RR_loc <- which(likvec == min(likvec))[1]
          piRR_H0 <- pivec[theta_RR_loc]
          theta_RR_H0 <- thetavec[theta_RR_loc]
          
          
          if(theta_RR_H0  > max(piRR_H0, 1-piRR_H0)){
            
            thetavec <- c()
            pivec <- c()
            likvec <- c()
            for(TH in thetaTRY){
              AObj <- tryCatch( {alabama::auglag(c(gtools::logit(ifelse(piRR_H0 >= 1-SE, 1-SE-eps, ifelse(piRR_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR_H0)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_RR, gr = GradComp_RR, 
                                                 hin = eval_g_f_RR, ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, pim = 1-SESlack*SE, M = 1-SE, 
                                                 control.outer = list("trace" = FALSE))}, error = function(e) NULL)
              if(is.null(AObj)){
                AObj <- alabama::auglag(c(gtools::logit(ifelse(piRR_H0 >= 1-SE, 1-SE-eps, ifelse(piRR_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR_H0)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_RR, 
                                        hin = eval_g_f_RR, ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, pim = 1-SESlack*SE, M = 1-SE, 
                                        control.outer = list("trace" = FALSE))
              }
              
              pivec <- c(pivec, gtools::inv.logit(AObj$par[1], max = 1-SE, min = 1-SESlack*SE))
              thetavec <- c(thetavec, exp(AObj$par[2]))
              likvec <- c(likvec, LogLikComp_RR(parvec = c(AObj$par[1], AObj$par[2]), ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, pim = 1-SESlack*SE, M = 1-SE))
            }
            theta_RR_loc <- which(likvec == min(likvec))[1]
            piRR_H0 <- pivec[theta_RR_loc]
            theta_RR_H0 <- thetavec[theta_RR_loc]
          }
          
          
          
          
          thetavec <- c()
          pivec <- c()
          likvec <- c()
          for(TH in thetaTRY){
            OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piVV_H0 >= 1-SE, 1-SE-eps, ifelse(piVV_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV_H0)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_VV, gr = GradComp_VV, 
                                       method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv_H0, pim = 1-SESlack*SE, M = 1-SE)},
                                error = function(e) NULL)
            if(is.null(OptObj)){
              OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piVV_H0 >= 1-SE, 1-SE-eps, ifelse(piVV_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV_H0)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_VV, 
                              method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv_H0, pim = 1-SESlack*SE, M = 1-SE)},
                              error = function(e) NULL)
              if(is.null(OptObj)){
                OptObj <- optim(par = c(gtools::logit(ifelse(piVV_H0 >= 1-SE, 1-SE-eps, ifelse(piVV_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV_H0)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_VV, 
                                method = "Nelder-Mead", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv_H0, pim = 1-SESlack*SE, M = 1-SE)
              }
            }
            
            pivec <- c(pivec, gtools::inv.logit(OptObj$par[1], max = 1-SE, min = 1-SESlack*SE))
            thetavec <- c(thetavec, exp(OptObj$par[2]))
            likvec <- c(likvec, LogLikComp_VV(parvec = c(OptObj$par[1], OptObj$par[2]), ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv_H0, pim = 1-SESlack*SE, M = 1-SE))
          }
          theta_VV_loc <- which(likvec == min(likvec))[1]
          piVV_H0 <- pivec[theta_VV_loc]
          theta_VV_H0 <- thetavec[theta_VV_loc]
          
          
          if(theta_VV_H0  > max(piVV_H0, 1-piVV_H0)){
            
            thetavec <- c()
            pivec <- c()
            likvec <- c()
            for(TH in thetaTRY){
              AObj <- tryCatch( {alabama::auglag(c(gtools::logit(ifelse(piVV_H0 >= 1-SE, 1-SE-eps, ifelse(piVV_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV_H0)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_VV, gr = GradComp_VV, 
                                                 hin = eval_g_f_VV, ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv_H0, pim = 1-SESlack*SE, M = 1-SE, 
                                                 control.outer = list("trace" = FALSE))}, error = function(e) NULL)
              if(is.null(AObj)){
                AObj <- alabama::auglag(c(gtools::logit(ifelse(piVV_H0 >= 1-SE, 1-SE-eps, ifelse(piVV_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV_H0)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_VV, 
                                        hin = eval_g_f_VV, ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv_H0, pim = 1-SESlack*SE, M = 1-SE, 
                                        control.outer = list("trace" = FALSE))
              }
              
              pivec <- c(pivec, gtools::inv.logit(AObj$par[1], max = 1-SE, min = 1-SESlack*SE))
              thetavec <- c(thetavec, exp(AObj$par[2]))
              likvec <- c(likvec, LogLikComp_VV(parvec = c(AObj$par[1], AObj$par[2]), ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv_H0, pim = 1-SESlack*SE, M = 1-SE))
            }
            theta_VV_loc <- which(likvec == min(likvec))[1]
            piVV_H0 <- pivec[theta_VV_loc]
            theta_VV_H0 <- thetavec[theta_VV_loc]
          }
          
        } else{
          
          thetavec <- c()
          pivec <- c()
          likvec <- c()
          for(TH in thetaTRY){
            OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piRR_H0 >= 1-SE, 1-SE-eps, ifelse(piRR_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR_H0)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_RR, gr = GradComp_RR, 
                                       method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, pim = 1-SESlack*SE, M = 1-SE)},
                                error = function(e) NULL)
            if(is.null(OptObj)){
              OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piRR_H0 >= 1-SE, 1-SE-eps, ifelse(piRR_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR_H0)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_RR, 
                              method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, pim = 1-SESlack*SE, M = 1-SE)},
                              error = function(e) NULL)
              if(is.null(OptObj)){
                OptObj <- optim(par = c(gtools::logit(ifelse(piRR_H0 >= 1-SE, 1-SE-eps, ifelse(piRR_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piRR_H0)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_RR, 
                                method = "Nelder-Mead", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, pim = 1-SESlack*SE, M = 1-SE)
              }
            }
            
            pivec <- c(pivec, gtools::inv.logit(OptObj$par[1], max = 1-SE, min = 1-SESlack*SE))
            thetavec <- c(thetavec, exp(OptObj$par[2]))
            likvec <- c(likvec, LogLikComp_RR(parvec = c(OptObj$par[1], OptObj$par[2]), ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spr = spr_H0, pim = 1-SESlack*SE, M = 1-SE))
          }
          theta_RR_loc <- which(likvec == min(likvec))[1]
          piRR <- pivec[theta_RR_loc]
          theta_RR <- thetavec[theta_RR_loc]
          
          thetavec <- c()
          pivec <- c()
          likvec <- c()
          for(TH in thetaTRY){
            OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piVV_H0 >= 1-SE, 1-SE-eps, ifelse(piVV_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV_H0)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_VV, gr = GradComp_VV, 
                                       method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv_H0, pim = 1-SESlack*SE, M = 1-SE)},
                                error = function(e) NULL)
            if(is.null(OptObj)){
              OptObj <- tryCatch( {optim(par = c(gtools::logit(ifelse(piVV_H0 >= 1-SE, 1-SE-eps, ifelse(piVV_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV_H0)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_VV, 
                              method = "BFGS", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv_H0, pim = 1-SESlack*SE, M = 1-SE)},
                              error = function(e) NULL)
              if(is.null(OptObj)){
                OptObj <- optim(par = c(gtools::logit(ifelse(piVV_H0 >= 1-SE, 1-SE-eps, ifelse(piVV_H0 <= 1-SESlack*SE, 1-SESlack*SE+eps, piVV_H0)), max = 1-SE, min = 1-SESlack*SE), log(TH)), fn = LogLikComp_VV, 
                                method = "Nelder-Mead", ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv_H0, pim = 1-SESlack*SE, M = 1-SE)
              }
            }
            
            pivec <- c(pivec, gtools::inv.logit(OptObj$par[1], max = 1-SE, min = 1-SESlack*SE))
            thetavec <- c(thetavec, exp(OptObj$par[2]))
            likvec <- c(likvec, LogLikComp_VV(parvec = c(OptObj$par[1], OptObj$par[2]), ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, spv = spv_H0, pim = 1-SESlack*SE, M = 1-SE))
          }
          theta_VV_loc <- which(likvec == min(likvec))[1]
          piVV_H0 <- pivec[theta_VV_loc]
          theta_VV_H0 <- thetavec[theta_VV_loc]
          
          
        }
        
      }
      
      #compareDerivatives(LogLikComp_hom, GradComp_hom, HessComp_hom, t0 = 0.3, SE=SE, ref_counts=ref_counts, var_counts=var_counts, spr=spr, spv=spv)
      #compareDerivatives(LogLikComp_het, GradComp_het, t0 = c(gtools::logit(1-0.002), log(0.5)), ref_counts=ref_counts, var_counts=var_counts, sprv = sprv)
      
      if(NoSplitHet){
        thetavec <- c()
        likvec <- c()
        for(TH in thetaTRY){
          OptObj <- tryCatch( {optim(par = log(TH), probshift=0.5, fn = LogLikComp_het_H0, gr = GradComp_het_H0, method = "BFGS",
                                     ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0)},
                              error = function(e) NULL)
          if(is.null(OptObj)){
            OptObj <- optim(par = log(TH), probshift=0.5, fn = LogLikComp_het_H0, method = "BFGS",
                            ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0)
          }
          
          thetavec <- c(thetavec, exp(OptObj$par[1]))
          likvec <- c(likvec, LogLikComp_het_H0(theta=OptObj$par[1], probshift=0.5, ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, sprv=sprv_H0))
        }
        
        theta_het_loc <- which(likvec == min(likvec))[1]
        theta_het_H0 <- thetavec[theta_het_loc]

        if(theta_het_H0  > max(0.5, 1-0.5)){
          thetavec <- c()
          likvec <- c()
          for(TH in thetaTRY){
            AObj <- tryCatch( {alabama::auglag(log(TH), probshift=0.5, fn = LogLikComp_het_H0, gr = GradComp_het_H0, hin = eval_g_f_2_H0,
                                               ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0, control.outer = list("trace" = FALSE))},
                              error = function(e) NULL)
            if(is.null(AObj)){
              AObj <- alabama::auglag(log(TH), probshift=0.5, fn = LogLikComp_het_H0, hin = eval_g_f_2_H0,
                                      ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0, control.outer = list("trace" = FALSE))
            }
            
            thetavec <- c(thetavec, exp(AObj$par[1]))
            likvec <- c(likvec, LogLikComp_het_H0(theta=AObj$par[1], probshift=0.5, ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, sprv=sprv_H0))
          }
          theta_het_loc <- which(likvec == min(likvec))[1]
          theta_het_H0 <- thetavec[theta_het_loc]
        }
        
      } else{
        thetavec <- c()
        likvec <- c()
        for(TH in thetaTRY){
          OptObj <- tryCatch( {optim(par = log(TH), probshift=0.5, fn = LogLikComp_het_H0, gr = GradComp_het_H0, method = "BFGS",
                                     ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0)},
                              error = function(e) NULL)
          if(is.null(OptObj)){
            OptObj <- optim(par = log(TH), probshift=0.5, fn = LogLikComp_het_H0, method = "BFGS",
                            ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv_H0)
          }
          
          thetavec <- c(thetavec, exp(OptObj$par[1]))
          likvec <- c(likvec, LogLikComp_het_H0(theta=OptObj$par[1], probshift=0.5, ref_counts=data_counts$ref_count, var_counts=data_counts$var_count, sprv=sprv_H0))
        }
        
        theta_het_loc <- which(likvec == min(likvec))[1]
        theta_het_H0 <- thetavec[theta_het_loc]
        
      }
      
    }
    
    Q <- sum(ifelse(pr_H0>0, spr_H0*log(pr_H0), 0) + spr_H0*dBetaBinom(data_counts$ref_count, data_counts$ref_count+data_counts$var_count, pi = piRR_H0, theta = theta_RR_H0, LOG = TRUE) +
               ifelse(prv_H0>0, sprv_H0*log(prv_H0), 0) + sprv_H0*dBetaBinom(data_counts$ref_count, data_counts$ref_count+data_counts$var_count, pi = 0.5, theta = theta_het_H0, LOG = TRUE) +
               ifelse(pv_H0>0, spv_H0*log(pv_H0), 0) + spv_H0*dBetaBinom(data_counts$var_count, data_counts$ref_count+data_counts$var_count, pi = piVV_H0, theta = theta_VV_H0, LOG = TRUE))
    
    dlta <- abs(Qold-Q)
    
    allelefreq_H0 <- mean(spr_H0) + mean(sprv_H0)/2
  }
  
  
  if(MyFlag == "!"){
    lrtstat <- NA
    pval <- NA
  } else{
    dmixh0 <- pmf_betabinomMix_SESlack(data_counts$ref_count, data_counts$var_count, probshift = 0.5, piRR_H0, piVV_H0,  mean(spr_H0), mean(spv_H0), mean(sprv_H0), 
                                         theta_RR = theta_RR_H0, theta_VV = theta_VV_H0, theta_het = theta_het_H0)
    dmixase <- pmf_betabinomMix_SESlack(data_counts$ref_count, data_counts$var_count, probshift = probshift, piRR, piVV, mean(spr), mean(spv), mean(sprv), 
                                          theta_RR = theta_RR, theta_VV = theta_VV, theta_het = theta_het)
    lrtstat <- -2 * (sum(log(dmixh0)) - sum(log(dmixase)))
    pval <- pchisq(lrtstat, df = 1, lower.tail = F)
  }
  
  data_counts$prr_H0 <- spr_H0
  data_counts$pvv_H0 <- spv_H0
  data_counts$prv_H0 <- sprv_H0
  
  results <- list(AB = probshift, AB_lrt = lrtstat, AB_p = pval, GOF = logLikelihood, GOFaltMEAN = GOFaltMEAN_S, GOFaltMEDIAN = GOFaltMEDIAN_S, 
                  GOFaltPERDIST = GOFaltPERDIST_S, GOFaltONLYHET = GOFaltONLYHET_S, GOFexactMEAN = GOFexactMEAN_S, GOFexactMEANLOG = GOFexactMEANLOG_S, 
                  nrep = nrep, quality = quality, rho_rr = pr, rho_vv = pv, rho_rv = prv, rho_rr_H0 = pr_H0, rho_vv_H0 = pv_H0, rho_rv_H0 = prv_H0, 
                  data_hash = data_counts, theta_RR = theta_RR, theta_VV = theta_VV, theta_het = theta_het, theta_RR_NoShift = theta_RR_H0,
                  theta_VV_NoShift = theta_VV_H0, theta_het_NoShift = theta_het_H0, piRR = piRR, piVV = piVV, piRR_NoShift = piRR_H0, piVV_NoShift = piVV_H0)
  return(results)
}




