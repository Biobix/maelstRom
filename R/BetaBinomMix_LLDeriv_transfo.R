#' Beta-binomial mixture model log-likelihood derivatives to its parameters involved in optimization
#'
#' @description \code{BetaBinomMix_LLDeriv} Lorem ipsum
#'
#' @param ref_counts Number or Numeric vector Reference count(s).
#' @param var_counts Number or Numeric vector. Variant count(s).
#' @param isCase Number or Numeric vector. Equals 1 for cases and 0 for controls.
#' @param probshift Number. The reference allele fraction in heterozygotes, indicating allelic bias when deviating from 0.5
#' @param SE Number. Sequencing error rate.
#' @param prr Number. Reference homozygote genotype probability of the locus.
#' @param pvv Number. Variant homozygote genotype probability of the locus.
#' @param prv Number. Heterozygote genotype probability of the locus.
#' @param theta_hom Number. The dispersion parameter of the homozygous peaks.
#' @param Der1 String. Derivate with respect to this parameter.
#' @param Der2 String. Derivate with respect to this parameter as well; optional.
#' @return A log-likelihood derivative
#' @export

BetaBinomMix_LLDeriv_transfo <- function(ref_counts, var_counts, probshift, SE, prr, pvv, prv, theta_hom, theta_het, Der1, Der2 = NULL){
  
  Lik_Complete <- maelstRom::pmf_betabinomMix(ref_counts, var_counts, probshift, SE, prr, pvv, prv, theta_hom, theta_het)
  
  PiLogit <- gtools::logit(probshift)
  #LogThetaHom <- log(theta_hom)
  #LogThetaHet <- log(theta_het)
  
  PiCorr <- (exp(PiLogit)/(1+exp(PiLogit))^2)
  #ThetaCorrs are just exp(log()), so the thetas themselves
  
  PiCorr2 <- -exp(PiLogit)*(-1+exp(PiLogit))/(1+exp(PiLogit))^3
  
  if(is.null(Der2)){
    
    if(Der1 == "probshift"){
      
      PreVec <- ((prv*maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het, LOG = FALSE))/Lik_Complete)
      
      if(any(is.infinite(PreVec))){
        rrterm <- (log(prr)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, 1-SE, theta_hom, LOG = TRUE))
        rvterm <- (log(prv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het, LOG = TRUE))
        vvterm <- (log(pvv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom, LOG = TRUE))
        PreVec[is.infinite(PreVec) & rvterm > rrterm & rvterm > vvterm] <- 1
        PreVec[is.infinite(PreVec)] <- 0
      }
      
      OUT <- PreVec*(maelstRom::grad_pi(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het) * PiCorr)
      
    } else if (Der1 == "theta_hom"){
      
      PreVec1 <- ((prr*maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, 1-SE, theta_hom, LOG = FALSE))/Lik_Complete)
      PreVec2 <- ((pvv*maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom, LOG = FALSE))/Lik_Complete)
      
      if(any(is.infinite(PreVec1))){
        rrterm <- (log(prr)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, 1-SE, theta_hom, LOG = TRUE))
        rvterm <- (log(prv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het, LOG = TRUE))
        vvterm <- (log(pvv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom, LOG = TRUE))
        PreVec[is.infinite(PreVec1) & rrterm > rvterm & rrterm > vvterm] <- 1
        PreVec[is.infinite(PreVec1)] <- 0
      }
      if(any(is.infinite(PreVec2))){
        rrterm <- (log(prr)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, 1-SE, theta_hom, LOG = TRUE))
        rvterm <- (log(prv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het, LOG = TRUE))
        vvterm <- (log(pvv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom, LOG = TRUE))
        PreVec[is.infinite(PreVec2) & vvterm > rvterm & vvterm > rrterm] <- 1
        PreVec[is.infinite(PreVec2)] <- 0
      }
      
      OUT <- PreVec1*(maelstRom::grad_theta(ms = ref_counts, ns = ref_counts+var_counts, 1-SE, theta_hom) * theta_hom) +
        PreVec2*(maelstRom::grad_theta(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom) * theta_hom)
      
    } else if (Der1 == "theta_het"){
      
      PreVec <- ((prv*maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het, LOG = FALSE))/Lik_Complete)

      if(any(is.infinite(PreVec))){
        rrterm <- (log(prr)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, 1-SE, theta_hom, LOG = TRUE))
        rvterm <- (log(prv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het, LOG = TRUE))
        vvterm <- (log(pvv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom, LOG = TRUE))
        PreVec[is.infinite(PreVec) & rvterm > rrterm & rvterm > vvterm] <- 1
        PreVec[is.infinite(PreVec)] <- 0
      }
      
      OUT <- PreVec*(maelstRom::grad_theta(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het) * theta_het)
      
    } else if (Der1 == "prr"){
      
      PreVec <- ((prr*maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, 1-SE, theta_hom, LOG = FALSE))/Lik_Complete)
      
      if(any(is.infinite(PreVec))){
        rrterm <- (log(prr)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, 1-SE, theta_hom, LOG = TRUE))
        rvterm <- (log(prv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het, LOG = TRUE))
        vvterm <- (log(pvv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom, LOG = TRUE))
        PreVec[is.infinite(PreVec) & rrterm > rvterm & rrterm > vvterm] <- 1
        PreVec[is.infinite(PreVec)] <- 0
      }
      
      PreVec2 <- ((prr*maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom, LOG = FALSE))/Lik_Complete)
      
      if(any(is.infinite(PreVec2))){
        rrterm <- (log(prr)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, 1-SE, theta_hom, LOG = TRUE))
        rvterm <- (log(prv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het, LOG = TRUE))
        vvterm <- (log(pvv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom, LOG = TRUE))
        PreVec2[is.infinite(PreVec2) & vvterm > rvterm & vvterm > rrterm] <- prr/pvv
        PreVec2[is.infinite(PreVec2)] <- 0
      }
      
      OUT <- PreVec*(1/prr) - PreVec2*(1/prr)
      
    } else if (Der1 == "prv"){
      
      PreVec <- ((prv*maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het, LOG = FALSE))/Lik_Complete)
      
      if(any(is.infinite(PreVec))){
        rrterm <- (log(prr)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, 1-SE, theta_hom, LOG = TRUE))
        rvterm <- (log(prv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het, LOG = TRUE))
        vvterm <- (log(pvv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom, LOG = TRUE))
        PreVec[is.infinite(PreVec) & rvterm > rrterm & rvterm > vvterm] <- 1
        PreVec[is.infinite(PreVec)] <- 0
      }
      
      PreVec2 <- ((prv*maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom, LOG = FALSE))/Lik_Complete)
      
      if(any(is.infinite(PreVec2))){
        rrterm <- (log(prr)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, 1-SE, theta_hom, LOG = TRUE))
        rvterm <- (log(prv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het, LOG = TRUE))
        vvterm <- (log(pvv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom, LOG = TRUE))
        PreVec2[is.infinite(PreVec2) & vvterm > rvterm & vvterm > rrterm] <- prv/pvv
        PreVec2[is.infinite(PreVec2)] <- 0
      }
      
      OUT <- PreVec*(1/prv) - PreVec2*(1/prv)
      
    }
    
  }else{
    
    Aterm <- - BetaBinomMix_LLDeriv_transfo(ref_counts, var_counts, probshift, SE, prr, pvv, prv, theta_hom, theta_het, Der1 = Der1, Der2 = NULL) *
      BetaBinomMix_LLDeriv_transfo(ref_counts, var_counts, probshift, SE, prr, pvv, prv, theta_hom, theta_het, Der1 = Der2, Der2 = NULL)
    
    #if(Der1 == "theta_hom"){
    #  CorrCorr <- theta_hom
    #  if(Der2 == "theta_hom"){
    #    CorrCorrD <- theta_hom
    #  }else{
    #    CorrCorrD <- 0
    #  }
    #}else if(Der1 == "theta_het"){
    #  CorrCorr <- theta_het
    #  if(Der2 == "theta_het"){
    #    CorrCorrD <- theta_het
    #  }else{
    #    CorrCorrD <- 0
    #  }
    #}else if(Der1 == "probshift"){
    #  CorrCorr <- PiCorr
    #  if(Der2 == "probshift"){
    #    CorrCorrD <- PiCorr2
    #  }else{
    #    CorrCorrD <- 0
    #  }
    #}else{
    #  CorrCorr <- 1
    #  CorrCorrD <- 0
    #}
    #
    #Aterm <- Aterm*CorrCorr + CorrCorrD/Lik_Complete
    
    PreVec_rr <- ((prr*maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, 1-SE, theta_hom, LOG = FALSE))/Lik_Complete)
    PreVec_rv <- ((prv*maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het, LOG = FALSE))/Lik_Complete)
    PreVec_vv <- ((pvv*maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het, LOG = FALSE))/Lik_Complete)

    if(any(is.infinite(PreVec_rr))){
      rrterm <- (log(prr)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, 1-SE, theta_hom, LOG = TRUE))
      rvterm <- (log(prv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het, LOG = TRUE))
      vvterm <- (log(pvv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom, LOG = TRUE))
      PreVec_rr[is.infinite(PreVec_rr) & rrterm > rvterm & rrterm > vvterm] <- 1
      PreVec_rr[is.infinite(PreVec_rr)] <- 0
    }
    if(any(is.infinite(PreVec_rv))){
      rrterm <- (log(prr)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, 1-SE, theta_hom, LOG = TRUE))
      rvterm <- (log(prv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het, LOG = TRUE))
      vvterm <- (log(pvv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom, LOG = TRUE))
      PreVec_rv[is.infinite(PreVec_rv) & rvterm > rrterm & rvterm > vvterm] <- 1
      PreVec_rv[is.infinite(PreVec_rv)] <- 0
    }
    if(any(is.infinite(PreVec_vv))){
      rrterm <- (log(prr)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, 1-SE, theta_hom, LOG = TRUE))
      rvterm <- (log(prv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het, LOG = TRUE))
      vvterm <- (log(pvv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom, LOG = TRUE))
      PreVec_vv[is.infinite(PreVec_vv) & vvterm > rrterm & vvterm > rvterm] <- 1
      PreVec_vv[is.infinite(PreVec_vv)] <- 0
    }
    
    
    PreVec2_rr <- ((prr*maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom, LOG = FALSE))/Lik_Complete)
    PreVec2_rv <- ((prv*maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom, LOG = FALSE))/Lik_Complete)
    
    if(any(is.infinite(PreVec2_rr))){
      rrterm <- (log(prr)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, 1-SE, theta_hom, LOG = TRUE))
      rvterm <- (log(prv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het, LOG = TRUE))
      vvterm <- (log(pvv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom, LOG = TRUE))
      PreVec2_rr[is.infinite(PreVec2_rr) & vvterm > rvterm & vvterm > rrterm] <- prr/pvv
      PreVec2_rr[is.infinite(PreVec2_rr)] <- 0
    }
    if(any(is.infinite(PreVec2_rv))){
      rrterm <- (log(prr)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, 1-SE, theta_hom, LOG = TRUE))
      rvterm <- (log(prv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het, LOG = TRUE))
      vvterm <- (log(pvv)+maelstRom::dBetaBinom(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom, LOG = TRUE))
      PreVec2_rv[is.infinite(PreVec2_rv) & vvterm > rvterm & vvterm > rrterm] <- prv/pvv
      PreVec2_rv[is.infinite(PreVec2_rv)] <- 0
    }
    
    
    if(Der1 == "probshift"){
      
      if(Der2 == "probshift"){
        
        Bterm <- PreVec_rv*((maelstRom::grad_pi(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het)*PiCorr)^2 +
                             (grad_pi_pi(ms = ref_counts, ns = ref_counts+var_counts, pi = probshift, theta = theta_het)*PiCorr^2 +
                              maelstRom::grad_pi(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het)*PiCorr2))
        
      } else if (Der2 == "theta_hom"){
        
        Bterm <- 0
        
      } else if (Der2 == "theta_het"){
        
        Bterm <- PreVec_rv*((maelstRom::grad_theta(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het)*theta_het *
                             maelstRom::grad_pi(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het)*PiCorr) +
                             grad_pi_theta(ms = ref_counts, ns = ref_counts+var_counts, pi = probshift, theta = theta_het)*PiCorr*theta_het)
        
      } else if (Der2 == "prr"){
        
        Bterm <- 0
        
      } else if (Der2 == "prv"){
        
        Bterm <- PreVec_rv*(((1/prv) *
                             maelstRom::grad_pi(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het)*PiCorr))
        
      }
      
    } else if (Der1 == "theta_hom"){
      
      if(Der2 == "probshift"){
        
        Bterm <- 0
        
      } else if (Der2 == "theta_hom"){
        
        Bterm <- PreVec_rr * ((maelstRom::grad_theta(ms = ref_counts, ns = ref_counts+var_counts, 1-SE, theta_hom)*theta_hom)^2 +
                               (grad_theta_theta(ms = ref_counts, ns = ref_counts+var_counts, pi = 1-SE, theta = theta_hom)*theta_hom^2 +
                                maelstRom::grad_theta(ms = ref_counts, ns = ref_counts+var_counts, 1-SE, theta_hom)*theta_hom)) +
                 PreVec_vv * ((maelstRom::grad_theta(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom)*theta_hom)^2 +
                               (grad_theta_theta(ms = ref_counts, ns = ref_counts+var_counts, pi = SE, theta = theta_hom)*theta_hom^2 +
                                maelstRom::grad_theta(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom)*theta_hom))
        
      } else if (Der2 == "theta_het"){
        
        Bterm <- 0
        
      } else if (Der2 == "prr"){
        
        Bterm <- PreVec_rr * (((1/prr)*maelstRom::grad_theta(ms = ref_counts, ns = ref_counts+var_counts, 1-SE, theta_hom)*theta_hom)) -
          PreVec2_rr * (((1/prr)*maelstRom::grad_theta(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom)*theta_hom))
        
      } else if (Der2 == "prv"){
        
        Bterm <- - PreVec2_rv * (((1/prv)*maelstRom::grad_theta(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom)*theta_hom))
        
      }
      
    } else if (Der1 == "theta_het"){
      
      if(Der2 == "probshift"){
        
        Bterm <- PreVec_rv*((maelstRom::grad_theta(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het)*theta_het *
                               maelstRom::grad_pi(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het)*PiCorr) +
                              grad_pi_theta(ms = ref_counts, ns = ref_counts+var_counts, pi = probshift, theta = theta_het)*PiCorr*theta_het)
        
      } else if (Der2 == "theta_hom"){
        
        Bterm <- 0
        
      } else if (Der2 == "theta_het"){
        
        Bterm <- PreVec_rv * ((maelstRom::grad_theta(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het)*theta_het)^2 +
                               (grad_theta_theta(ms = ref_counts, ns = ref_counts+var_counts, pi = probshift, theta = theta_het)*theta_het^2 +
                                  maelstRom::grad_theta(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het)*theta_het))
        
      } else if (Der2 == "prr"){
        
        Bterm <- 0
        
      } else if (Der2 == "prv"){
        
        Bterm <- PreVec_rv * ((1/prv)*(maelstRom::grad_theta(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het)*theta_het))
        
      }
      
    } else if (Der1 == "prr"){
      
      if(Der2 == "probshift"){
        
        Bterm <- 0
        
      } else if (Der2 == "theta_hom"){
        
        Bterm <- PreVec_rr * (((1/prr)*maelstRom::grad_theta(ms = ref_counts, ns = ref_counts+var_counts, 1-SE, theta_hom)*theta_hom)) -
          PreVec2_rr * (((1/prr)*maelstRom::grad_theta(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom)*theta_hom))
        
      } else if (Der2 == "theta_het"){
        
        Bterm <- 0
        
      } else if (Der2 == "prr"){
        
        Bterm <- 0
        
      } else if (Der2 == "prv"){
        
        Bterm <- 0
        
      }
      
    } else if (Der1 == "prv"){
      
      if(Der2 == "probshift"){
        
        Bterm <- PreVec_rv*(((1/prv) *
                               maelstRom::grad_pi(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het)*PiCorr))
        
        
      } else if (Der2 == "theta_hom"){
        
        Bterm <- - PreVec2_rv * (((1/prv)*maelstRom::grad_theta(ms = ref_counts, ns = ref_counts+var_counts, SE, theta_hom)*theta_hom))
        
      } else if (Der2 == "theta_het"){
        
        Bterm <- PreVec_rv * ((1/prv)*(maelstRom::grad_theta(ms = ref_counts, ns = ref_counts+var_counts, probshift, theta_het)*theta_het))
        
      } else if (Der2 == "prr"){
        
        Bterm <- 0
        
      } else if (Der2 == "prv"){
        
        Bterm <- 0
        
      }
      
    }
    OUT <- Aterm + Bterm
    
  }
  
  return(OUT)
  
}
