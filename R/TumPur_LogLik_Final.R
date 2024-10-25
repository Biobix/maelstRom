TumPur_LogLik_Final <- function(optpars, ref_counts, var_counts, tumpur, weights, SCPthreshold = 0.999){
  
  a <- exp(optpars[1])
  b <- exp(optpars[2])
  q <- exp(optpars[3])
  
  SampLiks <- rep(0, length(ref_counts))
  
  for(k in 1:length(SampLiks)){
    #print(k)
    
    RC <- ref_counts[k]
    VC <- var_counts[k]
    TC <- RC+VC
    
    TP <- tumpur[k]
    
    SampCompProbs <- dbinom(x = 0:TC, size = TC, prob = TP)
    
    SCP_sort <- sort(SampCompProbs, decreasing = TRUE)
    SCP_order <- order(SampCompProbs, decreasing = TRUE)
    
    SCP_cutoff <- (which(cumsum(SCP_sort) >= SCPthreshold))[1]
    
    ToKeep <- SCP_order[1:SCP_cutoff]
    SCP_oi <- SampCompProbs[ToKeep]
    TumReads_oi <- sort((0:TC)[ToKeep])
    
    #SampLiks[k] <- TumPurHelpFun_Final(RC, TC, TumReads_oi, a, b, q, SCP_oi) # The likelihood result for one sample...
    RCcol <- 0:RC
    SampLiks[k] <- TumPurHelpFun_CPP(RCcol, TC, TumReads_oi, a, b, q, SCP_oi, n)
  }
  
  return(-(sum(weights*SampLiks)))
  
}
