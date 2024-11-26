#' Similar to MomentEst_MixedBetaBinom, but allowing for slack on the SE parameter (see EMfit_betabinom_SEslack)

MomentEst_MixedBetaBinom_SESlack <- function(ref_counts, var_counts, spr, spv, sprv, pi_RR_fix = NULL, pi_VV_fix = NULL, pi_het_fix = NULL){
  
  #homWeights <- rep(1, length(ref_counts)*2)
  #hetWeigths <- rep(1, length(ref_counts))
  
  tot_counts <- ref_counts + var_counts
  
  Pref <- ref_counts/tot_counts
  Pvar <- var_counts/tot_counts
  
  if(is.null(pi_RR_fix)){
    pi_RR_ME <- sum( spr * Pref ) / sum( spr )
  } else{
    pi_RR_ME <- pi_RR_fix
  }
  
  if(is.null(pi_VV_fix)){
    pi_VV_ME <- sum( spv * Pvar ) / sum( spv )
  } else{
    pi_VV_ME <- pi_VV_fix
  }
  
  #if(is.null(pi_hom_fix)){
  #  pi_hom_ME <- sum( c(spr, spv) * c(Pref, Pvar) ) / sum( c(spr, spv) )
  #} else{
  #  pi_hom_ME <- pi_hom_fix
  #}
  
  if(is.null(pi_het_fix)){
    pi_het_ME <- sum(  sprv * Pref ) / sum( sprv )
  } else{
    pi_het_ME <- pi_het_fix
  }
  
  SS_RR <- sum( spr * (pi_RR_ME - Pref)^2) / sum(spr)
  SS_VV <- sum( spv * (pi_VV_ME - Pvar)^2) / sum(spv)
  #SS_hom <- sum( c(spr, spv) * (pi_hom_ME - c(Pref, Pvar))^2) / sum(c(spr, spv))
  SS_het <- sum( sprv * (pi_het_ME - Pref)^2) / sum(sprv)
  
  #z_hom <- 1 - ( c(spr, spv) / sum(c(spr, spv)) )
  #z_het <- 1 - (sprv / sum(sprv))
  
  z_RR <- 1 - (rep(1, length(spr)) / sum(rep(1, length(spr))))
  z_VV <- 1 - (rep(1, length(spv)) / sum(rep(1, length(spv))))
  #z_hom <- 1 - ( rep(1, length(c(spr,spv))) / sum(rep(1, length(c(spr,spv)))) )
  z_het <- 1 - (rep(1, length(sprv)) / sum(rep(1, length(sprv))))
  
  theta_RR_ME <- 1 / ( (pi_RR_ME * (1-pi_RR_ME) * ((sum(spr * z_RR) - sum(spr * z_RR / tot_counts))/sum(spr)) /
                           SS_RR - pi_RR_ME * (1-pi_RR_ME) * (sum(spr * z_RR / tot_counts)/sum(spr)) ) - 1)
  theta_VV_ME <- 1 / ( (pi_VV_ME * (1-pi_VV_ME) * ((sum(spv * z_VV) - sum(spv * z_VV / tot_counts))/sum(spv)) /
                          SS_VV - pi_VV_ME * (1-pi_VV_ME) * (sum(spv * z_VV / tot_counts)/sum(spv)) ) - 1)
  #theta_hom_ME <- 1 / ( (pi_hom_ME * (1-pi_hom_ME) * ((sum(c(spr, spv) * z_hom) - sum(c(spr, spv) * z_hom / rep(tot_counts,2)))/sum(c(spr, spv))) /
  #                         SS_hom - pi_hom_ME * (1-pi_hom_ME) * (sum(c(spr, spv) * z_hom / rep(tot_counts,2))/sum(c(spr, spv))) ) - 1)
  theta_het_ME <- 1 / ( (pi_het_ME * (1-pi_het_ME) * ((sum(sprv * z_het) - sum(sprv * z_het / tot_counts))/sum(sprv)) /
                           SS_het - pi_het_ME * (1-pi_het_ME) * (sum(sprv * z_het / tot_counts)/sum(sprv)) ) - 1)
  
  OutVec <- c(pi_RR_ME, pi_VV_ME, theta_RR_ME, theta_VV_ME, pi_het_ME, theta_het_ME)
  names(OutVec) <- c("piRR", "piVV", "theta_RR", "theta_VV", "pi_het", "theta_het")
  
  return(OutVec)
  
}
