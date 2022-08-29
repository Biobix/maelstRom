#' Calculates moment estimates for the beta-binomial parameters given observed data
#'
#' \code{MomentEst_MixedBetaBinom} calculates moment estimates for the pi- and theta-parameters of the beta-binomial data,
#' given data with a varying number of observations. This estimator was previously described by \insertCite{kleinman}{MAGE},
#' but was adapted into a weighted (according to how likely a sample is to be homozygous reference, homozygous variant, or
#' heterozygous) estimator to be used in an expectation-maximization-framework. Due to the varying number of observations, this moment estimate
#' is actually rather tricky should incorporate an additional per-sample weight that is to be tuned via an iterative procedure,
#' but given this function's aim is just obtaining a rough estimate as starting point for numerical expectation maximization,
#' these weights are all set at one (which corresponds to the ideal weights in case the true theta equals zero).
#' 
#' @importFrom Rdpack reprompt
#' @param ref_counts Numeric vector. Reference counts.
#' @param var_counts Numeric vector. Variant counts.
#' @param spr Numeric vector. Chances for samples to be reference homozygotes
#' @param spv Numeric vector. Chances for samples to be variant homozygotes
#' @param sprv Numeric vector. Chances for samples to be heterozygotes
#' @param pi_hom_fix Number. If set to a value, fixes the pi-parameter of the homozygous peaks instead of generating a moment estimate;
#'    this is an option because the moment estimate of the theta-parameter depends on the pi-parameter.
#' @param pi_het_fix Number. If set to a value, fixes the pi-parameter of the heterozygous peak instead of generating a moment estimate;
#'    this is an option because the moment estimate of the theta-parameter depends on the pi-parameter.
#' @export
#' @return A named vector containing moment estimates of the pi- and theta-parameters for both the homozygous and heterozygous peaks.
#' @references
#'    \insertAllCited{}


MomentEst_MixedBetaBinom <- function(ref_counts, var_counts, spr, spv, sprv, pi_hom_fix = NULL, pi_het_fix = NULL){
  
  #homWeights <- rep(1, length(ref_counts)*2)
  #hetWeigths <- rep(1, length(ref_counts))
  
  tot_counts <- ref_counts + var_counts
  
  Pref <- ref_counts/tot_counts
  Pvar <- var_counts/tot_counts
  
  if(is.null(pi_hom_fix)){
    pi_hom_ME <- sum( c(spr, spv) * c(Pref, Pvar) ) / sum( c(spr, spv) )
  } else{
    pi_hom_ME <- pi_hom_fix
  }
  
  if(is.null(pi_het_fix)){
    pi_het_ME <- sum(  sprv * Pref ) / sum( sprv )
  } else{
    pi_het_ME <- pi_het_fix
  }
  
  SS_hom <- sum( c(spr, spv) * (pi_hom_ME - c(Pref, Pvar))^2) / sum(c(spr, spv))
  SS_het <- sum( sprv * (pi_het_ME - Pref)^2) / sum(sprv)
  
  #z_hom <- 1 - ( c(spr, spv) / sum(c(spr, spv)) )
  #z_het <- 1 - (sprv / sum(sprv))
  z_hom <- 1 - ( rep(1, length(c(spr,spv))) / sum(rep(1, length(c(spr,spv)))) )
  z_het <- 1 - (rep(1, length(sprv)) / sum(rep(1, length(sprv))))
  
  theta_hom_ME <- 1 / ( (pi_hom_ME * (1-pi_hom_ME) * ((sum(c(spr, spv) * z_hom) - sum(c(spr, spv) * z_hom / rep(tot_counts,2)))/sum(c(spr, spv))) /
                           SS_hom - pi_hom_ME * (1-pi_hom_ME) * (sum(c(spr, spv) * z_hom / rep(tot_counts,2))/sum(c(spr, spv))) ) - 1)
  theta_het_ME <- 1 / ( (pi_het_ME * (1-pi_het_ME) * ((sum(sprv * z_het) - sum(sprv * z_het / tot_counts))/sum(sprv)) /
                           SS_het - pi_het_ME * (1-pi_het_ME) * (sum(sprv * z_het / tot_counts)/sum(sprv)) ) - 1)
  
  OutVec <- c(pi_hom_ME, theta_hom_ME, pi_het_ME, theta_het_ME)
  names(OutVec) <- c("pi_hom", "theta_hom", "pi_het", "theta_het")
  
  return(OutVec)
  
}
