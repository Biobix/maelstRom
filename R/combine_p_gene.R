#' Combine p-values of SNPs per gene.
#'
#' \code{combine_p_gene} combines per-SNP p-values into per-gene p-values;
#' the method for doing so can be chosen via input
#'
#' @param p_values Numeric vector. P-values to combine.
#' @param weights Numeric vector. Weights given to each observation (corresponds to how much an observation occurs, but can be non-integer as well).
#' @param method String. Specifies the method for combining p-values, the options being "arithmetic" (mean), "geometric" (mean),
#' "harmonic" (mean), "minimum" (p-value), "maximum" (p-value) or "GMP" (Generalized Mean P-value).
#' In case of the latter, a value for r must be specified, and the GMP is then calculated as (sum(weights * p_values^r)/sum(weights))^(1/r), see \insertCite{wilson2020generalized}{MAGE}.
#' @param r Numeric. r-value to be used in the calculation of the GMP, should that option be chosen as the method.
#' @export
#' @return Combined p-value of a gene.
#' @examples
#' combine_p_gene(c(0.015, 0.414, 8.47E-03, 7.43E-03, 0.574, 0.837))
#' combine_p_gene(c(0.063, 0.725, 0.657, 0.378, 0.291))
#' @references
#'    \insertAllCited{}

combine_p_gene <- function (p_values, weights = NULL, method = "geometric", r = NULL) {
  
  if(is.null(weights)){
    weights <- rep(1, length(p_values))
  }
  
  if(method == "geometric"){
    combined_p <- exp(sum(weights * log(p_values))/sum(weights))
  } else if(method == "harmonic"){
    combined_p <- sum(weights) / sum(weights / p_values)
  } else if(method == "arithmetic"){
    combined_p <- sum(weights*p_values)/sum(weights)
  } else if(method == "minimum"){
    combined_p <- min(p_values)
  } else if(method == "maximum"){
    combined_p <- max(p_values)
  } else if(method == "GMP"){
    if(is.null(r) | !is.numeric(r)){
      stop("Please assign an input value to r to use in the GMP calculation")
    }
    if(r == 0){
      combined_p <- exp(sum(weights * log(p_values))/sum(weights))
    } else{
      combined_p <- (sum(weigths * p_values^r)/sum(weigths))^(1/r)
    }
  }
  return(combined_p)
}
