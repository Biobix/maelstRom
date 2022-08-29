#' Calculate multinomial coefficient(s).
#'
#'\code{multinomial_coeff} calculates multinomial coefficients given locus data, mainly for internal use by \code{imprinting_est}.
#'
#' @param ref_counts Numeric vector. Reference count(s).
#' @param var_counts Numeric vector. Variant count(s).
#' @return The multinomial coefficient(s) of \code{ref_counts} and \code{var_counts}.
#' @export
#' @examples
#' multinomial_coeff(5, 11)
#' multinomial_coeff(c(5, 8, 10, 3, 5, 6, 23), c(8, 8, 6, 4, 4, 10, 0))

multinomial_coeff <- function(ref_counts, var_counts) {
  x_m <- gmp::factorialZ(ref_counts + var_counts) / (gmp::factorialZ(ref_counts) * gmp::factorialZ(var_counts))
  return(x_m)
}
