#' @name TumPur_LogLik_CPPs
#' @title Various implementations of the tumor-purity-accounting (negative) log likelihood (sum of completely correlated beta-binomials)
#'
#' @description UNUSED. These different implementations of, in essence, the same log-likelihood (to be used in maximum likelihood estimation) use
#' different parameters transformations in order to find out which one is optimal for numerical optimization.
#'
#' @param optpars Numeric vector. To be optimized parameters, usually consisting of the beta-binomial alpha (a) and beta (b) parameter,
#' and a q-parameters denoting their difference (multiplicative factor) in cases compared to controls (same pi assumed, different theta;
#' complete correlation within a sample assumed as well).
#' @param ref_counts Numeric vector. Vector of reference allele counts.
#' @param var_counts Numeric vector. Vector of variant allele counts.
#' @param tumpur Numeric vector. Vector of tumor purities (or an estimate thereof) at expression-level.
#' 
#' @param weights Numeric vector. Weights to be used for returning a weighted log-likelihood (sample-wise multiplication of log-likelihood with these weights).
#' @param SCPthreshold Number. Every sample's corresponding tumor purity is used to model binomial sampling of tumor- and control reads making up the tumor observations.
#' More specifically, every possible sampling's probability is tracked, and its corresponding log-likelihood (of the sum of completely correlated beta-binomials) calculated.
#' This might be a bit overkill however; when SCPthreshold is specified, only the most likely binomial read sampling are taken into account, up until their total binomial
#' probability sums to SCPthreshold. When SCPthreshold is set to 0, only the single most likely sampling is taken into account.
#' @param n Number. Number of intervals to be used during calculations of sum-of-completely-correlated-beta-binomial likelihood via numerical integration,
#' when the used integration method (NumIntMethod) is either the modified trapezoidal method using Gregory weights ("Gregory") or Tanh-sinh quadrature ("TanhSinhQuad")
#' @param NumIntMethod String. Used numeric integration method (see above); either the modified trapezoidal method using Gregory weights ("Gregory"),
#' Tanh-sinh quadrature ("TanhSinhQuad"), Newton Cotes ("NewtonCotes"), or Gaussian quadrature ("GaussianQuad").
#' @param prec Number. A low number (e.g. 0.0001) acting as stop criterion when NumIntMethod == "TanhSinhQuad"
#' (its iterative calculations stops if the current calculation differs less then prec from the previous one).
#' @param Wvec Numeric vector. Vector of numeric integration weights used when NumIntMethod == "NewtonCotes", or NumIntMethod == "GaussianQuad".
#' @param Nvec Numeric vector. Vector of node values used when NumIntMethod == "GaussianQuad".
#' @return Negative log-likelihood.

NULL

