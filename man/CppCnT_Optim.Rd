% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{CppCnT_Optim}
\alias{CppCnT_Optim}
\title{Optimizer for pi, theta_control, and theta_case, for heterozygous samples in the differential control-case beta-binomial mixture fit using expectation-maximization}
\usage{
CppCnT_Optim(
  StartVals,
  ref_counts,
  var_counts,
  isCase,
  sprv,
  MemLim = 2048L,
  Xtra = 7L,
  step_size = 0.01,
  tol = 0.1,
  epsabs = 0.001
)
}
\description{
Optimizer for pi, theta_control, and theta_case, for heterozygous samples in the differential control-case beta-binomial mixture fit using expectation-maximization
}
