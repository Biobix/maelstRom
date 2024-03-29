% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LikelyDistsHet.R
\name{LikelyDistsHet}
\alias{LikelyDistsHet}
\title{Returns parameter- and likelihood-distances based on leave-on-out (case-deletion) MLE's of \code{EMfit_betabinom_robust}'s fitted model}
\usage{
LikelyDistsHet(
  ref_counts,
  var_counts,
  sprv,
  parvec_cur,
  NoSplitHet,
  ResetThetaMin,
  ResetThetaMax,
  SE,
  ReEstPars = FALSE
)
}
\arguments{
\item{ref_counts}{Numeric vector. reference counts.}

\item{var_counts}{Numeric vector. variant counts.}

\item{sprv}{Numeric vector. Each sample's EM-weight reflecting its likelihood to be part of the heterozygous population.}

\item{parvec_cur}{Numeric vector. Pi and theta (in that order) of the heterozyous peak of the full-data fit.}

\item{NoSplitHet}{Logical. If TRUE, don't allow the beta-binomial fit for heterozygotes to be bimodal}

\item{ResetThetaMin}{Number. Initial theta values in numeric optimization get capped at this minimum (e.g. in case the moment estimate is even lower)}

\item{ResetThetaMax}{Number. Initial theta values in numeric optimization get capped at this maximum (e.g. in case the moment estimate is even higher)}

\item{SE}{Number. Sequencing error rate.}

\item{ReEstPars}{Logical. If TRUE, re-estimates \code{parvec_cur} given \code{ref_counts} and \code{var_counts}. This is useless if these are the actual
counts of the full dataset, but are useful for an emperical approach in which "expected" parameter- and likelihood-distances if the assumed model is 100%
correct are simulated by drawing \code{ref_counts} and \code{var_counts} from this assumed model (see \code{EMfit_betabinom_robust})}
}
\value{
A list containing the following components:
\item{LikDists}{A vector containing likelihood distances per sample (2 times full-data log-likelihood minus re-fitted log-likelihood leaving out the sample).}
\item{PiDists}{A vector containing pi distances per sample (leave-sample-out refitted pi minus full-data pi).}
\item{ThetaDists}{A vector containing theta distances per sample (leave-sample-out refitted theta minus full-data theta).}
}
\description{
\code{LikelyDistsHet} is mainly for internal use by \code{EMfit_betabinom_robust}. It assists in robustifying the function's EM-fit by iteratively re-fitting
the model on the entire input dataset except for one point, after which the difference in heterozygous pi- and theta-estimates and likelihoods is logged.
The difference with their respective full-data fit counterparts is a measure for the left-out data point's influence on the model fit;
if either one is sufficiently high the data point could be considered an outlier.
}
