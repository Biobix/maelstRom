% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmf_betabinomMix.R
\name{pmf_betabinomMix}
\alias{pmf_betabinomMix}
\title{Probability Mass Function of the beta-binomial mixture distribution modeling population-level RNAseq data}
\usage{
pmf_betabinomMix(
  ref_counts,
  var_counts,
  probshift,
  SE,
  pr,
  pv,
  prv,
  theta_hom,
  theta_het
)
}
\arguments{
\item{ref_counts}{Number or Numeric vector Reference count(s).}

\item{var_counts}{Number or Numeric vector. Variant count(s).}

\item{probshift}{Number. The reference allele fraction in heterozygotes, indicating allelic bias when deviating from 0.5}

\item{SE}{Number. Sequencing error rate.}

\item{pr}{Number. Reference homozygote genotype probability of the locus.}

\item{pv}{Number. Variant homozygote genotype probability of the locus.}

\item{prv}{Number. Heterozygote genotype probability of the locus.}

\item{theta_hom}{Number. The dispersion parameter of the homozygous peaks.}

\item{theta_het}{Number. The dispersion parameter of the heterozygous peak.}
}
\value{
Probability of observing \code{ref_counts} and \code{var_counts}
}
\description{
\code{pmf_betabinomMix} calculates the probability of observing given population-level RNAseq data (i.e. both reference- and variant counts of one or more samples)
assuming a beta-binomial mixture model with parameter values as determined by the input. More specifically, the formula used is (using maelstRom's \code{dBetaBinom} function):

pr * dBetaBinom(ref_counts, ref_counts + var_counts, pi = 1 - SE, theta = theta_hom, LOG = FALSE) +

pv * dBetaBinom(var_counts, ref_counts + var_counts, pi = 1 - SE, theta = theta_hom, LOG = FALSE) +

prv * dBetaBinom(ref_counts, ref_counts + var_counts, pi = probshift, theta = theta_het, LOG = FALSE)
}
