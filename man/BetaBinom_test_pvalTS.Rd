% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BetaBinom_test_pvalTS.R
\name{BetaBinom_test_pvalTS}
\alias{BetaBinom_test_pvalTS}
\title{Performs an exact beta-binomial p-test}
\usage{
BetaBinom_test_pvalTS(m, n, pi, theta)
}
\arguments{
\item{m}{Number. Number of successes}

\item{n}{Number. Number of trials}

\item{pi}{Number. Probability of success; \code{0 >= pi >= 1}}

\item{theta}{Number. Overdispersion parameter; \code{0 >= theta > +Inf}}
}
\value{
The probability to make an observation equally or less likely than the input data
}
\description{
\code{BetaBinom_test_pvalTS} performs an exact test,
returning a p-value that reflects the probability of the observed data \code{m}
assuming a beta-binomial PMF with given parameters \code{m}, \code{pi} and \code{theta}.
See \code{dBetaBinom} for more information on the used parameterization.
}
