% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{check.ll}
\alias{check.ll}
\title{Get a sample of log-likelihood estimates for given parameter}
\usage{
check.ll(p0, ll, runs = 200)
}
\arguments{
\item{p0}{parameter vector}

\item{ll}{log-likelihood estimate function}

\item{runs}{number of iterations}
}
\value{
samples of log-likelihood estimate
}
\description{
This is done in parallel if possible. It is for checking log-likelihood variance.
}
