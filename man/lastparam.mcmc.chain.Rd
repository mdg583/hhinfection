% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{lastparam.mcmc.chain}
\alias{lastparam.mcmc.chain}
\title{Get parameter values at the end of a chain}
\usage{
\method{lastparam}{mcmc.chain}(chain, n = 1)
}
\arguments{
\item{chain}{mcmc.chain}

\item{n}{number of end samples to average over}
}
\value{
parameter vector of average parameters over n last samples
}
\description{
If n is provides, parameter values for each chain are averaged over n final samples
}
