% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{mcmc.marginal.lik}
\alias{mcmc.marginal.lik}
\title{MCMC Function method suitable for estimated likelihood function}
\usage{
mcmc.marginal.lik(posterior, ll, p0, sigma, lp, psuedo.ll = TRUE)
}
\arguments{
\item{posterior}{A matrix of posterior samples. Columns are parameters, rows are samples.}

\item{ll}{log likelihood function, taking as arguments parameters p0}

\item{p0}{any paramter value, one of high posterior density}

\item{sigma}{covariance matrix for multivariate gaussian update kernel}

\item{lp}{function of parameters p0 giving the prior log-probability for each parameter}

\item{psuedo.ll}{indicates if ll is psuedo-likelihood, in which case it must be recalculated each time}
}
\value{
log marginal likelihood estimate
}
\description{
MCMC Function method suitable for estimated likelihood function
}
