% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{plot.mcmc.chains}
\alias{plot.mcmc.chains}
\title{Create a plot for a set of MCMC chains}
\usage{
\method{plot}{mcmc.chains}(chains, type = "trace", ...)
}
\arguments{
\item{chains}{mcmc.chains}

\item{type}{type of plot to plot}

\item{...}{arguments for plotting}
}
\description{
Plotting uses ggplot and the functions provided by bayesplot: \code{\link[bayesplot:MCMC-traces]{bayesplot::mcmc_trace_data()}}
\code{\link[bayesplot:MCMC-diagnostics]{bayesplot::mcmc_acf()}} \code{\link[bayesplot:MCMC-intervals]{bayesplot::mcmc_intervals()}} \code{\link[bayesplot:MCMC-distributions]{bayesplot::mcmc_violin()}} \code{\link[bayesplot:MCMC-distributions]{bayesplot::mcmc_dens()}}
}
\details{
Plot options:
\itemize{
\item trace: trace of parameters for chains
\item likelihood: trace of likelihood for chains
\item acf: autocorrelation plots for parameters
\item intervals: intervals plots for parameters
\item violins: violins plots for parameters
\item density: density plots for parameters
\item desnity_2d: Two-dimensional density plot for 2 parameters
}
}
