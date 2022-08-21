#' MCMC Function method suitable for psuedo-likelihood
#'
#' @param n length of chain
#' @param ll log likelihood function, taking as arguments parameters p0
#' @param p0 initial value of parameters
#' @param sigma covariance matrix for multivariate gaussian update kernel
#' @param lp function of parameters p0 giving the prior log-probability for each parameter
#' @param psuedo.ll indicates if ll is psuedo-likelihood, in which case it must be recalculated each time
#' @return mcmc.chain object, containing chain samples, acceptance variables and computed log-likelihoods
#' @export
#' @importFrom MASS mvrnorm
mcmc.chain = function(n,ll,p0,sigma,lp,psuedo.ll=TRUE){
  updates = floor(n/10) # for console output
  chain = matrix(NA,nrow=n,ncol=length(p0))
  colnames(chain) = names(p0)
  ac = rep(NA,n) # acceptance vector, for acceptance rate
  li = rep(NA,n) # estimated log-likelihoods during chain construction
  l0 = ll(p0) # initial log-likelihood estimate
  for(i in 1:n){
    delta.p = mvrnorm(1,rep(0,length(p0)),sigma)
    p1 = p0 + delta.p
    l1 = ll(p1) # New log-likelihood
    # If using psuedo-likelihood, liklihood function needs to be re-calculated even if parameters unchanged
    if(psuedo.ll){
      a = l1 - ll(p0) + sum(lp(p1) - lp(p0))
    }else{
      a = l1 - l0 + sum(lp(p1) - lp(p0))
    }
    # I think nan happens for Inf-Inf (both likelihoods 0)
    ac[i] = if(is.nan(a)) FALSE else a >= 0 || rbinom(1,1,exp(a))==1
    if(ac[i]){
      p0 = p1
      l0 = l1
    }
    chain[i,] = p0
    li[i] = l0
    if(i %% updates == 0){
      system(sprintf('echo "%s"', paste0(i," : ", round(l0,3))))
    }
  }
  rval = list(chain=chain,ac=ac,li=li)
  class(rval) = "mcmc.chain"
  rval
}

#' Build a set of MCMC chains
#'
#' @param num.chains number of chains to build
#' @param n length of chains
#' @param ll log likelihood function, taking as arguments parameters p0
#' @param p0.list list, initial value of parameters for each chain
#' @param sigma covariance matrix for multivariate gaussian update kernel
#' @param lp function of parameters p0 giving the prior log-probability for each parameter
#' @param psuedo.ll indicates if ll is psuedo-likelihood, in which case it must be recalculated each time
#' @return mcmc.chains object, list of chain objects
#' @export
#' @importFrom MASS mvrnorm
mcmc.chains = function(num.chains,n,ll,p0.list,sigma,lp,psuedo.ll=TRUE){
  require("parallel")
  cores = getOption("mc.cores")
  tim = system.time({
    chains = mclapply(1:num.chains,function(i){
      mcmc.chain(n,ll,p0.list[[i]],sigma,lp,psuedo.ll)
    })
  })
  message(paste0("Time to build chains: ", tim["elapsed"]))
  class(chains) = c(class(chains),"mcmc.chains")
  chains
}

#' MCMC Function method suitable for estimated likelihood function
#'
#' @param posterior A matrix of posterior samples. Columns are parameters, rows are samples.
#' @param ll log likelihood function, taking as arguments parameters p0
#' @param p0 any paramter value, one of high posterior density
#' @param sigma covariance matrix for multivariate gaussian update kernel
#' @param lp function of parameters p0 giving the prior log-probability for each parameter
#' @param psuedo.ll indicates if ll is psuedo-likelihood, in which case it must be recalculated each time
#' @return log marginal likelihood estimate
#' @export
#' @importFrom matrixStats logSumExp
#' @importFrom mvtnorm rmvnorm dmvnorm
mcmc.marginal.lik = function(posterior,ll,p0,sigma,lp,psuedo.ll=TRUE){
  # Log prob density of multivariate norm, with mean 0 and given sigma
  lg.dmvnorm = function(p0, sigma){
    dmvnorm(p0,mean=rep(0,length(p0)),sigma=sigma,log=TRUE)
  }

  library(parallel)
  cores = getOption("mc.cores")
  if(is.null(cores)) cores = 1

  # Skip log-likelihood of p0 each time if not using psuedo-likelihood
  if(!psuedo.ll){
    lp0 = ll(p0)
    ll.p0 = function() lp0
  }else{
    ll.p0 = function() ll(p0)
  }

  # Numerator
  updates = floor(nrow(posterior)/10) # for console output
  f.numer = function(i){
    if(updates && i %% updates == 0){
      system(sprintf('echo "num %d"', i))
    }
    p.i = as.numeric(posterior[i,])
    a = min(0,ll.p0() - ll(p.i) + sum(lp(p0) - lp(p.i))) + lg.dmvnorm(p0-p.i,sigma)
    if(is.nan(a)) -Inf else a
  }

  # Generate samples for denomenator
  prop.s = rmvnorm(nrow(posterior),mean=p0,sigma=sigma)

  # Denominator
  f.denom = function(i){
    if(updates && i %% updates == 0){
      system(sprintf('echo "den %d"', i))
    }
    p.i = as.numeric(prop.s[i,])
    a = min(0,ll(p.i) - ll.p0() + sum(lp(p.i) - lp(p0)))
    if(is.nan(a)) -Inf else a
  }

  t.num = system.time({
    x.numer = if(cores > 1){
      x.1 = do.call(c, mclapply(1:nrow(posterior), f.numer))
      logSumExp(x.1)
    }else{
      x.1 = vapply(1:nrow(posterior),f.numer,0)
      logSumExp(x.1)
    }
  })["elapsed"]

  message(paste0("Time to generate numerator: ", round(t.num,4)))

  t.den = system.time({
    x.denom = if(cores > 1){
      x.2 = do.call(c, mclapply(1:nrow(posterior), f.denom))
      logSumExp(x.2)
    }else{
      x.2 = vapply(1:nrow(posterior),f.denom,0)
      logSumExp(x.2)
    }
  })["elapsed"]
  message(paste0("Time to generate denominator: ", round(t.den,4)))

  f.mle = function(i){
    if(updates && i %% updates == 0){
      system(sprintf('echo "mle %d"', i))
    }
    ll(p0)
  }

  t.mle = system.time({
    x.mle = if(cores > 1){
      x.3 = do.call(c, mclapply(1:nrow(posterior), f.mle))
      # sum(x.3) / nrow(posterior)
      logSumExp(x.3) - log(nrow(posterior))
    }else{
      x.3 = vapply(1:nrow(posterior),f.mle,0)
      # sum(x.3) / nrow(posterior)
      logSumExp(x.3) - log(nrow(posterior))
    }
  })["elapsed"]
  message(paste0("Time to generate ll(mle): ", round(t.mle,4)))

  # Estimated marginal likelihood
  x.mle + sum(lp(p0)) - (x.numer - x.denom)
}

#' @export
lastparam = function(chain,n){
  UseMethod("lastparam")
}

#' Get parameter values at the end of a chain
#'
#' If n is provides, parameter values for each chain are averaged over n final samples
#'
#' @param chain mcmc.chain
#' @param n number of end samples to average over
#' @return parameter vector of average parameters over n last samples
#' @export
lastparam.mcmc.chain = function(chain,n=1){
  l = nrow(chain$chain)
  apply(chain$chain[(l+1-n):l,,drop=FALSE],2,mean)
}

#' Get parameter values at the end of each chain in mcmc.chains object
#'
#' If n is provides, parameter values for each chain are averaged over n final samples
#'
#' @param chains mcmc.chains
#' @param n number of end samples to average over
#' @return list of parameter vectors of average parameters over n last samples
#' @export
lastparam.mcmc.chains = function(chains,n=1){
  lapply(chains,function(chain){
    lastparam(chain,n)
  })
}

#' Get all chain samples from the chains in an mcmc.chains objects as a matrix
#'
#' @param chains mcmc.chains
#' @return a matrix with each column a parameter
combine.chains = function(chains){
  chains_list = lapply(chains,function(ch) ch$chain)
  do.call(rbind,chains_list)
}

#' Get all samples for a given parameter across chains
#'
#' @param chains mcmc.chains object
#' @param param name of parameter to get samples for
chains_param = function(chains,param){
  do.call(cbind,lapply(chains,function(ch) ch$chain[,param]))
}

#' Create a new sigma based on parameter chain
#'
#' Update uses cov(chain) * 2.38^2 / k for k parameters
#'
#' @param chains An mcmc.chain or mcmc.chains
#' @return A suggested sigma function for multivariate normal proposal function
#' @export
chains.sigma = function(chains){
  if("mcmc.chains" %in% class(chains)){
    param.chain = combine.chains(chains)
  }else{
    param.chain = chain$chain
  }
  cov(param.chain)*(2.38^2)/ncol(param.chain)
}

#' Print convergence stats for MCMC chain
#'
#' @param chain mcmc.chain
#' @export
print.mcmc.chain = function(chain){
  ch.len = nrow(chain$chain)
  cat(paste0("MCMC chain of length ",ch.len,".\n"))
  cat("\n")
  cat(paste0("Acceptance Rate: ", 100*round(mean(chain$ac),3), "%\n"))
  cat("Final Parameters:\n")
  print(as.data.frame(chain$chain[ch.len,]))
  cat("\n")
}

#' Print convergence stats set of for MCMC chains
#'
#' Includes convergence stats R-hat and effective sample size (block and tail),
#' as computed by rstan functions [rstan::Rhat()] [rstan::ess_bulk()] [rstan::ess_tail()].
#'
#' @param chains mcmc.chains
#' @export
#' @importFrom rstan Rhat ess_bulk ess_tail
print.mcmc.chains = function(chains){
  n = length(chains)
  ch.len = nrow(chains[[1]]$chain)
  cat(paste0("Set of ", n, " MCMC chains of length ",ch.len,".\n"))
  cat("\n")
  cat("Acceptance Rates\n")
  acc.df = as.data.frame(do.call(cbind,lapply(chains, function(chain){
    mean(chain$ac)
  })))
  colnames(acc.df) = 1:length(chains)
  print(acc.df,row.names=FALSE)
  cat("\n")
  # cat("Initial Parameters\n")
  # print(as.data.frame(do.call(rbind, lapply(chains, function(chain){
  #   chain$chain[1,]
  # }))))
  cat("Final Parameters\n")
  pf.df = as.data.frame(do.call(rbind, lapply(chains, function(chain){
    chain$chain[ch.len,]
  })))
  print(pf.df)
  cat("\n")

  cat("Evidence of Convergence\n")
  param.names = colnames(chains[[1]]$chain)
  cv.df = do.call(rbind,lapply(param.names,function(param){
    ch.p = chains_param(chains, param)
    rhat = Rhat(ch.p)
    ess.b = ess_bulk(ch.p)
    ess.t = ess_tail(ch.p)
    if(rhat < 1.01 && ess.b > 400 && ess.t > 400){
      conv="*"
    }else{
      conv=""
    }
    data.frame(R.hat=rhat,ess.b=ess.b,ess.t=ess.t,conv=conv)
  }))
  colnames(cv.df) = c("R.hat","ess.b","ess.t","")
  rownames(cv.df) = param.names
  print(cv.df)
}

#' Create a plot for a set of MCMC chains
#'
#' Plotting uses ggplot and the functions provided by bayesplot: [bayesplot::mcmc_trace_data()]
#' [bayesplot::mcmc_acf()] [bayesplot::mcmc_intervals()] [bayesplot::mcmc_violin()] [bayesplot::mcmc_dens()]
#'
#' Plot options:
#' - trace: trace of parameters for chains
#' - likelihood: trace of likelihood for chains
#' - acf: autocorrelation plots for parameters
#' - intervals: intervals plots for parameters
#' - violins: violins plots for parameters
#' - density: density plots for parameters
#' - desnity_2d: Two-dimensional density plot for 2 parameters
#'
#' @param chains mcmc.chains
#' @param type type of plot to plot
#' @param ... arguments for plotting
#' @export
#' @importFrom bayesplot mcmc_trace_data mcmc_acf mcmc_intervals mcmc_violin mcmc_dens
#' @importFrom dplyr `%>%`
#' @import ggplot2
plot.mcmc.chains = function(chains, type="trace",...){
  chains_list = lapply(chains,function(ch) ch$chain)
  if(type=="trace"){
    mcmc_trace_data(chains_list) %>%
      ggplot(aes(x=iteration,y=value)) +
      geom_line(size=.1) +
      facet_grid(parameter~chain) +
      labs(title="Parameter Trace")
  }else if(type=="likelihood"){
    lik.df = as.data.frame(do.call(rbind,lapply(1:length(chains), function(i){
      data.frame(
        chain=i,
        i=1:length(chains[[i]]$li),
        li=chains[[i]]$li
      )
    })))
    ggplot(lik.df,aes(x=i,y=li)) +
      geom_line() +
      facet_wrap(~chain) +
      labs(title="likelihood")
  }else if(type=="acf"){
    mcmc_acf(chains_list,...) +
      facet_grid(Parameter ~ Chain) + labs(title="Parameter ACF")
  }else if(type=="intervals"){
    mcmc_intervals(chains_list,...) + labs(title="Parameter Intervals")
  }else if(type=="violins"){
    mcmc_violin(chains_list,...)
  }else if(type=="density"){
    mcmc_dens(chains_list,...)
  }else if(type=="density_2d"){
    combined = as.data.frame(do.call(rbind, chains_list))
    fargs = list(...)
    if("pars" %in% names(fargs)){
      pars = fargs[["pars"]]
      fargs["pars"] = NULL
    }else{
      pars = colnames(combined)[1:2]
    }
    if("transformations" %in% names(fargs) && is.function(fargs[["transformations"]])){
      combined = as.data.frame(apply(combined, c(1,2), fargs[["transformations"]]))
    }
    ggplot(combined,aes_string(x=pars[[1]],y=pars[[2]])) +
      stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white")
  }
}

#' Initialize covariance matrix of proposal function
#'
#' Covariance matrix is diagonal with standard deviation 1/10 sd of priors
#'
#' @param p.var vector of variances of parameter priors
#' @return a nxn diagonal matrix
#' @export
init.sig = function(p.var){
  # We set the initial proposal sd as 1/10 the sd of the priors
  sqrt(diag(p.var))/10
}

#' Get initial parameter values that give a finite likelihood
#'
#' @param param.names names of parameters (also indicates number of parameters)
#' @param ll log-likelihood estimate function
#' @param chains number of chains.
#' @return If a non-0 likelihood is found within 1000 uniform random samples, the corresponding parameter vector is returned
#' @export
init.params = function(param.names,ll,chains=1,prior){
  if(chains > 1){
    lapply(1:chains, function(i){init.params(param.names,ll,1,prior)})
  }else{
    nvars = length(param.names)
    p0 = prior()
    tries = 0
    inf.nan = function(x) any(is.nan(x) | is.infinite(x))
    while(inf.nan(ll(p0))){
      tries = tries + 1
      p0 = prior()
      if(tries > 1000) stop("Failed to find parameter values in 1000 iterations")
    }
    names(p0) = param.names
    p0
  }
}

#' Get a sample of log-likelihood estimates for given parameter
#'
#' This is done in parallel if possible. It is for checking log-likelihood variance.
#'
#' @param p0 parameter vector
#' @param ll log-likelihood estimate function
#' @param runs number of iterations
#' @return samples of log-likelihood estimate
#' @export
#' @import parallel
#' @importFrom dplyr `%>%`
check.ll = function(p0,ll,runs=200){
  # Function to cut a vector into list items, for assigning to threads. cut.list(1:5,4)
  cut.list = function(x,...){
    y = cut(x,...)
    lapply(levels(y),function(fct) x[y==fct])
  }
  cores = getOption("mc.cores")
  if(is.null(cores)) cores = 1
  if(cores > 1 && runs > cores*10){
    tmp.x = mclapply(cut.list(1:runs,cores), function(run.i){
      vapply(run.i,function(i){
        ll(p0)
      },0)
    }) %>% do.call(c,.)
  }else{
    tmp.x = vapply(1:runs,function(i){
      if(i %% floor(runs/4) == 0) message(i)
      ll(p0)
    },0)
  }
  tmp.x[!is.infinite(tmp.x)] #%>% var()
}
