#' Probability of observed outcomes across multiple households
#'
#' Probability computer by iterating or simulating primary infections:
#' * T < 0: always iterate
#' * T = 0: always simulate (M iterations)
#' * T > 0: iterate for households up to size T, then simulate
#'
#' @param obs vector of observed cases. 1=positive for infection, 0=negative for infection, NA=unobserved individual
#' @param prim probabilities of primary infection
#' @param hh_sizes household sizes
#' @param q probability of secondary infection
#' @param Sn test sensitivity
#' @param Sp test specificity
#' @param Thresh household size threshold for use of simulation for primary infection
#' @param M number of iterations for primary infection simulation
#' @return log-probability of observation
#' @export
prob_obs_hhs = function(obs,prim,hh_sizes,q,Sn=1,Sp=1,Thresh=-1,M=0){
  obs = ifelse(is.na(obs),-1,as.integer(obs)) # Convert NA to -1
  prim = as.numeric(prim)
  hh_sizes = as.integer(hh_sizes)
  q=as.numeric(q)
  Sn = as.numeric(Sn)
  Sp = as.numeric(Sp)
  Thresh = as.integer(Thresh)
  M = as.integer(M)

  # Recycle obs and prim
  df.rec = data.frame(obs,prim)
  obs = df.rec$obs
  prim = df.rec$prim

  # Recycle hh_sizes and q
  df.rec = data.frame(hh_sizes,q)
  hh_sizes = df.rec$hh_sizes
  q = df.rec$q

  # Check that hh_sizes are consistent with observations
  if(length(obs) != sum(hh_sizes)) stop("Sum of hh sizes must match len of observations")

  .Rlprob_obs_hhs(obs,prim,hh_sizes,q,length(hh_sizes),Sn,Sp,Thresh,M)
}

#' Probability of observed household outcome
#'
#' Probability computer by iterating or simulating primary infections:
#' * T < 0: always iterate
#' * T = 0: always simulate (M iterations)
#' * T > 0: iterate for households up to size T, then simulate
#'
#' @param obs observed household cases. 1=infected, 0=not infection, NA=no observation
#' @param prim probabilities of primary infection
#' @param q probability of secondary infection
#' @param Sn test sensitivity
#' @param Sp test specificity
#' @param Thresh household size threshold for use of simulation for primary infection
#' @param M number of iterations for primary infection simulation
#' @return log-probability of observation
#' @export
prob_obs_hh = function(obs,prim,q,Sn=1,Sp=1,Thresh=-1,M=0){
  obs = ifelse(is.na(obs),-1,as.integer(obs)) # Convert NA to -1
  prim = as.numeric(prim)
  q=as.numeric(q)
  Sn = as.numeric(Sn)
  Sp = as.numeric(Sp)
  Thresh = as.integer(Thresh)
  M = as.integer(M)

  # Recycle obs and prim
  df.rec = data.frame(obs,prim)
  obs = df.rec$obs
  prim = df.rec$prim

  hh_size = length(obs)

  .Rlprob_obs(obs,prim,hh_size,q,Sn,Sp,Thresh,M)
}

#' Probability of observed household outcome given primary cases
#'
#' @param obs observed household cases (vector) NA=missed
#' @param pri primary household infections (vector). 0=uninfected, 1=infected
#' @param q probability of secondary infection
#' @param Sn test sensitivity
#' @param Sp test specificity
#' @return log-probability of observation
#' @export
prob_obs_pri = function(obs,pri,q,Sn=1,Sp=1){
  obs = ifelse(is.na(obs),-1,as.integer(obs)) # Convert NA to -1
  pri = as.integer(pri)
  q=as.numeric(q)
  Sn = as.numeric(Sn)
  Sp = as.numeric(Sp)

  # Recycle obs and prim
  df.rec = data.frame(obs,pri)
  obs = df.rec$obs
  pri = df.rec$pri

  hh_size = length(obs)

  .Rlprob_obs_pri(obs,pri,hh_size,q,Sn,Sp)
}

#' Probability of observation among potential secondary infections
#'
#' @param obs number of observed positive infections among potential secondary infections.
#' @param mis number of missed observation among potential secondary infections.
#' @param pri number of primary infections, which are potential sources of secondary infection.
#' @param m number of potential secondary household infections.
#' @param q probability of secondary infection
#' @param Sn test sensitivity
#' @param Sp test specificity
#' @return log-probability of observation
#' @export
prob_sec_obs_pri = function(obs,mis,pri,m,q,Sn=1,Sp=1){
  obs = as.integer(obs)
  mis = as.integer(mis)
  pri = as.integer(pri)
  m = as.integer(m)
  q = as.numeric(q)
  Sn = as.numeric(Sn)
  Sp = as.numeric(Sp)
  .Rlprob_sec_obs_pri(obs, mis, pri, m, q, Sn, Sp)
}
