#' Probability of observed outcomes across multiple households
#'
#' Probability is computed by simulating primary infections
#' T < 0: always iterate
#' T = 0: always simulate (M iterations)
#' T > 0: iterate for households up to size T, then simulate
#'
#' @param obs vector of observed cases. 1=positive for infection, 0=negative for infection, NA=unobserved individual
#' @param prim probabilities of primary infection
#' @param hh_sizes household sizes
#' @param q probability of secondary infection
#' @param Sn test sensitivity
#' @param Sp test specificity
#' @param Thresh household size threshold for use of simulation for primary infection
#' @param M number of iterations for primary infection simulation
#' @param it whether to use iteration for secondary cases (for performance measure only)
#' @return probability of observation
#' @export
prob_obs_hhs = function(obs,prim,hh_sizes,q,Sn=1,Sp=1,Thresh=-1,M=0,it=FALSE){
  # Recycle obs and prim
  df.rec = data.frame(
    obs=ifelse(is.na(obs),-1,as.integer(obs)),
    prim=as.numeric(prim)
  )
  obs = df.rec$obs
  prim = df.rec$prim

  # Recycle hh_sizes and q
  df.rec = data.frame(
    hh_sizes=as.integer(hh_sizes),
    q=as.numeric(q)
  )
  hh_sizes = df.rec$hh_sizes
  q = df.rec$q

  # Check that hh_sizes are correct
  if(length(obs) != sum(hh_sizes)) stop("Sum of hh sizes must match len of observations")

  Sn = as.numeric(Sn)
  Sp = as.numeric(Sp)
  Thresh = as.integer(Thresh)
  M = as.integer(M)

  if(it){
    .Rlprob_obs_hhs_it(obs,prim,hh_sizes,q,length(hh_sizes),Sn,Sp,Thresh,M)
  }else{
    .Rlprob_obs_hhs(obs,prim,hh_sizes,q,length(hh_sizes),Sn,Sp,Thresh,M)
  }
}

#' Probability of observed household outcome
#' Probability for primary infection is either computed by iteration or by simulation, depending on T.
#'
#' @param obs observed household cases. 1=infected, 0=not infection, NA=no observation
#' @param prim probabilities of primary infection
#' @param q probability of secondary infection
#' @param Sn test sensitivity
#' @param Sp test specificity
#' @param Thresh household size threshold for use of simulation for primary infection
#' @param M number of iterations for primary infection simulation
#' @param it whether to use iteration for secondary cases (for performance comparison only)
#' @return probability of observation
#' @export
prob_obs_hh = function(obs,prim,q,Sn=1,Sp=1,Thresh=-1,M=0,it=FALSE){
  df.rec = data.frame(
    obs = ifelse(is.na(obs),-1,as.integer(obs)),
    prim = as.numeric(prim)
  )
  obs = df.rec$obs
  prim = df.rec$prim
  hh_size = length(obs)

  q = as.numeric(q)
  Sn = as.numeric(Sn)
  Sp = as.numeric(Sp)
  Thresh = as.integer(Thresh)
  M = as.integer(M)

  if(it){
    .Rlprob_obs_it(obs,prim,hh_size,q,Sn,Sp,Thresh,M)
  }else{
    .Rlprob_obs(obs,prim,hh_size,q,Sn,Sp,Thresh,M)
  }
}

#' Probability of observed household outcome given primary cases
#' @param obs observed household cases (vector) NA=missed
#' @param pri primary household infections (vector)
#' @param q probability of secondary infection
#' @param Sn test sensitivity
#' @param Sp test specificity
#' @param it whether to use iteration for secondary cases
#' @return probability of observation
#' @export
lprob_obs_pri = function(obs,pri,q,Sn=1,Sp=1,it=FALSE){
  df.rec = data.frame(
    obs = ifelse(is.na(obs),-1,as.integer(obs)),
    pri = as.integer(pri)
  )
  obs = df.rec$obs
  pri = df.rec$pri
  hh_size = length(obs)

  q = as.numeric(q)
  Sn = as.numeric(Sn)
  Sp = as.numeric(Sp)
  if(it){
    .Rlprob_obs_pri_it(obs,pri,hh_size,q,Sn,Sp)
  }else{
    .Rlprob_obs_pri(obs,pri,hh_size,q,Sn,Sp)
  }
}

#' Probability of observed cases and missed observations among potential secondary infections given the number of primary infections.
#'
#' Secondary infections are computed using a generating functions and depend only on the number of secondary infections, independant of individuals infected.
#' Secondary infection probabilities are computed by iterating over possible numbers of true positives, and counting corresponding potential secondary infections
#'
#' @param obs number of observed positive infections among potential secondary infections.
#' @param mis number of missed observation among potential secondary infections.
#' @param pri number of primary infections, which are potential sources of secondary infection.
#' @param m number of potential secondary household infections.
#' @param q probability of secondary infection
#' @param Sn test sensitivity
#' @param Sp test specificity
#' @return probability of observation
#' @export
lprob_sec_obs_pri = function(obs,mis,pri,m,q,Sn=1,Sp=1){
  obs = as.integer(obs)
  mis = as.integer(mis)
  pri = as.integer(pri)
  m = as.integer(m)
  q = as.numeric(q)
  Sn = as.numeric(Sn)
  Sp = as.numeric(Sp)
  .Rlprob_sec_obs_pri(obs, mis, pri, m, q, Sn, Sp)
}

#' Probability of observed cases and missed observations among potential secondary infections given the number of primary infections, using iteration for secondary cases.
#'
#' Using an iterative method for secondary infections, for performance comparison only
#'
#' Secondary infections are computed using a generating functions and depend only on the number of secondary infections, independant of individuals infected.
#' Secondary infection probabilities are computed by iterating over possible numbers of true positives, and counting corresponding potential secondary infections
#'
#' @param obs vector of observations, -1=missed, 0=negative, 1=positive
#' @param mis number of missed observation among potential secondary infections.
#' @param pri number of primary infections, which are potential sources of secondary infection.
#' @param m number of potential secondary household infections.
#' @param q probability of secondary infection
#' @param Sn test sensitivity
#' @param Sp test specificity
#' @return probability of observation
#' @export
lprob_sec_obs_pri_it = function(obs,pri,m,q,Sn=1,Sp=1){
  obs = as.integer(obs)
  pri = as.integer(pri)
  m = as.integer(m)
  q = as.numeric(q)
  Sn = as.numeric(Sn)
  Sp = as.numeric(Sp)
  .Rlprob_sec_obs_pri_it(obs, k, m, q, Sn, Sp)
}
