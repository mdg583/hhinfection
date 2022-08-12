#' Precompute secondary infection probabilities
#'
#' Probabilities are generated for all household sizes from 1 to max_hh,
#' and for all possible numbers of primary infection from 0 to household size.
#'
#' Probabilities for a specific household are retrieved using sec_prob.
#'
#' @param max_hh Maximum household size for which to compute probabilities
#' @param sar Secondary attack rate (probability of infection between any infectious and susceptible)
#' @return list object containing pre-computed probabilities and max household size
#' @export
sec_prob_init = function(max_hh,sar){
  dat = .Rgdat(max_hh,sar)
  list(
    data=dat,
    max_hh=max_hh
  )
}

#' Get secondary household infection probabilities
#'
#' This uses pre-computed probabilities generated using sec_prob_innit
#'
#' @param sec_prob_data Maximum household size for which to compute probabilities
#' @param k Number of primary infections
#' @param m Number of susceptibles (household size - k)
#' @return vector of probabilities of infection for household, from 0 ... m
#' @export
sec_prob = function(sec_prob_data,k,m){
  .Rg(sec_prob_data$data, sec_prob_data$max_hh, k, m)
}
