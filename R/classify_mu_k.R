#' Refine the boundary between sub-communities.
#'
#' This function computes a model-driven testable boundary between dispersal vanguards and laggards.
#' It assigns every microbe within the whole community into either mu-wing (i.e., the sub-community of dispersal vanguards) 
#' or k-wing (i.e., the sub-community of dispersal laggards). This is the same function as
#' classify_vag_lag. This function is for paper-peer-review only.
#'
#' @param x  The n by m matrix of microbial abundance.
#' @param mle_trio Then n by 3 matrix of MLE estimates.
#'
#' @return A list of microbes' id for dispersal vanguards and laggards.
#'
#' @examples
#' # import data
#' data(hzmicrobe)
#' param_trio_bay <- calc_mle_trio(mic_bay, n_sample=10, replicates=3)
#' id_vag_lag <- classify_vag_lag(mic_bay, param_trio_bay)
#' # print the ids for dispersal vanguards
#' print(id_vag_lag$vanguards)
#' # print the ids for dispersal laggards
#' print(id_vag_lag$laggards)
#'
#' @export
classify_mu_k <- function(x, mle_trio){
  # coarse over-dispersion group and non-over-dispersion group
  od_test <- apply(x, 1, is_overdispersion)
  id_od  <- names(od_test)[od_test==TRUE]
  id_nod <- names(od_test)[od_test==FALSE]
  k_od  <- mle_trio[id_od,]$k
  k_nod <- mle_trio[id_nod,]$k
  l_k_od  <- log(k_od[k_od>0])
  l_k_nod <- log(k_nod[k_nod>0])
  
  calc_mad <- function(v){
    # calculate median absolute deviation
    return(median(abs(v-median(v))))
  }
  
  # b1: the lower bound for the temporary transition region
  b1_logk <- median(l_k_od) + 1.4826*calc_mad(l_k_od)
  # b2: the upper bound for the temporary transition region
  b2_logk <- median(l_k_nod) - 1.4826*calc_mad(l_k_nod)
  # assign microbes with log(k) less than b1 to mu-wing, i.e., dispersal vanguards
  id_mu <- rownames(mle_trio[mle_trio$k<exp(b1_logk),])
  # assign microbes with log(k) greater than b2 to k-wing, i.e., dispersal laggards
  id_k <- rownames(mle_trio[mle_trio$k>exp(b2_logk),])
  
  # assign microbes with log(k) between b1 and b2 to the temporary transition zone (ttz)
  id_ttz <- rownames(mle_trio[(mle_trio$k>=exp(b1_logk))&(mle_trio$k<=exp(b2_logk)),])
  mic_ttz <- x[id_ttz,]
  
  # refine the sub-community boundary among ttz
  ttz_vag_lag <- classify_ttz(mic_ttz)
  
  id_mu <- c(id_mu, ttz_vag_lag$vanguards)
  id_k  <- c(id_k,  ttz_vag_lag$laggards)
  
  lis_out <- list(id_mu, id_k)
  names(lis_out) <- c('vanguards','laggards')
  return(lis_out)
}





