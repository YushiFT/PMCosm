#' Refine the boundary between sub-communities.
#'
#' This function computes a model-driven testable boundary between dispersal vanguards and laggards.
#' It assigns microbes in the temporary transition region into either mu-wing or k-wing.
#'
#' @param x  The $n$ by $m$ matrix of microbial abundance.
#'
#' @return A list of microbes' id for dispersal vanguards and laggards.
#'
#' @examples
#' # import data
#' data(hzmicrobe)
#' id_lis <- refine_mu_k(mic_bay)
#'
#' @export
#'
refine_mu_k <- function(x){
  id_mu <- c()
  id_k  <- c()
  lis_out <- list(id_mu, id_k)
  names(lis_out) <- c('mu','k')
  return(lis_out)
}



