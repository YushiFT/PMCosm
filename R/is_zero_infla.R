#' Zero-inflation test
#'
#' This function detects the existence of inflated zeros.
#'
#' @param v  An m-length vector, denoting the abundance of a specific microbe.
#' @param n_sample The numeric number of sample size.
#' @param replicates The numeric number of replicates per sample site.
#' @param max_prop_zero The maximum acceptable value of zero proportions per sample site.
#'
#'
#' @return True or False.
#'
#' @examples
#' # import data
#' data(hzmicrobe)
#' is_zero_infla(as.numeric(mic_bay[1,]), n_sample=10, replicates=3)
#' is_zero_infla(as.numeric(mic_era[12,]), n_sample=12, replicates=3)
#'
#'
#' @export
#'
is_zero_infla <- function(v, n_sample=NULL, replicates=NULL, max_prop_zero=2/3){
  # evaluate zero inflation
  # v is an m-length vector, denoting the abundance of a specific microbe
  replicates <- ifelse(is.null(replicates), 3, replicates)
  n_sample <- ifelse(is.null(n_sample), length(v)/replicates, n_sample)
  record <- FALSE
  max_zero <- replicates * max_prop_zero
  sudo_m <- diag(n_sample)
  for(k in 1:ncol(sudo_m)){
    sudo_v <- rep(sudo_m[,k], each=replicates)
    num_0  <- sum((v * sudo_v) == 0) - sum(sudo_v==0)
    record <- ifelse(num_0>max_zero,TRUE,record)
  }
  return(record)
}



