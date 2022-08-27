#' Dispersion test
#'
#' This function examines the dispersion pattern.
#'
#' @param v  The m-length vector, denotying the abundance of a specific microbe.
#'
#' @return True or False.
#'
#' @examples
#' # import data
#' data(hzmicrobe)
#' is_dispersion(as.numeric(mic_bay[1,]))
#' is_dispersion(as.numeric(mic_era[12,]))
#'
#'
#' @export
#'
is_dispersion <- function(v){
  # evaluate dispersion
  # v is an m-length vector, denoting the abundance of a specific microbe
  return(var(v)>mean(v))
}


