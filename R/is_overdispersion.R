#' Over-dispersion test
#'
#' This function diagnose the over-dispersion generally for microbial abundance data.
#'
#' @param v  An m-length vector, denoting the abundance of a specific microbe.
#'
#' @return True or False.
#'
#' @examples
#' # import data
#' data(hzmicrobe)
#' is_overdispersion(as.numeric(mic_bay[1,]))
#' is_overdispersion(as.numeric(mic_era[12,]))
#'
#'
#' @export
is_overdispersion <- function(v, alpha=0.05){
  test_dat <- data.frame(count = as.numeric(v))
  fit_test <- glm(count~., data=test_dat, family=poisson)
  disp_test <- AER::dispersiontest(fit_test, trafo=2, alternative='greater')
  return(disp_test$p.value <= alpha)
}



