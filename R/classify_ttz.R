#' Identify Dispersal Vanguards and Laggards in the Temporary Transition Zone
#'
#' This function computes a model-driven testable boundary between dispersal vanguards and laggards
#' in the Temporary Transition Zone (TTZ). It assigns every microbe within TTZ into either mu-wing 
#' (i.e., the sub-community of dispersal vanguards) or k-wing (i.e., the sub-community of dispersal laggards).
#'
#' @param x  The matrix of microbial abundance.
#'
#' @return A list of microbes' id for dispersal vanguards and laggards.
#'
#' @examples
#' # import data
#' data(hzmicrobe)
#' param_trio_bay <- calc_mle_trio(mic_bay, n_sample=10, replicates=3)
#' # classify_ttz is wrapped up in function classify_vag_lag
#' id_vag_lag <- classify_vag_lag(mic_bay, param_trio_bay)
#' # print the ids for dispersal vanguards
#' print(id_vag_lag$vanguards)
#' # print the ids for dispersal laggards
#' print(id_vag_lag$laggards)
#'
#' @export
classify_ttz <- function(x){
  # loglikelihood function for poisson model
  calc_ll_poisson <- function(v, theta){
    return(sum(- theta + v * log(theta) - lgamma(v + 1)))
  }
  id_mu <- c()
  id_k  <- c()
  for(i in 1:nrow(x)){
    test_dat <- data.frame(count = as.numeric(x[i,]))
    # over-dispersion test
    # to mu-wing if observing significant over-dispersion
    fit_p <- glm(count~., data=test_dat, family=poisson)
    disp_test <- AER::dispersiontest(fit_p, trafo=2, alternative="greater")
    p_od  <- disp_test$p.value
    if(p_od <= 0.05){
      id_mu <- c(id_mu, rownames(x)[i])
    } else {
      # conduct likelihood-ratio-test
      if (is_zero_infla(test_dat$count)){
        # fit zero-inflated gamma poisson
        fit_zgp <- pscl::zeroinfl(count~1|1, link='logit', dist='negbin',data=test_dat)
        ll_zgp <- fit_zgp$loglik
        # fit zero-inflated poisson
        fit_zp <- pscl::zeroinfl(count~1|1, link='logit', dist='poisson',data=test_dat)
        ll_zp <- fit_zp$loglik
        # h0 prefer zero-inflated poisson
        # h0 prefer k-wing
        # mu-wing if p_value less than or equal to 0.05
        # k-wing if p_value greater than 0.05
        p_lrt <- pchisq(2*(ll_zgp-ll_zp), df = 3-2, lower.tail = FALSE)
        if(p_lrt <= 0.05) {
          id_mu <- c(id_mu, rownames(x)[i])
        } else {
          id_k  <- c(id_k,  rownames(x)[i])
        }
      } else {
        # fit gamma poisson
        fit_gp <- VGAM::vglm(count~1, VGAM::negbinomial(zero=-2), data=test_dat)
        # gamma poisson log-likelihood
        ll_gp <- fit_gp@criterion$loglikelihood
        # fit poisson
        fit_p <- glm(count~1, data=test_dat, family=poisson(), trace=FALSE)
        # poisson log-likelihood
        mu  <- exp(fit_p$coefficients[['(Intercept)']])
        ll_p <- calc_ll_poisson(test_dat$count, theta=mu)
        # h0 prefer poisson
        # h0 prefer k-wing
        # mu-wing if p_value less than or equal to 0.05
        # k-wing if p_value greater than 0.05
        p_lrt <- pchisq(2*(ll_gp-ll_p), df = 2-1, lower.tail = FALSE)
        if(p_lrt <= 0.05) {
          id_mu <- c(id_mu, rownames(x)[i])
        } else {
          id_k  <- c(id_k,  rownames(x)[i])
        }
      }
    }
  }
  lis_out <- list(id_mu, id_k)
  names(lis_out) <- c('vanguards','laggards')
  return(lis_out)
}



