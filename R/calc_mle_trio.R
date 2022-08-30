#' MLE Parameters Trio
#'
#' This function computes the MLE of prior parameters
#' in the Gamma-Poisson model or Zero-inflated Gamma-Poisson model.
#'
#' @param x  The n by m matrix of microbial abundance.
#' @param n_sample The numeric number of sample size.
#' @param replicates The numeric number of replicates per sample site.
#' @param max_prop_zero The maximum acceptable value of zero proportions per sample site.
#' @param trace If TRUE, print out optimization details. Default is FALSE.
#' @param nbinomial_source The source package of fitting Negative Binomial distribution. Can be either VGAM or MASS. Default is MASS. 
#'
#' @return A trio of two estimates of gamma-poisson parameters and the coefficient of zero-inflated model.
#'
#' @examples
#' # import data
#' data(hzmicrobe)
#' param_bay <- calc_mle_trio(mic_bay, n_sample=10, replicates=3)
#'
#'
#' @export
#'
calc_mle_trio <- function(x, n_sample=NULL, replicates=NULL, max_prop_zero=2/3, trace=FALSE,
                          nbinomial_source='MASS'){
  param_trio <- data.frame()
  oldw <- getOption("warn")
  if(!trace){options(warn=-1)}
  for(i in 1:nrow(x)){
    if (is_zero_infla(as.numeric(x[i,]), n_sample, replicates, max_prop_zero=2/3)){
      # fit zero-inflated model
      test_dat <- data.frame(count = as.numeric(x[i,]))
      tryCatch({
        fit_test <- pscl::zeroinfl(count~1|1, link='logit', dist='negbin',data=test_dat)
        mu  <- exp(fit_test$coefficients$count[['(Intercept)']])
        k   <- fit_test$theta
        logit_pi <- fit_test$coefficients$zero[['(Intercept)']]
        pi0 <- exp(logit_pi)/(1+exp(logit_pi))
      }, error=function(e){
        cat("Fine-tuning recommended for microbe ", rownames(x)[i], "\n")
      })
    } else if (is_overdispersion(as.numeric(x[i,]), alpha=0.05)) {
      # fit negative binomial
      test_dat <- data.frame(count = as.numeric(x[i,]))
      tryCatch({
        if(nbinomial_source=='VGAM'){
          fit_test <- VGAM::vglm(count~1, negbinomial(deviance=TRUE), data=test_dat)
          mu  <- exp(fit_test@coefficients[['(Intercept):1']])
          k   <- exp(fit_test@coefficients[['(Intercept):2']])
          pi0 <- 0
        } else {
          fit_test <- MASS::glm.nb(count~1, data=test_dat,
                                   control=glm.control(maxit=100,trace=FALSE))
          mu  <- exp(fit_test$coefficients[['(Intercept)']])
          k   <- fit_test$theta
          pi0 <- 0
        }
      }, error=function(e){
        cat("Fine-tuning recommended for microbe ", rownames(x)[i], "\n")
      })
    } else {
      # fit poisson
      test_dat <- data.frame(count = as.numeric(x[i,]))
      fit_test <- glm(count~1, data=test_dat, family=poisson(), trace=FALSE)
      mu  <- exp(fit_test$coefficients[['(Intercept)']])
      k   <- Inf
      pi0 <- 0
    }
    param_trio <- rbind(param_trio,
                        data.frame(mu  = mu,
                                   k   = k,
                                   pi0 = pi0))
    print(paste0('process ',round(i/nrow(x)*100, digits=2),'%' ))
  }
  options(warn=0)
  rownames(param_trio) <- rownames(x)
  print('completed!!!')
  return(param_trio)
}





