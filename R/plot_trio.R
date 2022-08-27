#' Plot Parameters Trio to Observe the Two-Wing Pattern
#'
#' This function provides graphic display of the MLE parameter trio.
#'
#' @param param_trio  The n by 3 matrix of MLE parameter trio estimates.
#'
#' @return 
#' Nothing of interest.
#' 
#' @examples
#' # import data
#' data(hzmicrobe)
#' param_trio_bay <- calc_mle_trio(mic_bay, n_sample=10, replicates=3)
#' plot_trio(param_trio_bay)
#'
#'
#' @export
plot_trio <- function(param_trio){
  ggplot() +
    geom_point(data=param_trio, 
               aes(x=log(k), y=log(mu), color=pi0),alpha=0.6,size=0.24) +
    geom_vline(data=param_plt, aes(xintercept = mu_up), 
               linetype="dotdash", color = "black", size=0.12) +
    geom_vline(data=param_plt, aes(xintercept = k_low), 
               linetype="dotdash", color = "black", size=0.12) +
    ylab(TeX('\\log{(\\hat{\\mu})}')) +
    xlab(TeX('\\log{(\\hat{k})}')) +
    ggtitle(TeX('')) +
    theme_bw() +
    theme(aspect.ratio=1) +
    scale_color_continuous(type = "viridis") +
    facet_wrap(Procedure~., scale='free', nrow=2) +
    labs(color=TeX('$\\hat{\\pi}_0$')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 9, vjust = 0.5, hjust = 0),
          axis.text.y = element_text(size=9),
          strip.text = element_text(size = 12),
          legend.title = element_text(color = "black", size = 12),
          legend.text = element_text(color = "black", size = 9))
  
}


