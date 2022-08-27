#' Plot Parameters Trio to Observe the Two-Wing Pattern
#'
#' This function provides graphic display of the MLE parameter trio.
#'
#' @param param_trio  The n by 3 matrix of MLE parameter trio estimates.
#' @param point_size  The numeric number between 0 and 1 for the point size. 
#' @param a The numeric number between 0 and 1 for the alpha transparency level.
#' @param zoom_in If TRUE, zoom in the majority of points. Default is FALSE.
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
plot_trio <- function(param_trio, point_size=0.8, a=0.6, zoom_in=FALSE){
  if(zoom_in){
    ggplot() +
      geom_point(data=param_trio, 
                 aes(x=log(k), y=log(mu), color=pi0),alpha=a,size=point_size) +
      ylab(TeX('\\log{(\\hat{\\mu})}')) +
      xlab(TeX('\\log{(\\hat{k})}')) +
      xlim(-5,20) +
      ylim(-5,10) +
      ggtitle(TeX('')) +
      theme_bw() +
      theme(aspect.ratio=1) +
      scale_color_continuous(type = "viridis") +
      labs(color=TeX('$\\hat{\\pi}_0$')) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12),
            axis.text.x = element_text(size = 9),
            axis.text.y = element_text(size=9),
            strip.text = element_text(size = 12),
            legend.title = element_text(color = "black", size = 12),
            legend.text = element_text(color = "black", size = 9))
  } else {
    ggplot() +
      geom_point(data=param_trio, 
                 aes(x=log(k), y=log(mu), color=pi0),alpha=a,size=point_size) +
      ylab(TeX('\\log{(\\hat{\\mu})}')) +
      xlab(TeX('\\log{(\\hat{k})}')) +
      ggtitle(TeX('')) +
      theme_bw() +
      theme(aspect.ratio=1) +
      scale_color_continuous(type = "viridis") +
      labs(color=TeX('$\\hat{\\pi}_0$')) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12),
            axis.text.x = element_text(size = 9),
            axis.text.y = element_text(size=9),
            strip.text = element_text(size = 12),
            legend.title = element_text(color = "black", size = 12),
            legend.text = element_text(color = "black", size = 9))
  }
  
  
}


