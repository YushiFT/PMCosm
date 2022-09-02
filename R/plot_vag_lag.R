#' Visualize the Two-Wing admixed structure
#'
#' This function provides graphic display of the admixed structure of dispersal vanguards and laggards.
#'
#' @param param_trio  The n by 3 matrix of MLE parameter trio estimates.
#' @param point_size  The numeric number for the point size. 
#' @param a The numeric number between 0 and 1 for the alpha transparency level.
#' @param zoom_in If TRUE, zoom in the majority of points. Default is FALSE.
#' @param show_inf If TRUE, list Poisson distributed data -- i.e., microbes with infinite k -- on the right. Defalt is FALSE.
#'
#' @return 
#' Nothing of interest.
#' 
#' @examples
#' # import data
#' data(hzmicrobe)
#' param_trio_bay <- calc_mle_trio(mic_bay, n_sample=10, replicates=3)
#' id_vag_lag <- classify_vag_lag(mic_bay, param_trio_bay)
#' plot_vag_lag(param_trio_bay, id_vag_lag)
#'
#'
#' @export
plot_vag_lag <- function(param_trio, id_category,
                         point_size=0.8, a=0.6, zoom_in=FALSE, show_inf=FALSE){
  if(!show_inf) {
    param_trio <- param_trio[param_trio$k!=Inf,]
  }
  
  param_trio$Category <- "None"
  param_trio[id_category$vanguards,]$Category <- "Dispersal Vanguards"
  param_trio[id_category$laggards,]$Category <- "Dispersal Laggards"
  param_trio$Category <- factor(param_trio$Category,
                                levels=c("Dispersal Vanguards",
                                         "Dispersal Laggards"))
  if(zoom_in){
    ggplot() +
      geom_point(data=param_trio, 
                 aes(x=log(k), y=log(mu), color=Category),alpha=a,size=point_size) +
      ylab(TeX('\\log{(\\hat{\\mu})}')) +
      xlab(TeX('\\log{(\\hat{k})}')) +
      xlim(-5,20) +
      ggtitle(TeX('')) +
      theme_bw() +
      theme(aspect.ratio=1) +
      scale_colour_manual(values = c('#984ea3','#4daf4a')) +
      labs(color=TeX('Category')) +
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
                 aes(x=log(k), y=log(mu), color=Category),alpha=a,size=point_size) +
      ylab(TeX('\\log{(\\hat{\\mu})}')) +
      xlab(TeX('\\log{(\\hat{k})}')) +
      ggtitle(TeX('')) +
      theme_bw() +
      theme(aspect.ratio=1) +
      scale_colour_manual(values = c('#984ea3','#4daf4a')) +
      labs(color=TeX('Category')) +
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


