#' plot_proportion_of_variance
#'
#' performs pca on the sample raw matrix
#'
#' @param PCAObj_Summary A list of summary of prcomp function.
#' @return A plotly object
#' @examples 
#' plot_proportion_of_variance(PCAObj_Summary)
#' @export
plot_proportion_of_variance <- function(PCAObj_Summary){
  require(ggplot2)
  require(plotly)
  
  pc_prop_var_df <- PCAObj_Summary$importance[2,]
  data <- data.frame(pc= (names(pc_prop_var_df)), variance=as.numeric(pc_prop_var_df)*100)
  
  p <- ggplot(data, aes(x=pc, y=variance, text=pc, fill=pc)) + 
    labs(title="Proportion of Variance", x="Principal Component", y="Variance (%)", fill='PC') +
    geom_bar(stat = "identity") + scale_x_discrete(limits=data$pc) + coord_cartesian(ylim=c(0,100))+
    scale_y_continuous(breaks = pretty(0:100, n=10))+
    theme(plot.title = element_text(hjust=0.5, size=18, face = "bold"),
          legend.position = "none", legend.direction = "none", # legend positioned at the bottom, horizantal direction,
          axis.line = element_line(size=1, colour = "black"),	# axis line of size 1 inch in black color
          panel.grid.major = element_blank(),	# major grids included
          panel.grid.minor = element_blank(),	# no minor grids
          panel.border = element_blank(), panel.background = element_blank(), # no borders and background color
          axis.title = element_text(colour="black", size = 15, face = "bold"), # axis title 
          axis.text.x = element_text(colour="black", size = 10, margin=unit(c(0.5,0.5,0.1,0.1), "cm"), face = "bold"), # x-axis text in fontsize 10
          axis.text.y = element_text(colour="black", size = 10, margin=unit(c(0.5,0.5,0.1,0.1), "cm"), face = "bold"), # y-axis text in fontsize 10
          legend.text = element_text(size = 10, face = "bold"),
          legend.title = element_text(colour="black", size=12, face="bold"),
          axis.ticks.length = unit(-0.25, "cm")) # ticks facing inward with 0.25cm length    
  
  p <- ggplotly(p, tooltip = c("variance")) %>%
    plotly::config(displaylogo = FALSE,
                   modeBarButtons = list(list("zoomIn2d"), 
                                         list("zoomOut2d"), 
                                         list('toImage')), 
                   mathjax = 'cdn')
  
  return(p)    
}