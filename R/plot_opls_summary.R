#' plot_opls_summary
#' 
#' Plots the summary statistics of the opls object.
#' 
#' @param opl An opls object
#' @param title_label Title of the plot
#' @param title_label_size Size of title label
#' @param axis_title_size Size of axis title
#' @param opls_cohort_text_format set text format of legends
#' @param opls_legend_text_align align legends
#' @param opls_legend_title_size set font size of cohort title
#' @param opls_plot_axis_format set axis format
#' @param opls_plot_axis_text_size set axis text size
#' @returns a ggplot2 object
#' @examples 
#' plot_opls_summary(OPLS_object,title_label="My Graph")
#' @import ggplot2 ropls reshape2
#' @export
plot_opls_summary <- function(opl,
                              title_label = "Model Summary Diagnostics",title_label_size = 18 ,
                              axis_title_size = 14,opls_cohort_text_format= 'bold' ,
                              opls_legend_text_align= "right" ,opls_legend_title_size= 18 ,
                              opls_plot_axis_format= 'bold' ,
                              opls_plot_axis_text_size= 14){
  
  require(ggplot2)
  require(ropls)
  require(reshape2)
  
  if(class(opl)!="opls"){
    stop("Object not of 'opls' class.")
  }
  if(opl@typeC!="OPLS-DA"){
    stop("You are not using a valid 'OPLS-DA' computed on the 'opls' object. Check the 'typeC' subclass in the opls s4 object.")
  }
  if(is.null(opl)){
    stop("Scores cannot be computed on NULL object.")
  }
  
  ### Data pre-processing before plotting.
  summ = opl@modelDF
  summ = summ[!(rownames(summ)=="sum"),!(names(summ)=="Signif.")]
  rownames(summ) = c("Score Component",paste("orthoComponent",1:(nrow(summ)-1),sep='-'))
  summ[nrow(summ)+1,] = names(summ)
  rownames(summ) = c(rownames(summ)[1:(nrow(summ))-1],"stat") # Adding a stat column to be used as x axis while plotting
  
  suppressWarnings({
    e = melt(as.data.frame(t(summ)),
             id.vars = "stat")
  })
  e$value = as.numeric(e$value) # melt converted the statitics values to strings
  #We have removed the signif. col and summ row before melt on summ
  
  p <- ggplot(data = e,mapping = aes(x = stat, y = value, fill = variable)) + 
    geom_col(position = position_dodge())+
    ggtitle(title_label)+
    
    labs(x = "Statistics", y = "Value", fill = "Component") + # x and y axis labels
    theme(legend.position = opls_legend_text_align, legend.direction = "vertical", # legend positioned at the bottom, horizantal direction,
          axis.line = element_line(size = 1, colour = "black"), # axis line of size 1 inch in black color
          plot.title = element_text(colour="black", size = title_label_size, face = "plain", hjust=0.5),
          axis.title = element_text(colour="black", size = axis_title_size, face = opls_plot_axis_format), # axis title
          panel.grid.major = element_blank(), # major grids included
          panel.grid.minor = element_blank(), # no minor grids
          panel.border = element_blank(), panel.background = element_blank(), # no borders and background color
          legend.title = element_text(colour="black", size= opls_legend_title_size, face= opls_cohort_text_format),
          axis.ticks.length = unit(0.25, "cm"))
  p <- plotly::ggplotly(p) %>%
    plotly::config(displaylogo = FALSE, 
                     modeBarButtons = list(list("zoomIn2d"), list("zoomOut2d"), 
                                           list("toImage")), mathjax = "cdn")
  p <- update_plotly_config(p)
  
  return(p)
  
}