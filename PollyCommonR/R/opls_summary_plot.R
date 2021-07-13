#' opls_summary_plot
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
#' 
#' @returns a ggplot2 object
#' 
#' @examples 
#' opla_summary_plot(OPLS_object,title_label="My Graph")
#' 
#' @import ggplot2 ropls
#' @export

opls_summary_plot = function(opl,
                             title_label = "Model Summary Diagnostics",title_label_size = 18 ,
                             axis_title_size = 14,opls_cohort_text_format= 'bold' ,
                             opls_legend_text_align= "right" ,opls_legend_title_size= 18 ,
                             opls_plot_axis_format= 'bold' ,
                             opls_plot_axis_text_size= 14){
  
  require(ggplot2)
  require(ropls)
  
  if(class(opl)!="opls"){
    stop("Object not of 'opls' class.")
  }
  if(opl@typeC!="OPLS-DA"){
    stop("You are not using a valid 'OPLS-DA' computed on the 'opls' object. Check the 'typeC' subclass in the opls s4 object.")
  }
  if(is.null(opl)){
    stop("Scores cannot be computed on NULL object.")
  }
  
  summ = as.data.frame(opl@modelDF[-3,][,-7])
  #summ$type = rownames(summ)
  summ = as.data.frame(t(summ))
  summ$stat = rownames(summ)
  summ = as.data.frame(rbind(as.matrix(summ[,c(3,1)]),as.matrix(summ[,c(3,2)])))
  p1 = "Model Score"
  o1 = "Orthogonal Score"
  summ$Component = c(p1,p1,p1,p1,p1,p1,o1,o1,o1,o1,o1,o1)
  #summ
  
  p <- ggplot(summ,aes(x=stat,y=p1,fill=Component))+
    geom_bar(stat="identity", position=position_dodge())+
    ggtitle(title_label)+
    
    labs(x = "Statistics", y = "Value") + # x and y axis labels
    theme(legend.position = opls_legend_text_align, legend.direction = "vertical", # legend positioned at the bottom, horizantal direction,
          axis.line = element_line(size = 1, colour = "black"), # axis line of size 1 inch in black color
          plot.title = element_text(colour="black", size = title_label_size, face = "plain", hjust=0.5),
          axis.title = element_text(colour="black", size = axis_title_size, face = opls_plot_axis_format), # axis title
          panel.grid.major = element_blank(), # major grids included
          panel.grid.minor = element_blank(), # no minor grids
          panel.border = element_blank(), panel.background = element_blank(), # no borders and background color
          legend.title = element_text(colour="black", size= opls_legend_title_size, face= opls_cohort_text_format),
          axis.ticks.length = unit(0.25, "cm"))
  
  return(p)
  
}