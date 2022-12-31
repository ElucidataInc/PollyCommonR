#' plot_opls_distance
#' 
#' Plot the distance outliers plot.
#' 
#' @param dist_data A dataframe containing the distance data.
#' @param metadata Dataframe containing samples to cohort mapping information
#' @param condition A metadata column where cohorts are present
#' @param significance A boolean indicating whether to plot significance level lines.
#' @param significance_data A numeric of length 2 containing the significance level distances of both axes.
#' @param interactive make plot interactive (default is TRUE)
#' @param title_label Title of the plot
#' @param x_title The title for x-axis
#' @param y_title The title for y-axis
#' @param marker_size The size of marker point
#' @param title_label_size Size of title label
#' @param axis_title_size Size of axis title
#' @param opls_cohort_text_format set text format of cohort legends
#' @param opls_cohort_text_align align cohort legends
#' @param opls_cohort_title_size set font size of cohort title
#' @param opls_cohort_sample_size set font size of cohorts
#' @param opls_plot_axis_format set axis format
#' @param opls_plot_axis_text_size set axis text size
#' @returns a plotly object
#' @examples 
#' plot_opls_distance(dist_data,Metadata,"condition",significance=TRUE, significance_data)
#' @import ggplot2 plotly
#' @export
plot_opls_distance = function(dist_data, metadata=NULL, condition = NULL,
                              significance = FALSE, significance_data = NULL,
                              interactive = TRUE,
                              title_label = "", x_title = "Score Component Distance", y_title = "Orthogonal Component Distance",
                              marker_size = 6,title_label_size = 18 ,
                              axis_title_size = 14,opls_cohort_text_format= 'bold' ,
                              opls_cohort_text_align= "right" ,opls_cohort_title_size= 18 ,
                              opls_cohort_sample_size= 15 ,opls_plot_axis_format= 'bold' ,
                              opls_plot_axis_text_size= 14){
  require(ggplot2)
  require(plotly)
  
  if (identical(dist_data, NULL)){
    message("dist_data is NULL")
    
    return (NULL)
  }
  
  if(is.null(metadata)){
    stop("Metadata object cannot be NULL.")
  }
  
  if(is.null(condition)){
    stop("Condition must be a column in metadata and cannot be NULL.")
  }
  if (!(condition %in% colnames(metadata))) {
    warning(c(condition, " is not a valid cohort column, please choose from metadata colnames"))
    return(NULL)
  }
  if (length(unique(metadata[condition])) < 1) {
    warning("The number of cohorts should be greater than or equal to 1")
    return(NULL)
  }
  
  #if(length(dist_data)>2){dist_data = dist_data[,1:2]}
  
  names(dist_data) <- c("sdsVn","dsVn")
  
    p <- ggplot2::ggplot(dist_data, aes(x = sdsVn,
                                      y = dsVn,
                                      text = rownames(dist_data),
                                      group = metadata[, condition],
                                      fill = metadata[, condition]
    )) + # calls the ggplot function with dose on the x-axis and len on the y-axis
    geom_point(shape = 21, size = marker_size, alpha = 0.7)  # scatter plot function with shape of points defined as 21 scale.
  
    if(significance==TRUE){
      if(class(significance_data)=='numeric' && length(significance_data)==2){
    
      p  <- p + geom_vline(xintercept = significance_data[1], linetype = "dashed", alpha = 0.7) + 
        geom_hline(yintercept=significance_data[2],linetype = "dashed", alpha = 0.7) #Adds the significance lines to the plot
      }
      else{ warning("Significance data must be a numeric vector of length 2.")}
    
    }
  
  
  p <- p + ggtitle(title_label)+
    labs(x = x_title, y = y_title, fill = condition) + # x and y axis labels
    theme(legend.position = opls_cohort_text_align, legend.direction = "vertical", # legend positioned at the bottom, horizantal direction,
          axis.line = element_line(size = 1, colour = "black"), # axis line of size 1 inch in black color
          panel.grid.major = element_blank(), # major grids included
          panel.grid.minor = element_blank(), # no minor grids
          panel.border = element_blank(), panel.background = element_blank(), # no borders and background color
          plot.title = element_text(colour="black", size = title_label_size, face = "plain", hjust=0.5),
          axis.title = element_text(colour="black", size = axis_title_size, face = opls_plot_axis_format), # axis title
          axis.text.x = element_text(colour="black", size = opls_plot_axis_text_size, margin=unit(c(0.2,0.2,0.1,0.1), "cm"), face = opls_plot_axis_format), # x-axis text in fontsize 10
          axis.text.y = element_text(colour="black", size = opls_plot_axis_text_size, margin=unit(c(0.2,0.2,0.1,0.1), "cm"), face = opls_plot_axis_format), # y-axis text in fontsize 10
          legend.text = element_text(size = opls_cohort_sample_size, face = opls_cohort_text_format),
          legend.title = element_text(colour="black", size= opls_cohort_title_size, face= opls_cohort_text_format),
          axis.ticks.length = unit(0.25, "cm"))
  
  if(interactive==TRUE){
  
  p <- p + theme(legend.title = element_blank())
  p <- plotly::ggplotly(p, tooltip = "text") %>% layout(hovermode = TRUE) %>%
    add_annotations(text=condition, xref="paper", yref="paper",
                    x=1.03, xanchor="left",
                    y=0.97, yanchor="bottom",
                    font = list(size = (23.91034/18)*opls_cohort_title_size),
                    legendtitle=TRUE, showarrow=FALSE ) %>% plotly::config(displaylogo = FALSE,
                                                                           modeBarButtons = list(list("zoomIn2d"), 
                                                                                                 list("zoomOut2d"), 
                                                                                                 list('toImage')), 
                                                                           mathjax = 'cdn')
  }
  
  return (p)
  
  
}