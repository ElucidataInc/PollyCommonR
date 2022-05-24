#' plot_opls_score
#' 
#' Plots the model scores vs orthogonal score for each sample.
#'  
#' @param scores_data A dataframe containing the scores data.
#' @param metadata Dataframe containing samples to cohort mapping information. The first column must be ssame as the rownames in in the scores_data.
#' @param condition A metadata column where cohorts are present
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
#' @returns A plotly object
#' @examples plot_opls_score(scores_data,metadat)
#' @import ggplot2 plotly
#' @export 
plot_opls_score <- function(scores_data, metadata=NULL, condition = NULL,
                            show_ellipse = FALSE,interactive = TRUE,
                            title_label = "", x_title = "Score", y_title = "Orthogonal Score",
                            marker_size = 6,title_label_size = 18,
                            axis_title_size = 14,opls_cohort_text_format= 'bold' ,
                            opls_cohort_text_align= "right" ,opls_cohort_title_size= 18 ,
                            opls_cohort_sample_size= 15 ,opls_plot_axis_format= 'bold' ,
                            opls_plot_axis_text_size= 14){
  
  
  
  if (identical(scores_data, NULL)){
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
  if(length(scores_data)>2){scores_data = scores_data[,1:2]}
  
  # Data Processing for unequivocal call to columns.
  names(scores_data) = c("p1",paste("o",1:(ncol(scores_data)-1),sep=''))
  
  # Plotting now
  p <- ggplot2::ggplot(scores_data, aes(x = p1, y = o1,
                                        text = rownames(scores_data),#paste(, !! sym(condition), sep = "<br>"),
                                        group = metadata[, condition],
                                        fill = metadata[, condition]
  )) + # calls the ggplot function with dose on the x-axis and len on the y-axis
    geom_point(shape = 21, size = marker_size, alpha = 0.7) + # scatter plot function with shape of points defined as 21 scale.
    ggtitle(title_label)+
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
  if (show_ellipse == TRUE){
    p <- p + stat_ellipse(geom = "polygon", alpha = 1/6, aes(fill = metadata[,condition]))
  }
  
  if(interactive ==TRUE){
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
  return(p)
}
