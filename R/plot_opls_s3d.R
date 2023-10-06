#' plot_opls_s3d
#' 
#' Plots the S-plot for OPLS-DA loadings.
#' 
#' @param s_data A dataframe of size (n x 2) with the two columns corresponding to the Loadings and Correlation. 
#' @param x_label The label for x axis.
#' @param y_label The label for y axis.
#' @param title_label The title for the plot
#' @return A plotly plot
#' @import ropls plotly
#' @export
plot_opls_s3d =  function(s_data=NULL, x_label = "Loadings", y_label = "Correlation", title_label = "S-Plot for OPLS-DA"){
  
  if (identical(s_data, NULL)){
    message("s_data is NULL")
    
    return (NULL)
  }
  if (is.null(as.data.frame(s_data))){
    warning("s_data could not be converted to a data frame. Use a valid input.")
    return(NULL)
  }
  s_data = as.data.frame(s_data)
  
  if(ncol(s_data)>2){
    warning("More than 2 columns detected. Using the first two for plotting.")
    s_data = s_data[,1:2]
  }
  
  if(ncol(s_data)<2){
    warning("Cannot plot on a single column.")
    return(NULL)
  }
  if (length(unique(rownames(s_data))) != nrow(s_data)){
    warning("rownames must be unique. Using number for labels.")
    rownames(s_data) = paste("Row ", 1:nrow(s_data))
  }
  names(s_data) = c("Loading","Correlation")
  plot = plotly::plot_ly(data = s_data, x = ~Loading, y= ~Correlation,
                         type="scatter",
                         mode = "markers",
                         color = ~Correlation,
                         showlegend = FALSE,
                         text = rownames(s_data),
                         hovertemplate = paste('<b>ID</b> :<i>%{text}</i>',
                                               '<b><b>Loading</b>: %{x: .4f} </b>',
                                               '<b><b>Correlation</b>: %{y: .4f} </b>',
                                               sep="\n",
                                               "<extra></extra>")
  ) 
  plot <- plot %>% layout(title = list(text = title_label, xref = "paper", yref = "paper"),
                         yaxis = list(title = y_label),
                         xaxis = list(title = x_label))
  
  plot <- plot %>% plotly::config(displaylogo = FALSE,
                                 modeBarButtons = list(list("zoom2d"),
                                                       list("select2d"),
                                                       list("lasso2d"),
                                                       list("autoScale2d"),
                                                       list("resetScale2d"),
                                                       list("pan2d"),
                                                       list("zoomIn2d"), 
                                                       list("zoomOut2d"),
                                                       list("hoverClosestCartesian"),
                                                       list('toImage')),
                                 mathjax = 'cdn')
  plot <- update_plotly_config(plot)
  
  return(plot)
}