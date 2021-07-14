#' opls_s_plot
#' 
#' @param opl The computed "opls" object (ropls library) on which the plot has to be made.
#' @param x_label The label for x axis.
#' @param y_label The label for y axis.
#' @param title_label The title for the plot
#' 
#' @return A plotly plot
#' 
#' @import ropls plotly
#' @export

opls_s_plot =  function(opl, x_label = "Loadings", y_label = "Correlation", title_label = "S-Plot for OPLS-DA"){
  if(class(opl)!="opls"){
    stop("Object not of 'opls' class.")
  }
  if(opl@typeC!="OPLS-DA"){
    warning("You are not using a valid OPLS-DA computed on the 'opls' object. Check the 'typeC' subclass in the opls s4 object.")
  }
  if(is.null(opl)){
    stop("NULL object cannot be plotted.")
  }
  
  require("plotly")
  require("ropls")
  
  plot_data = as.data.frame(cbind(opl@loadingMN,opl@loadingMN/(sd(opl@scoreMN)*opl@xSdVn)))
  names(plot_data) = c("Loading","Correlation")
  
  plot = plotly::plot_ly(data = plot_data, x = ~Loading, y= ~Correlation,
                         type="scatter",
                         mode = "markers",
                         color = ~Correlation,
                         showlegend = FALSE,
                         text = rownames(plot_data),
                         hovertemplate = paste('<b>Gene</b> :<i>%{text}</i>',
                                               '<b><b>Loading</b>: %{x: .4f} </b>',
                                               '<b><b>Correlation</b>: %{y: .4f} </b>',
                                               sep="\n",
                                               "<extra></extra>")
  ) 
  plot = plot %>% layout(title = list(text = title_label, xref = "paper", yref = "paper"),
                         yaxis = list(title = y_label),
                         xaxis = list(title = x_label))
  
  plot = plot %>% plotly::config(displaylogo = FALSE,
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
  
  return(plot)
}