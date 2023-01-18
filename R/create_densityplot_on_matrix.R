#' create_densityplot_on_matrix
#'
#' Makes densityplot on matrix
#'
#' @param sample_raw_mat sample_raw_mat matrix/dataframe containing raw values
#' @param data_type use datafrom whole matrix or a column ("all" or "column")
#' @param col_name A column name from sample_raw_mat
#' @param x_label Label x-axis
#' @param y_label Label y-axis
#' @param title_label Title of the plot
#' @param interactive Make plot interactive using plotly
#' @return ggplot object
#' @examples
#' create_densityplot_on_matrix(sample_raw_mat = NULL, data_type = "all", col_name = NULL)
#' @import ggplot2 plotly
#' @export
create_densityplot_on_matrix <- function(sample_raw_mat = NULL, data_type = "all", col_name = NULL,
                                         x_label = NULL, y_label = NULL, title_label = NULL, 
                                         interactive = FALSE){
  require(ggplot2)
  require(plotly)
  message("Create Densityplot On Matrix Started...")
  
  if (!any(data_type %in% c("column", "all"))){
    warning("Please input valide data_type ('column' or 'all')")
    return (NULL)
  }
  
  if (data_type == "column"){
    if (identical(col_name, NULL)==TRUE){
      warning("The col_name parameter is missing, please provide valid column name from sample_raw_mat")
      return (NULL)
    }
    if (!(col_name %in% colnames(sample_raw_mat))){
      warning(c(col_name, " is not present in sample_raw_mat columns"))
      return (NULL)
    }
    mat_df <- sample_raw_mat
    col_value <- col_name
  }
  if (data_type == "all"){
    mat_df <- data.frame(value = as.vector(t(as.matrix(sample_raw_mat))))
    col_value <- "value"
  }
  
  if (identical(y_label, NULL)){
    y_label <- "Density"
  }  
  
  p <- ggplot(mat_df, aes_string(x=col_value))+
    geom_density()+
    labs(title = title_label, x = x_label, y = y_label)+
    theme(axis.line = element_line(size = 1, colour = "black"), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_blank(), panel.background = element_blank(), 
          plot.title = element_text(colour = "black", size = 18, 
                                    face = "plain", hjust=0.5), 
          axis.title = element_text(colour = "black",
                                    size = 14, face = "plain"), 
          axis.text.x = element_text(colour = "black", size = 10, angle = 90, 
                                     margin = unit(c(0.2, 0.2, 0.1, 0.1), "cm"), face = "plain"),
          axis.text.y = element_text(colour = "black", size = 10, 
                                     margin = unit(c(0.2, 0.2, 0.1, 0.1), "cm"), face = "plain"), 
          axis.ticks.length = unit(0.25, "cm"))
  
  if (interactive == TRUE) {
    p <- plotly::ggplotly(p) %>% 
      plotly::config(displaylogo = FALSE, 
                     modeBarButtons = list(list("zoomIn2d"), list("zoomOut2d"), 
                                           list("toImage")), mathjax = "cdn")
  }
  
  message("Create Densityplot On Matrix Completed...")
  
  return (p)
}