#' create_densityplot_on_matrix
#'
#' Makes densityplot on matrix
#'
#' @param sample_raw_mat sample_raw_mat matrix/dataframe containing raw values
#' @param data_type use datafrom whole matrix or a column ("all" or "column")
#' @param col_name A column name from sample_raw_mat
#' @return ggplot object
#' @examples
#' create_densityplot_on_matrix(sample_raw_mat = NULL, data_type = "all", col_name = NULL)
#' @import ggplot2
#' @export
create_densityplot_on_matrix <- function(sample_raw_mat = NULL, data_type = "all", col_name = NULL){
  require(ggplot2)
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
  
  p <- ggplot(mat_df, aes_string(x=col_value))+
    geom_density()+
    labs(x = NULL, y = "Density")+
    theme(axis.line = element_line(size = 1, colour = "black"), # axis line of size 1 inch in black color
          panel.grid.major = element_blank(), # major grids included
          panel.grid.minor = element_blank(), # no minor grids
          panel.border = element_blank(), panel.background = element_blank(), # no borders and background color
          axis.title = element_text(colour="black", size = 18, face = "plain"), # axis title
          axis.text.x = element_text(colour="black", size = 18,# angle = 90,
                                     margin=unit(c(0.2,0.2,0.1,0.1), "cm"),
                                     face = "plain"), # x-axis text in fontsize 10
          axis.text.y = element_text(colour="black", size = 18,
                                     margin=unit(c(0.2,0.2,0.1,0.1), "cm"), 
                                     face = "plain"), # y-axis text in fontsize 10
          axis.ticks.length = unit(0.25, "cm"))
  
  message("Create Densityplot On Matrix Completed...")
  
  return (p)
}