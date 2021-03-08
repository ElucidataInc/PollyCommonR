#' create_boxplot_on_matrix
#'
#' Makes boxplot on matrix
#'
#' @param sample_raw_mat sample_raw_mat matrix/dataframe containing raw values
#' @param x_label Label x-axis
#' @param y_label Label y-axis
#' @param title_label Title of the plot
#' @param flip_coord Flip Coordinates
#' @param plot_axis_format font face ("plain", "italic", "bold", "bold.italic")
#' @param plot_axis_text_size set axis text size
#' @return ggplot object
#' @examples
#' create_boxplot_on_matrix(sample_raw_mat = NULL, x_label = "Sample", y_label = "Raw Intensity", plot_axis_format = 'plain', plot_axis_text_size = 10)
#' @import ggplot2
#' @export
create_boxplot_on_matrix <- function(sample_raw_mat = NULL, x_label = "", y_label = "",
                                     title_label = "", flip_coord = FALSE, 
                                     plot_axis_format = 'plain', plot_axis_text_size = 10){
  require(ggplot2)
  message("Create Boxplot On Matrix Started...")
  
  p <- ggplot(stack(sample_raw_mat), aes(x = ind, y = values, fill=ind))+
    geom_boxplot(show.legend = FALSE)+
    ggtitle(title_label)+
    labs(x = x_label,
         y = y_label)+
    ggsci::scale_color_aaas() + # filling the point colors
    theme(axis.line = element_line(size = 1, colour = "black"), # axis line of size 1 inch in black color
          panel.grid.major = element_blank(), # major grids included
          panel.grid.minor = element_blank(), # no minor grids
          panel.border = element_blank(), panel.background = element_blank(), # no borders and background color
          plot.title = element_text(colour="black", size = 18, face = "plain", hjust=0.5),
          axis.title = element_text(colour="black", size = 14, face = plot_axis_format), # axis title
          axis.text.x = element_text(colour="black", size = plot_axis_text_size, angle = 90,
                                     hjust = 1, margin=unit(c(0.2,0.2,0.1,0.1), "cm"),
                                     face = plot_axis_format), # x-axis text in fontsize 10
          axis.text.y = element_text(colour="black", size = plot_axis_text_size,
                                     margin=unit(c(0.2,0.2,0.1,0.1), "cm"), 
                                     face = plot_axis_format), # y-axis text in fontsize 10
          axis.ticks.length = unit(0.25, "cm"))
  
  if (flip_coord == TRUE ){
    p <- p +
      coord_flip()+
      theme(axis.text.x = element_text(angle = 0))
  }
  
  message("Create Boxplot On Matrix Completed...")
  
  return(p)
}