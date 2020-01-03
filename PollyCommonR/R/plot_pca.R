#' plot_pca
#'
#' makes a plotly PCA plot
#'
#' @param pca_loadings_var_list list containing the PC loadings for every sample and the variances explained, returned by the compute_pca function.
#' This returns a plotly object of the PCA plot for all the samples.
#' @param metadata dataframe containing samples to cohort mapping information
#' @param condition A metadata column where cohorts are present
#' @param pc_x PC to keep on x-axis
#' @param pc_y PC to keep on y-axis
#'
#' @return plotly object
#' @export
plot_pca <- function(PCAObj_Summary, metadata, condition, pc_x, pc_y, interactive = TRUE, pca_cohort_text_format= 'bold', pca_cohort_text_align= "right",
                     pca_cohort_title_size= 18, pca_cohort_sample_size= 15, pca_plot_axis_format= 'bold', pca_plot_axis_text_size= 14) {
  message("Plot PCA Started...")
  require(stringr)
  require(plotly)
  require(ggsci)
  
  pca_var_df <- as.data.frame(PCAObj_Summary$x)
  pca_var_df$variable <- rownames(pca_var_df) 
  comb_pca_metadata <- merge(pca_var_df, metadata_df, by.x = 'variable', by.y = 1, sort = FALSE)
  
  p <- ggplot(comb_pca_metadata, aes(x = eval(parse(text=paste("PC", pc_x, sep = ""))),
                                     y = eval(parse(text=paste("PC", pc_y, sep = ""))),
                                     text = variable, group = comb_pca_metadata[,condition],
                                     fill = comb_pca_metadata[,condition]
  )) + # calls the ggplot function with dose on the x-axis and len on the y-axis
    geom_point(shape = 21, size = 6, alpha = 0.7) + # scatter plot function with shape of points defined as 21 scale.
    labs(x = paste("PC",pc_x, '(', round(PCAObj_Summary$importance[2,pc_x]*100, 2), '%)'),
         y = paste("PC",pc_y, '(', round(PCAObj_Summary$importance[2,pc_y]*100, 2), '%)')) + # x and y axis labels
    guides(color=FALSE, fill = guide_legend('Cohort'))+
    scale_color_aaas() + # filling the point colors
    theme(legend.position = pca_cohort_text_align, legend.direction = "vertical", # legend positioned at the bottom, horizantal direction,
          axis.line = element_line(size = 1, colour = "black"), # axis line of size 1 inch in black color
          panel.grid.major = element_blank(), # major grids included
          panel.grid.minor = element_blank(), # no minor grids
          panel.border = element_blank(), panel.background = element_blank(), # no borders and background color
          axis.title = element_text(colour="black", size = 18, face = pca_plot_axis_format), # axis title
          axis.text.x = element_text(colour="black", size = pca_plot_axis_text_size, margin=unit(c(0.5,0.5,0.1,0.1), "cm"), face = pca_plot_axis_format), # x-axis text in fontsize 10
          axis.text.y = element_text(colour="black", size = pca_plot_axis_text_size, margin=unit(c(0.5,0.5,0.1,0.1), "cm"), face = pca_plot_axis_format), # y-axis text in fontsize 10
          legend.text = element_text(size = pca_cohort_sample_size, face = pca_cohort_text_format),
          legend.title = element_text(colour="black", size= pca_cohort_title_size, face= pca_cohort_text_format),
          axis.ticks.length = unit(-0.25, "cm"))
  
  if (interactive == TRUE){
    p <- ggplotly(p, tooltip = "text") %>% layout(hovermode = "closest")
    for (cohort_index in 1:length(p$x$data)){
      name <- p$x$data[[cohort_index]]$name
      cohort <- strsplit(gsub(re, "\\1", stringr::str_extract_all(name, "\\(([^()]+)\\)")[[1]]),",")[[1]][1]
      p$x$data[[cohort_index]]$name <- cohort
    }
  }
  message("Plot PCA Completed...")
  
  return(p)
}

p <- plot_pca(pca_cal, metadata_df, 'Cohort', 1,2, interactive = TRUE)

