#' plot_pca
#'
#' makes a plotly PCA plot
#'
#' @param PCAObj_Summary A list of summary of prcomp function.
#' @param metadata dataframe containing samples to cohort mapping information
#' @param condition A metadata column where cohorts are present
#' @param pc_x PC to keep on x-axis
#' @param pc_y PC to keep on y-axis
#' @param show_ellipse show ellipse on plot (default is FALSE)
#' @param interactive make plot interactive (default is TRUE)
#' @param color_palette The vector of colors
#' @param title_label Title of the plot
#' @param marker_size The size of marker point
#' @param title_label_size Size of title label
#' @param axis_title_size Size of axis title
#' @param pca_cohort_text_format set text format of cohort legends
#' @param pca_cohort_text_align align cohort legends
#' @param pca_cohort_title_size set font size of cohort title
#' @param pca_cohort_sample_size set font size of cohorts
#' @param pca_plot_axis_format set axis format
#' @param pca_plot_axis_text_size set axis text size
#' @return plotly object
#' @examples
#' plot_pca(PCAObj_Summary, metadata, 'Cohort', pc_x = 1, pc_y = 2, interactive = TRUE, pca_cohort_text_format= 'bold', pca_cohort_text_align= "right",
#'          pca_cohort_title_size= 18, pca_cohort_sample_size= 15, pca_plot_axis_format= 'bold', pca_plot_axis_text_size= 14)
#' @import stringr ggplot2 plotly
#' @export
plot_pca <- function(PCAObj_Summary, metadata, condition, pc_x = 1, pc_y = 2, 
                     show_ellipse = FALSE, interactive = TRUE, color_palette = NULL,
                     title_label = "", marker_size = 6, title_label_size = 18, 
                     axis_title_size = 14, pca_cohort_text_format= 'bold', 
                     pca_cohort_text_align= "right", pca_cohort_title_size= 18, 
                     pca_cohort_sample_size= 15, pca_plot_axis_format= 'bold', 
                     pca_plot_axis_text_size= 14) {
  message("Plot PCA Started...")
  require(stringr)
  require(ggplot2)
  require(plotly)

  if (identical(PCAObj_Summary, NULL)){
    message("PCAObj_Summary is NULL")
    
    return (NULL)
  }
  
  if (!(condition %in% colnames(metadata))) {
    warning(c(condition, " is not a valid cohort column, please choose from metadata colnames"))
    return(NULL)
  }
  
  if (length(unique(metadata[[condition]])) < 1) {
    warning("The number of cohorts should be greater than or equal to 1")
    return(NULL)
  }
  
  pca_var_df <- as.data.frame(PCAObj_Summary$x)
  pca_var_df$variable <- rownames(pca_var_df)
  metadata[, condition] <- gsub(",", "_", as.character(metadata[, condition]), fixed = TRUE)
  comb_pca_metadata <- merge(pca_var_df, metadata, by.x = 'variable', by.y = 1, sort = FALSE)
  #print (head(comb_pca_metadata)) 
  p <- ggplot2::ggplot(comb_pca_metadata, aes(x = eval(parse(text=paste("PC", pc_x, sep = ""))),
                                              y = eval(parse(text=paste("PC", pc_y, sep = ""))),
                                              text = paste(variable, !! sym(condition), sep = "<br>"),
                                              group = comb_pca_metadata[, condition],
                                              fill = comb_pca_metadata[, condition]
  )) + # calls the ggplot function with dose on the x-axis and len on the y-axis
    geom_point(shape = 21, size = marker_size, alpha = 0.7) + # scatter plot function with shape of points defined as 21 scale.
    ggtitle(title_label) +
    labs(x = paste("PC",pc_x, '(', round(PCAObj_Summary$importance[2,pc_x]*100, 2), '%)'),
         y = paste("PC",pc_y, '(', round(PCAObj_Summary$importance[2,pc_y]*100, 2), '%)'),
         fill = condition) + # x and y axis labels
    theme(legend.position = pca_cohort_text_align, legend.direction = "vertical", # legend positioned at the bottom, horizantal direction,
          axis.line = element_line(size = 1, colour = "black"), # axis line of size 1 inch in black color
          panel.grid.major = element_blank(), # major grids included
          panel.grid.minor = element_blank(), # no minor grids
          panel.border = element_blank(), panel.background = element_blank(), # no borders and background color
          plot.title = element_text(colour="black", size = title_label_size, face = "plain", hjust=0.5),
          axis.title = element_text(colour="black", size = axis_title_size, face = pca_plot_axis_format), # axis title
          axis.text.x = element_text(colour="black", size = pca_plot_axis_text_size, margin=unit(c(0.2,0.2,0.1,0.1), "cm"), face = pca_plot_axis_format), # x-axis text in fontsize 10
          axis.text.y = element_text(colour="black", size = pca_plot_axis_text_size, margin=unit(c(0.2,0.2,0.1,0.1), "cm"), face = pca_plot_axis_format), # y-axis text in fontsize 10
          legend.text = element_text(size = pca_cohort_sample_size, face = pca_cohort_text_format),
          legend.title = element_text(colour="black", size= pca_cohort_title_size, face= pca_cohort_text_format),
          axis.ticks.length = unit(0.25, "cm"))
  if (show_ellipse == TRUE){
    p <- p + stat_ellipse(geom = "polygon", alpha = 1/6, aes(fill = comb_pca_metadata[,condition]))
  }
  
  if (!identical(color_palette, NULL)){
    if (length(unique(metadata[, condition])) != length(color_palette)){
      warning(c("The number of colors in color palette should be equal to number of cohort conditions which is ", length(unique(metadata[, condition]))))
      return(NULL)
    }
    
    color_bool <- function(pallete){ sapply(pallete, function(x) { tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE)})}  
    if (!all(color_bool(color_palette))){
      warning("The color_palette is not a valid color vector")
      return(NULL)  
    }
    
    try(p <- p + ggplot2::scale_fill_manual(values = color_palette), silent = TRUE)
  }
  
  if (interactive == TRUE){
    p <- p + theme(legend.title = element_blank())
    p <- plotly::ggplotly(p, tooltip = "text") %>% layout(hovermode = TRUE) %>%
      add_annotations(text=condition, xref="paper", yref="paper",
                      x=1.03, xanchor="left",
                      y=0.97, yanchor="bottom",
                      font = list(size = (23.91034/18)*pca_cohort_title_size),
                      legendtitle=TRUE, showarrow=FALSE ) %>% plotly::config(displaylogo = FALSE,
                                                                             modeBarButtons = list(list("zoomIn2d"), 
                                                                                                   list("zoomOut2d"), 
                                                                                                   list('toImage')), 
                                                                             mathjax = 'cdn')
    
    cohorts_vec <- unique(metadata[,condition])
    chr_size_ratio <- max(sapply(cohorts_vec, function(x) nchar(x) / nchar(condition)))
    for (cohort_index in 1:length(p$x$data)){
      name <- p$x$data[[cohort_index]]$name
      cohort <- name  
      if (!identical(name, "1")){
        re <- "\\(([^()]+)\\)"
        str_extracted <- stringr::str_extract_all(name, re)[[1]]
        if (length(str_extracted) > 0){
          cohort <- strsplit(gsub(re, "\\1", str_extracted), ",")[[1]][1] 
        }  
      } 
      if (chr_size_ratio < 1){
        cohort <- paste0(cohort, strrep(" ", (nchar(condition) + round(1/chr_size_ratio, 0)/1.5)))
      } 
      p$x$data[[cohort_index]]$name <- cohort 
    }
  }
  
  message("Plot PCA Completed...")
  
  return(p)
}