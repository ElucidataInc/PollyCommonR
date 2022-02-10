#' plot_pca
#'
#' makes a PCA plot in both interactive (plotly) and non interactive (ggplot) modes
#'
#' @param PCAObj_Summary A list of summary of prcomp function.
#' @param metadata dataframe containing samples to cohort mapping information
#' @param condition A metadata column where cohorts are present
#' @param pc_x PC to keep on x-axis
#' @param pc_y PC to keep on y-axis
#' @param show_ellipse show ellipse on plot (default is FALSE)
#' @param annotate_id A vector of ids (rownames of pca summary data) to be annotated
#' @param annotate_text_size The size of annotated text
#' @param interactive make plot interactive (default is TRUE)
#' @param color_palette The named vector with colors as values and cohorts as keys
#' @param group_shape The named vector with shapes as values (numeric values between [0, 25]) and cohorts as keys
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
#' @return plotly or ggplot object
#' @examples
#' plot_pca(PCAObj_Summary, metadata, 'Cohort', pc_x = 1, pc_y = 2, interactive = TRUE, pca_cohort_text_format= 'bold', pca_cohort_text_align= "right",
#'          pca_cohort_title_size= 18, pca_cohort_sample_size= 15, pca_plot_axis_format= 'bold', pca_plot_axis_text_size= 14)
#' @import stringr ggplot2 plotly
#' @export
plot_pca <- function(PCAObj_Summary, metadata, condition, pc_x = 1, pc_y = 2, 
                     show_ellipse = FALSE, annotate_id = NULL, annotate_text_size = 3,
                     interactive = TRUE, color_palette = NULL, group_shape = NULL, title_label = "",
                     marker_size = 6, title_label_size = 18, axis_title_size = 14, 
                     pca_cohort_text_format= 'bold', pca_cohort_text_align= "right", 
                     pca_cohort_title_size= 18, pca_cohort_sample_size= 15, 
                     pca_plot_axis_format= 'bold', pca_plot_axis_text_size= 14) {
  message("Plot PCA Started...")
  require(stringr)
  require(ggplot2)
  require(plotly)
  require(ggrepel)
  
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
  pca_var_df$id <- row.names(pca_var_df)
  metadata[, condition] <- gsub(",", "_", as.character(metadata[, condition]), fixed = TRUE)
  comb_pca_metadata <- merge(pca_var_df, metadata, by.x = 'id', by.y = 1, sort = FALSE)
  
  if (!identical(annotate_id, NULL)){
    common_annotate_id <- base::intersect(annotate_id, comb_pca_metadata$id)
    if (length(common_annotate_id) < 1){
      warning("The annotate ids are not matched with rownames of pca summary data")    
    }
    else {
      diff_annotate_id <- base::setdiff(annotate_id, comb_pca_metadata$id)
      if (length(diff_annotate_id) >= 1){
        warning(paste("The following annotate ids are not matched with rownames of pca summary data :", paste(sQuote(diff_annotate_id), collapse = ", "), "\n", collapse = " ")) 
      }  
    }       
  }    
  
  if (!identical(group_shape, NULL)) {
    p <- ggplot2::ggplot(comb_pca_metadata, aes(x = eval(parse(text=paste("PC", pc_x, sep = ""))),
                                                y = eval(parse(text=paste("PC", pc_y, sep = ""))),
                                                text = paste(id, !! sym(condition), sep = "<br>"),
                                                group = !! sym(condition),
                                                shape = !! sym(condition),
                                                fill = !! sym(condition),
                                                colour = !! sym(condition))) +
      geom_point(size = marker_size, alpha = 0.7)
  } else {
    p <- ggplot2::ggplot(comb_pca_metadata, aes(x = eval(parse(text=paste("PC", pc_x, sep = ""))),
                                                y = eval(parse(text=paste("PC", pc_y, sep = ""))),
                                                text = paste(id, !! sym(condition), sep = "<br>"),
                                                group = !! sym(condition),
                                                fill = !! sym(condition),
                                                colour = !! sym(condition))) +
      geom_point(shape = 21, size = marker_size, alpha = 0.7)
  }
  
  p <- p +
    ggtitle(title_label) +
    labs(x = paste("PC",pc_x, '(', round(PCAObj_Summary$importance[2,pc_x]*100, 2), '%)'),
         y = paste("PC",pc_y, '(', round(PCAObj_Summary$importance[2,pc_y]*100, 2), '%)'),
         colour = condition,
         shape = condition,
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
    p <- p + ggplot2::stat_ellipse(geom = "polygon", alpha = 1/6, aes(fill = !! sym(condition)), size = 0)
  }
  
  if (!identical(color_palette, NULL)){
    if (length(unique(metadata[, condition])) != length(color_palette)){
      warning(c("The number of colors in color palette should be equal to the number of cohort conditions which is ", length(unique(metadata[, condition]))))
      return(NULL)
    }
    
    if (!identical(names(color_palette), NULL)){
      names(color_palette) <- gsub(",", "_", names(color_palette), fixed = TRUE)
      diff_cohort <- base::setdiff(unique(metadata[, condition]), names(color_palette))
      if (length(diff_cohort) > 0){
        warning(paste0("The following cohorts are absent from color_palette names: ", paste0(diff_cohort, collapse = ", ")))
        return(NULL) 
      }     
    }
    
    color_bool <- function(pallete){ sapply(pallete, function(x) { tryCatch(is.matrix(grDevices::col2rgb(x)), error = function(e) FALSE)})}  
    if (!all(color_bool(color_palette))){
      warning("The color_palette is not a valid color vector")
      return(NULL)  
    }
    
    try(p <- p + ggplot2::scale_colour_manual(values = color_palette), silent = TRUE)
    try(p <- p + ggplot2::scale_fill_manual(values = color_palette), silent = TRUE)                                        
  }
  
  if (!identical(group_shape, NULL)){
    if (length(unique(metadata[, condition])) != length(group_shape)){
      warning(c("The number of shapes in group shape vector should be equal to the number of cohort conditions which is ", length(unique(metadata[, condition]))))
      return(NULL)
    }
    
    if (!identical(names(group_shape), NULL)){
      names(group_shape) <- gsub(",", "_", names(group_shape), fixed = TRUE)
      diff_cohort <- base::setdiff(unique(metadata[, condition]), names(group_shape))
      if (length(diff_cohort) > 0){
        warning(paste0("The following cohorts are absent from group_shape names: ", paste0(diff_cohort, collapse = ", ")))
        return(NULL) 
      }     
    }
    
    shape_bool <- function(g_shape){ tryCatch(suppressWarnings(as.numeric(g_shape)) %in% 0:25, error = function(e) FALSE)}                     
    if (!all(shape_bool(group_shape))){
      warning("The group_shape is not a valid shape vector")
      return(NULL)  
    }
    
    try(p <- p + ggplot2::scale_shape_manual(values = group_shape), silent = TRUE)
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
  
  if (!identical(annotate_id, NULL)) {
    if (identical(interactive, FALSE)) {
      p <- p + ggrepel::geom_text_repel(data = subset(comb_pca_metadata, id %in% annotate_id), aes_string(label = "id"), color = "black", size = annotate_text_size,
                                        box.padding = unit(0.6, 'lines'),
                                        point.padding = unit(0.6, 'lines'),
                                        segment.size = 0.5,
                                        arrow = arrow(length = unit(0.01, 'npc')),
                                        force = 0.5,
                                        max.iter = 3e3)
    } else {
      filtered_data <-   comb_pca_metadata[comb_pca_metadata$id %in% annotate_id, , drop = FALSE]
      annotate_data <- list(x = filtered_data[[paste("PC", pc_x, sep = "")]], y = filtered_data[[paste("PC", pc_y, sep = "")]],
                            text = filtered_data$id, xref = "x", yref = "y", showarrow = TRUE, arrowhead = 3.5, ax = 20, ay = -40,
                            font = list(color = "black", size = 4 * annotate_text_size))                               
      p <- p %>% layout(annotations = annotate_data)   
    }
  }                               
  
  message("Plot PCA Completed...")
  
  return(p)
}