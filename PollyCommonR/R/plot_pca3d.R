#' plot_pca3d
#'
#' makes a plotly PCA-3D plot
#'
#' @param PCAObj_Summary A list of summary of prcomp function.
#' @param metadata dataframe containing samples to cohort mapping information
#' @param condition A metadata column where cohorts are present
#' @param pc_x PC to keep on x-axis
#' @param pc_y PC to keep on y-axis
#' @param pc_z PC to keep on z-axis
#' @param color_palette The named vector with colors as values and cohorts as keys 
#' @param title_label Title of the plot
#' @param show_legend_title Show legend title (TRUE/FALSE)
#' @return plotly object
#' @examples
#' plot_pca3d(PCAObj_Summary, metadata, 'Cohort', pc_x = 1, pc_y = 2, pc_z= 3)
#' @import plotly
#' @export
plot_pca3d <- function(PCAObj_Summary = NULL, metadata = NULL, condition = NULL, color_palette = NULL,
                       pc_x = 1, pc_y = 2, pc_z= 3, title_label = "", show_legend_title = TRUE) {
  message("Plot PCA3D Started...")
  require(plotly)
  
  if (identical(PCAObj_Summary, NULL)){
    warning("The PCAObj_Summary is NULL")
    return (NULL)
  }
  
  if (identical(metadata, NULL)){
    warning("The metadata is NULL")
    return (NULL)
  }    
  
  if (!(condition %in% colnames(metadata))) {
    warning(c(condition, " is not a valid cohort column, please choose from metadata colnames"))
    return(NULL)
  }
  
  pca_var_df <- as.data.frame(PCAObj_Summary$x)
  pca_var_df$variable <- rownames(pca_var_df)
  metadata[, condition] <- gsub(",", "_", as.character(metadata[, condition]), fixed = TRUE)
  comb_pca_metadata <- merge(pca_var_df, metadata, by.x = 'variable', by.y = 1, sort = FALSE)
  
  if (!identical(color_palette, NULL)){
    if (length(unique(metadata[, condition])) != length(color_palette)){
      warning(c("The number of colors in color palette should be equal to number of cohort conditions which is ", length(unique(metadata[, condition]))))
      return(NULL)
    }
    
    if (!identical(names(color_palette), NULL)){
      names(color_palette) <- gsub(",", "_", names(color_palette), fixed = TRUE)
      diff_cohort <- base::setdiff(unique(metadata[, condition]), names(color_palette))
      if (length(diff_cohort) > 0){
        warning(paste0("The following cohorts are absent from color_palette names: ", paste0(diff_cohort, collapse = ", ")))
        return(NULL) 
      }
      color_palette <- color_palette[unique(comb_pca_metadata[, condition])]                                     
    }
    
    color_bool <- function(pallete){ sapply(pallete, function(x) { tryCatch(is.matrix(grDevices::col2rgb(x)), error = function(e) FALSE)})}  
    if (!all(color_bool(color_palette))){
      warning("The color_palette is not a valid color vector")
      return(NULL)  
    } 
  }    
  
  if (identical(show_legend_title, TRUE)) {
    legend_title = condition
  } else {
    legend_title = NULL
  }                                         
  
  p <- plot_ly(x = comb_pca_metadata[[paste("PC", pc_x, sep = "")]],
               y = comb_pca_metadata[[paste("PC", pc_y, sep = "")]],
               z = comb_pca_metadata[[paste("PC", pc_z, sep = "")]],
               color = comb_pca_metadata[[condition]], colors = color_palette, text = comb_pca_metadata[["variable"]]) %>%                                           
    add_markers() %>%
    layout(title = title_label,
           legend=list(title=list(text = legend_title)),
           scene = list(xaxis = list(title = paste("PC",pc_x, '(', round(PCAObj_Summary$importance[2,pc_x]*100, 2), '%)')),
                        yaxis = list(title = paste("PC",pc_y, '(', round(PCAObj_Summary$importance[2,pc_y]*100, 2), '%)')),
                        zaxis = list(title = paste("PC",pc_z, '(', round(PCAObj_Summary$importance[2,pc_z]*100, 2), '%)'))
           )) %>%
    plotly::config(displaylogo = FALSE,
                   modeBarButtons = list(list("zoomIn2d"),
                                         list("zoomOut2d"),
                                         list('toImage')), 
                   mathjax = 'cdn')
  
  message("Plot PCA3D Completed...")
  
  return(p)
}