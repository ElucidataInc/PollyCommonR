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
#' @return plotly object
#' @examples
#' plot_pca3d(PCAObj_Summary, metadata, 'Cohort', pc_x = 1, pc_y = 2, pc_z= 3)
#' @import plotly
#' @export
plot_pca3d <- function(PCAObj_Summary, metadata, condition, pc_x = 1, pc_y = 2, pc_z= 3) {
  message("Plot PCA3D Started...")
  require(plotly)
  
  if (identical(PCAObj_Summary, NULL)){
    message("PCAObj_Summary is NULL")
    
    return (NULL)
  }
  
  pca_var_df <- as.data.frame(PCAObj_Summary$x)
  pca_var_df$variable <- rownames(pca_var_df)
  metadata[,condition] <- gsub(",", "_", as.character(metadata[,condition]), fixed = TRUE)
  comb_pca_metadata <- merge(pca_var_df, metadata, by.x = 'variable', by.y = 1, sort = FALSE)
  
  p <- plot_ly(comb_pca_metadata, x = ~eval(parse(text=paste("PC", pc_x, sep = ""))),
               y = ~eval(parse(text=paste("PC", pc_y, sep = ""))),
               z = ~eval(parse(text=paste("PC", pc_z, sep = ""))),
               color = ~Cohort) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = paste("PC",pc_x, '(', round(PCAObj_Summary$importance[2,pc_x]*100, 2), '%)')),
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