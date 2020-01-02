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
plot_pca <- function(pca_loadings_var_list, metadata, condition, pc_x, pc_y) {
  message("Plot PCA Started...")
  require(plotly)

  cohort_vector <- metadata[, condition]

  xaxis = paste("PC", pc_x, " (", round(pca_loadings_var_list$var_exp[pc_x], 2), "%)", sep = "")
  yaxis = paste("PC", pc_y, " (", round(pca_loadings_var_list$var_exp[pc_y], 2), "%)", sep = "")

  p <- qplot(x = pca_loadings_var_list$loadings[, pc_x], y = pca_loadings_var_list$loadings[, pc_y], color = cohort_vector, main = "PCA Plot", xlab = xaxis, ylab = yaxis) +
  stat_ellipse(geom = "polygon", alpha = 1/2, aes(fill = cohort_vector))

  p <- ggplotly(p)

  message("Plot PCA Completed...")

  return(p)
}
