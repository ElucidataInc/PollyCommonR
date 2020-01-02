#' compute_pca
#'
#' performs pca on the sample raw matrix
#'
#' @param sample_raw_mat matrix/dataframe containing raw values.
#'
#' @return A list of loadings and pc values
#' @export
compute_pca <- function(sample_raw_mat) {
  message("Compute PCA Started...")
  require(stats)
  require(matrixStats)

  variancerow <- matrixStats::rowVars(as.matrix(sample_raw_mat))
  sample_raw_mat <- sample_raw_mat[!(variancerow <= 0.01),]
  sample_raw_mat <- sample_raw_mat[!apply(sample_raw_mat, 1, anyNA), ]
  data_mat <- as.data.frame(t(sample_raw_mat))
  pca_result <- stats::prcomp(data_mat, scale = T)

  loadings <- data.frame(pca_result$x,
    cohort = row.names(pca_result$x)
  )
  eigs <- pca_result$sdev^2
  var_exp <- eigs / sum(eigs)
  var_exp <- var_exp * 100

  message("Compute PCA Completed...")
  return(list(loadings = loadings, var_exp = var_exp))
}
