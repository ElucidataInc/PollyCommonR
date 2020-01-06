#' compute_pca
#'
#' performs pca on the sample raw matrix
#'
#' @param sample_raw_mat matrix/dataframe containing raw values.
#' @param num Number of genes/metabolites to be used for pca.
#' @return A list of loadings and pc values
#' @examples
#' compute_pca(sample_raw_mat, num = 1000)
#' @export
compute_pca <- function(sample_raw_mat, num = NULL) {
  message("Compute PCA Started...")
  require(stats)
  require(matrixStats)
  
  variancerow <- matrixStats::rowVars(as.matrix(sample_raw_mat))
  sample_raw_mat <- sample_raw_mat[!(variancerow == 0.01),]
  sample_raw_mat <- sample_raw_mat[!apply(sample_raw_mat, 1, anyNA), ]
  
  if (!identical(num, NULL)){
    temp_df <- apply(sample_raw_mat, 1, function(x) mad(x))
    sample_raw_mat <- as.data.frame(cbind(sample_raw_mat,temp_df))
    sample_raw_mat <- sample_raw_mat[order(-sample_raw_mat$temp_df),]
    data_mat <- as.data.frame(t(head(sample_raw_mat[, !(colnames(sample_raw_mat) %in% c("temp_df"))], num)))
  }else{
    data_mat <- as.data.frame(t(sample_raw_mat))
  }
  
  pca_result <- stats::prcomp(data_mat, scale = T)
  PCAObj_Summary <- summary(pca_result)
  
  message("Compute PCA Completed...")
  return (PCAObj_Summary)
}
