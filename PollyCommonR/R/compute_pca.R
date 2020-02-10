#' compute_pca
#'
#' performs pca on the sample raw matrix
#'
#' @param sample_raw_mat matrix/dataframe containing raw values.
#' @param num Number of genes/metabolites to be used for pca.
#' @param center a logical value indicating whether the variables should be 
#' shifted to be zero centered. Alternately, a vector of length equal the 
#' number of columns of x can be supplied. The value is passed to scale.
#' @param scale a logical value indicating whether the variables should be 
#' scaled to have unit variance before the analysis takes place. The default
#'  is FALSE for consistency with S, but in general scaling is advisable. 
#'  Alternatively, a vector of length equal the number of columns of x can 
#'  be supplied. The value is passed to scale (base function).
#' @return A list of loadings and pc values
#' @examples
#' compute_pca(sample_raw_mat, num = 1000, center = TRUE, scale = TRUE)
#' @import matrixStats stats
#' @export
compute_pca <- function(sample_raw_mat = NULL, num = NULL, center = TRUE, scale = TRUE) {
  message("Compute PCA Started...")
  require(matrixStats)
  require(stats)
 
  variancerow <- matrixStats::rowVars(as.matrix(sample_raw_mat))
  sample_raw_mat <- sample_raw_mat[!(variancerow == 0),]
  sample_raw_mat <- sample_raw_mat[!apply(sample_raw_mat, 1, anyNA), ]
  if (nrow(sample_raw_mat) == 0){
    warning("Not a valid matrix, have NANs or infs in the matrix")
    return (NULL)
  }
  if (!identical(num, NULL)){
    temp_df <- apply(sample_raw_mat, 1, function(x) mad(x))
    sample_raw_mat <- as.data.frame(cbind(sample_raw_mat,temp_df))
    sample_raw_mat <- sample_raw_mat[order(-sample_raw_mat$temp_df),]
    data_mat <- as.data.frame(t(head(sample_raw_mat[, !(colnames(sample_raw_mat) %in% c("temp_df"))], num)))
  }else{
    data_mat <- as.data.frame(t(sample_raw_mat))
  }
  
  pca_result <- stats::prcomp(data_mat, center = center, scale = scale)
  PCAObj_Summary <- summary(pca_result)
  
  message("Compute PCA Completed...")
  return (PCAObj_Summary)
}
