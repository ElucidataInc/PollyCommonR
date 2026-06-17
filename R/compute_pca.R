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
#' @import stats
#' @export
compute_pca <- function(sample_raw_mat = NULL, num = NULL, center = TRUE, scale = TRUE) {
  message("Compute PCA Started...")
  require(stats)

  mat <- data.matrix(sample_raw_mat)

  variancerow <- matrixStats::rowVars(mat, na.rm = FALSE)
  mat <- mat[!(variancerow == 0), , drop = FALSE]

  na_mask <- matrixStats::rowAnyNAs(mat)
  mat <- mat[!na_mask, , drop = FALSE]

  if (nrow(mat) == 0) {
    warning("Not a valid matrix, have NANs or infs in the matrix")
    return(NULL)
  }

  if (!identical(num, NULL)) {
    mad_vals <- matrixStats::rowMads(mat, na.rm = FALSE)
    top_idx  <- order(mad_vals, decreasing = TRUE)[seq_len(min(num, nrow(mat)))]
    data_mat <- t(mat[top_idx, , drop = FALSE])
  } else {
    data_mat <- t(mat)
  }

  pca_result     <- stats::prcomp(data_mat, center = center, scale. = scale)
  PCAObj_Summary <- summary(pca_result)

  message("Compute PCA Completed...")
  return(PCAObj_Summary)
}


