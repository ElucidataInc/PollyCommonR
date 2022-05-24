#' compute_fold_change
#'
#' Calculate fold change on dataframe/matrix where samples are in columns and features are in rows.
#'
#' @param data_mat A dataframe/matrix containing samples in columns and features in rows.
#' @param samples_a A vector of samples in cohort A.
#' @param samples_b A vector of samples in cohort B.
#' @param paired A logical indicating whether to use paired samples for fold change calculation.
#' @param log_flag A logical variable indicating whether data is log transformed or not. 
#' If TRUE then fold change will be calculated as B - A (logFC) otherwise B/A (FC).
#' @return A dataframe with calculated statistics.
#' @examples 
#' compute_fold_change(data_mat, samples_a, samples_b)
#' @export
compute_fold_change <- function(data_mat = NULL, samples_a = NULL, samples_b = NULL,
                                paired = FALSE, log_flag = FALSE){
  message("Compute Fold Change Started...")
  
  if (identical(data_mat, NULL)){
    warning("The data_mat is NULL")
    return (NULL)  
  }
  
  if (!identical(class(data_mat), "data.frame") && !identical(class(data_mat), "matrix")){
    warning("The data_mat is not a dataframe/matrix, please provide valid norm_data")
    return(NULL) 
  }
  
  common_samples_a <- base::intersect(samples_a, colnames(data_mat))
  if (length(common_samples_a) < 1){
    warning("No common samples found between samples_a and data matrix")
    return(NULL)  
  }    
  
  common_samples_b <- base::intersect(samples_b, colnames(data_mat))  
  if (length(common_samples_b) < 1){
    warning("No common samples found between samples_b and data matrix")
    return(NULL)  
  }
  
  samples_diff <- base::setdiff(c(samples_a, samples_b), colnames(data_mat)) 
  if (length(samples_diff) > 0){
    warning(paste0("The following samples are not present in the data matrix :", paste0(samples_diff, collapse = ", ")))
  }
  
  if (identical(paired, TRUE)){
    if (!(length(common_samples_a) == length(common_samples_b))){
      warning("The samples_a and samples_b should have equal number of samples for paired comparison")
      return(NULL)    
    }
  }
  
  fold_change_df <- NULL
  if (identical(paired, TRUE)){
    if (identical(log_flag, FALSE)){
      data_mat <- log2(data_mat)
    }
    
    samples_a_mat <- as.matrix(data_mat[, common_samples_a, drop = FALSE])
    samples_b_mat <- as.matrix(data_mat[, common_samples_b, drop = FALSE])
    logfc <- apply((samples_b_mat - samples_a_mat), 1, mean, na.rm = TRUE)
    
    if (identical(log_flag, FALSE)){
      fold_change_df <- data.frame(id = names(logfc), FC = 2^logfc, logFC = logfc, stringsAsFactors = FALSE, check.names = FALSE)
    }  
    else {
      fold_change_df <- data.frame(id = names(logfc), logFC = logfc, stringsAsFactors = FALSE, check.names = FALSE)
    }      
  }
  else {
    samples_a_mat <- as.matrix(data_mat[, common_samples_a, drop = FALSE])
    samples_b_mat <- as.matrix(data_mat[, common_samples_b, drop = FALSE])
    
    if (identical(log_flag, FALSE)){
      fold_change <- base::rowMeans(samples_b_mat, na.rm = TRUE) / base::rowMeans(samples_a_mat, na.rm = TRUE)
      fold_change_df <- data.frame(id = names(fold_change), FC = fold_change, logFC = log2(fold_change), stringsAsFactors = FALSE, check.names = FALSE)  
    }
    else {
      logfc <- base::rowMeans(samples_b_mat, na.rm = TRUE) - base::rowMeans(samples_a_mat, na.rm = TRUE)
      fold_change_df <- data.frame(id = names(logfc), logFC = logfc, stringsAsFactors = FALSE, check.names = FALSE)      
    }
  }
  
  message("Compute Fold Change Completed...")
  
  return (fold_change_df)  
}