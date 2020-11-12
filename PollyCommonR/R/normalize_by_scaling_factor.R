#' normalize_by_scaling_factor
#'
#' normaize data by samplewise scaling factor
#'
#' @param sample_raw_mat A dataframe/matrix containing samples with raw values.
#' @param normalization_agent A dataframe/matrix with samples as rownames and having numeric columns
#' @param scaling_factor_col A numeric column name which is present in normalization_agent dataframe
#' @return A dataframe/marix with normlazed data
#' @examples 
#' normalize_by_scaling_factor(sample_raw_mat, normalization_agent, scaling_factor_col = 1)
#' @export
normalize_by_scaling_factor <- function(sample_raw_mat, normalization_agent, scaling_factor_col = 1){
  message("Normalize By Sample Factor Started...")
  
  sample_raw_mat <- as.data.frame(sample_raw_mat)
  norm_agent_samples_vec <- rownames(normalization_agent)
  
  sample_raw_mat_cols <- colnames(sample_raw_mat)
  common_sample_cols <- base::intersect(sample_raw_mat_cols, norm_agent_samples_vec)
  if (length(common_sample_cols) == 0){
    message("No common samples found, please provide valid data, returning original data")
    return(sample_raw_mat)
  }
  
  if (identical(class(scaling_factor_col), "numeric")){
    if (ncol(normalization_agent) < scaling_factor_col){
      message("scaling_factor_col > number of normalization_agent columns, returning original data")
      return(sample_raw_mat)
    }    
  } else if (identical(class(scaling_factor_col), "character")){
    if (!(scaling_factor_col %in% colnames(normalization_agent))){
      message(c(scaling_factor_col, " column is not present in normalization_agent, returning original data"))
      return(sample_raw_mat)
    }
  } else {
    message("Not a valid scaling_factor_col, returning original data")
    return(sample_raw_mat)
  }
  
  if (!identical(class(normalization_agent[, scaling_factor_col]), "numeric")){
    message("Not not a numeric column present in normalization_agent, returning original data")
    return(sample_raw_mat)        
  }
  
  for (sample in common_sample_cols) {
    if (normalization_agent[sample, scaling_factor_col] != 0) {
      sample_raw_mat[, sample] <- sample_raw_mat[, sample] / normalization_agent[sample, scaling_factor_col]
    }
  }
  
  message("Normalize By Sample Factor Completed...")
  
  return(sample_raw_mat)  
}