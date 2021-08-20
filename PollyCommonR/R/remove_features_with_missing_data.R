#' remove_features_with_missing_data
#'
#' Remove features having missing values in samples above the specified cutoff  
#'
#' @param sample_raw_mat A dataframe/matrix containing samples in columns and features in rows
#' @param samples_cutoff The cutoff of the number of samples in percentage where missing values are allowed
#' @return The filtered sample_raw_mat
#' @examples 
#' remove_features_with_missing_data(sample_raw_mat, samples_cutoff)
#' @export
remove_features_with_missing_data <- function(sample_raw_mat = NULL, samples_cutoff = 50){
  message("Remove Features With Missing Data Started...")
  
  if (identical(sample_raw_mat, NULL)){
    warning("The sample_raw_mat is NULL")
    return(NULL) 
  }  
  
  if (!identical(class(sample_raw_mat), "data.frame") && !identical(class(sample_raw_mat), "matrix")){
    warning("The sample_raw_mat is not a dataframe/matrix, please provide valid sample_raw_mat")
    return(NULL) 
  }
  
  if (identical(samples_cutoff, NULL)){
    warning("The samples_cutoff is NULL, returning whole data")
  }
  else {
    samples_cutoff <- as.numeric(samples_cutoff)
    if (is.na(samples_cutoff)){
      warning("The samples_cutoff is not a numeric value") 
      return (NULL)
    }  
    correct_data_bool <- (apply(is.na(sample_raw_mat), 1, sum)/ncol(sample_raw_mat)) * 100 < samples_cutoff   
    missing_count <- length(correct_data_bool[correct_data_bool == FALSE])
    if (missing_count > 0){
      sample_raw_mat <- sample_raw_mat[correct_data_bool, ] 
      if (missing_count == 1){ features_str <- " feature is "}
      else { features_str <- " features are "}  
      
      message(paste0(missing_count, features_str, "removed from the data"))  
    }
    else {
      message("No features are removed from the data")
    }  
  }
  
  message("Remove Features With Missing Data Completed...")
  
  return (sample_raw_mat)  
}