#' remove_features_with_constant_data
#'
#' Remove features having contant values in all samples
#'
#' @param sample_raw_mat A dataframe/matrix containing samples in columns and features in rows
#' @return The filtered sample_raw_mat
#' @examples 
#' remove_features_with_constant_data(sample_raw_mat)
#' @export
remove_features_with_constant_data <- function(sample_raw_mat = NULL){
  message("Remove Features With Constant Data Started...")
  
  if (identical(sample_raw_mat, NULL)){
    warning("The sample_raw_mat is NULL")
    return(NULL) 
  }  
  
  if (!identical(class(sample_raw_mat), "data.frame") && !identical(class(sample_raw_mat), "matrix")){
    warning("The sample_raw_mat is not a dataframe/matrix, please provide valid sample_raw_mat")
    return(NULL) 
  }
  
  var_row <- apply(sample_raw_mat, 1, stats::var, na.rm = TRUE)
  const_row <- (var_row == 0 | is.na(var_row))
  const_num <- sum(const_row, na.rm = TRUE)  
  
  if(const_num > 0){
    sample_raw_mat <- sample_raw_mat[!const_row, , drop = FALSE]         
    if (const_num == 1){ features_str <- " feature is "}
    else { features_str <- " features are "}  
    message(paste0(const_num, features_str, "removed from the data")) 
  }
  else {
    message("No features are removed from the data")
  }
  
  message("Remove Features With Constant Data Completed...")
  
  return (sample_raw_mat)    
}