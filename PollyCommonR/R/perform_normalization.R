#' perform_normalization
#'
#' It normalizes the data using some generic methods
#'
#' @param sample_raw_mat matrix/dataframe containing raw values
#' @param norm_type Select normalization type from the following methods:
#' "auto_scale": mean-centered and divided by the standard deviation of each variable;
#' "pareto_scale": mean-centered and divided by the square root of the standard deviation of each variable;
#' "mean_center_scale": mean-centered only;  
#' "range_scale": mean-centered and divided by the range of each variable;
#' "log_trans": log transformation;
#' "log10_trans": log10 transformation;
#' "log2_trans": log2 transformation;
#' "glog_trans": generalized logarithm transformation, tolerant to 0 and negative values;
#' "glog10_trans": generalized logarithm transformation with base 10, tolerant to 0 and negative values;
#' "glog2_trans": generalized logarithm transformation with base 2, tolerant to 0 and negative values;
#' "asinh_trans": ArcSinh transformation;  
#' @return The normalized sample_raw_mat
#' @examples
#' perform_normalization(sample_raw_mat, norm_type)
#' @export
perform_normalization = function (sample_raw_mat = NULL, norm_type = NULL){
  message("Perform Normalization Started...")
  
  if (identical(sample_raw_mat, NULL)){
    warning("The sample_raw_mat is NULL")
    return(NULL)
  }
  
  norm_types <- c("auto_scale", "pareto_scale", "mean_center_scale", "range_scale",
                  "log_trans", "log10_trans", "log2_trans", "glog_trans", "glog10_trans",
                  "glog2_trans", "asinh_trans")
  
  if (!(norm_type %in% norm_types)){
    warning(paste0("Please select valid norm_type from the following methods: ", paste0(sQuote(norm_types), collapse = ", ")))
  }
  
  # Scale
  # normalize to zero mean and unit variance
  auto_scale <- function(x){ (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)}
  
  # normalize to zero mean but variance/SE
  pareto_scale <- function(x){ (x - mean(x, na.rm = TRUE))/sqrt(sd(x, na.rm = TRUE))}
  
  # normalize to zero mean but variance/SE
  mean_center_scale <- function(x){ x - mean(x, na.rm = TRUE)}
  
  # normalize to zero mean but variance/SE
  range_scale <- function(x){
    if(max(x, na.rm = TRUE) == min(x, na.rm = TRUE)){ x}
    else{ (x - mean(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))}
  }
  
  # Transformation
  # generalized log, tolerant to 0 and negative values
  glog_trans <- function(x, min_val){ log((x + sqrt(x^2 + min_val^2))/2)}
  
  # generalized log with base 10, tolerant to 0 and negative values
  glog10_trans <- function(x, min_val){ log10((x + sqrt(x^2 + min_val^2))/2)}
  
  # generalized log with base 2, tolerant to 0 and negative values
  glog2_trans <- function(x, min_val){ log2((x + sqrt(x^2 + min_val^2))/2)}
  
  # arc sinh transformation
  asinh_trans <- function(x){ log(x + sqrt(x^2 + 1))}
  
  if(identical(norm_type, 'auto_scale')){
    sample_raw_mat <- t(apply(sample_raw_mat, 1, auto_scale))
  }
  else if(identical(norm_type, 'pareto_scale')){
    sample_raw_mat <- t(apply(sample_raw_mat, 1, pareto_scale))
  }
  else if(identical(norm_type, 'mean_center_scale')){ 
    sample_raw_mat <- t(apply(sample_raw_mat, 1, mean_center_scale))
  }
  else if(identical(norm_type, 'range_scale')){
    sample_raw_mat <- t(apply(sample_raw_mat, 1, range_scale))
  }
  else if(identical(norm_type, 'log_trans')){
    sample_raw_mat <- log(sample_raw_mat)
  }
  else if(identical(norm_type, 'log10_trans')){
    sample_raw_mat <- log10(sample_raw_mat)
  }
  else if(identical(norm_type, 'log2_trans')){
    sample_raw_mat <- log2(sample_raw_mat)
  }
  else if(identical(norm_type, 'glog_trans')){
    numeric_values <- as.numeric(unlist(sample_raw_mat))
    min_val <- min(abs(numeric_values[is.finite(numeric_values) & numeric_values != 0]))/10
    sample_raw_mat <- t(apply(sample_raw_mat, 1, glog_trans, min_val))
  }
  else if(identical(norm_type, 'glog10_trans')){
    numeric_values <- as.numeric(unlist(sample_raw_mat))
    min_val <- min(abs(numeric_values[is.finite(numeric_values) & numeric_values != 0]))/10
    sample_raw_mat <- t(apply(sample_raw_mat, 1, glog10_trans, min_val))
  }
  else if(identical(norm_type, 'glog2_trans')){
    numeric_values <- as.numeric(unlist(sample_raw_mat))
    min_val <- min(abs(numeric_values[is.finite(numeric_values) & numeric_values != 0]))/10
    sample_raw_mat <- t(apply(sample_raw_mat, 1, glog2_trans, min_val))
  }
  else if(identical(norm_type, 'asinh_trans')){
    sample_raw_mat <- t(apply(sample_raw_mat, 1, asinh_trans))
  }
  else {
    message("\nNo normalization performed")
  }
  
  sample_raw_mat <- as.data.frame(sample_raw_mat, stringsAsFactors = FALSE, check.names = FALSE)
  
  message("Perform Normalization Completed...")    
  return (sample_raw_mat)    
}