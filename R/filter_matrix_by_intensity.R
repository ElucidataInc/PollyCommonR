#' filter_matrix_by_intensity
#'
#' Filter matrix by intensity present in at least x (in percent) samples across all samples or in at least one cohort. 
#'
#' @param sample_raw_mat A dataframe containing samples with raw values
#' @param algo An algo used to filter data by all samples (algo = "sample") or at least in one cohort (algo = "cohort")
#' @param metadata_df A dataframe with samples to cohort mapping
#' @param cohort_col A cohort column present in metadata
#' @param intensity_threshold The minimum intensity that should be present in samples
#' @param samples_threshold The number of samples in percentage where minimum intensity should be present
#' @return The filtered sample_raw_mat
#' @examples 
#' filter_matrix_by_intensity(sample_raw_mat, algo = "sample")
#' @import dplyr
#' @export
filter_matrix_by_intensity <- function(sample_raw_mat = NULL, algo = "sample", metadata_df = NULL, 
                                       cohort_col = NULL, intensity_threshold = NULL, 
                                       samples_threshold = NULL){
  message("Filter Matrix By Intensity Started...")
  require(dplyr)
  
  if (identical(sample_raw_mat, NULL)){
    warning("The sample_raw_mat is NULL")
    return(NULL) 
  }
  
  if (!identical(class(sample_raw_mat), "data.frame")){
    warning("The sample_raw_mat is not a dataframe, please provide valid sample_raw_mat")
    return(NULL) 
  }
  
  if (nrow(sample_raw_mat) < 1){
    warning("The sample_raw_mat is a blank dataframe, please provide valid sample_raw_mat")
    return(NULL) 
  }
  
  if (!(algo %in% c("sample", "cohort"))){
    warning("Please choose correct algo from 'sample' and 'cohort'")
    return (NULL)
  }
  
  if (algo %in% "cohort"){
    if (identical(metadata_df, NULL)){
      warning("The metadata_df is NULL")
      return(NULL) 
    }
    
    if (!identical(class(metadata_df), "data.frame")){
      warning("The metadata_df is not a dataframe, Please provide valid metadata_df")
      return(NULL) 
    }  
    
    if (nrow(metadata_df) < 1){
      warning("The metadata_df is a blank dataframe, please provide valid metadata_df")
      return(NULL) 
    }
    
    if (!(cohort_col %in% colnames(metadata_df))) {
      warning(c(cohort_col, " is not a valid cohort_col, please choose from metadata_df colnames"))
      return(NULL)
    }   
  }
  
  if (identical(intensity_threshold, NULL)){
    warning("The intensity_threshold is NULL")
    return(NULL) 
  }
  
  if (identical(samples_threshold, NULL)){
    warning("The samples_threshold is NULL")
    return(NULL) 
  }  
  
  intensity_threshold <- as.numeric(intensity_threshold)
  if (is.na(intensity_threshold)) {
    warning("The intensity_threshold is not a numeric value")
    return(NULL) 
  }
  
  samples_threshold <- as.numeric(samples_threshold)
  if (is.na(samples_threshold)) {
    warning("The samples_threshold is not a numeric value")
    return(NULL) 
  }
  
  features_vec <- vector()
  for (row_ind in row.names(sample_raw_mat)){ 
    if (algo %in% "sample"){ 
      feature_df <- as.data.frame(sample_raw_mat[row_ind, , drop = TRUE], stringsAsFactors = FALSE, check.names = FALSE)
      feature_df <- data.frame(intensity = t(feature_df), stringsAsFactors = FALSE, check.names = FALSE)
      feature_filtered <- dplyr::filter(feature_df, intensity >= intensity_threshold)
      feature_count_stats <- round(((nrow(feature_filtered) / nrow(feature_df)) * 100), 2)
      if (feature_count_stats >= samples_threshold){
        features_vec <- c(features_vec, row_ind)
      }
    }
    else if (algo %in% "cohort"){
      feature_df <- as.data.frame(sample_raw_mat[row_ind, , drop = TRUE], stringsAsFactors = FALSE, check.names = FALSE)
      feature_df <- data.frame(intensity = t(feature_df), stringsAsFactors = FALSE, check.names = FALSE)
      merged_sample_df <- base::merge(metadata_df, feature_df, by.x = 1, by.y = 0, sort = FALSE) 
      merged_sample_filtered <- dplyr::filter(merged_sample_df, intensity >= intensity_threshold)
      merged_sample_filtered <- dplyr::group_by_at(merged_sample_filtered, cohort_col) %>% summarise(count = dplyr::n())
      cohorts_count_df <- dplyr::group_by_at(merged_sample_df, cohort_col) %>% summarise(total_count = dplyr::n())
      feature_count_stats <- base::merge(merged_sample_filtered, cohorts_count_df, by = cohort_col, sort = FALSE)
      feature_count_stats$count_percent <- round(((feature_count_stats$count / feature_count_stats$total_count) * 100), digits = 2)
      feature_count_stats <- dplyr::filter(feature_count_stats, count_percent >= samples_threshold)
      if (nrow(feature_count_stats) >= 1){
        features_vec <- c(features_vec, row_ind)
      } 
    }   
  }
  
  sample_raw_mat <- sample_raw_mat[features_vec, , drop = FALSE]
  
  message("Filter Matrix By Intensity Completed...")
  
  return (sample_raw_mat)
}
