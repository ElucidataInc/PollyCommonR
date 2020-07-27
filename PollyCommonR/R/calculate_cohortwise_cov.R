#' calculate_cohortwise_cov
#'
#' calcualte coefficient of variation across all cohorts
#'
#' @param raw_matrix matrix/dataframe containing raw values.
#' @param metadata dataframe containing metadata information 
#' @param cohort_col A metadata column where cohorts are present
#' @return A dataframe with calculated cov
#' @examples
#' calculate_cohortwise_cov(raw_matrix, metadata, cohort_col='Cohort')
#' @import dplyr matrixStats
#' @export
calculate_cohortwise_cov <- function(raw_matrix, metadata, cohort_col='Cohort'){
  
  message("Calculate Coefficient of Variation Started...")
  require(dplyr)
  require(matrixStats)

  if (!(cohort_col %in% colnames(metadata))){
    warning(c(cohort_col, " column is not present in metadata"))
    return(NULL)
  }  

  raw_matrix_samples <- colnames(raw_matrix)
  metadata_samples <- as.character(metadata[, 1])
  common_samples <- intersect(metadata_samples, raw_matrix_samples)
  if (length(common_samples) == 0){
    warning("No common samples in matrix and metadata")
    return(NULL)
  }
  metadata <- dplyr::filter(metadata, !!(sym(colnames(metadata)[1])) %in% !!common_samples)
  id_cols_vec <- raw_matrix_samples[!raw_matrix_samples %in% common_samples]
  cohort_vec <- as.character(unique(metadata[[cohort_col]]))
  cohortwise_cov_list <- list()
  for (cohort in cohort_vec){
    interm_df = data.frame(raw_matrix[,id_cols_vec])
    interm_df[['id']] <- rownames(raw_matrix)
    interm_df[['cohort']] = cohort
    cohort_samples <- as.character(dplyr::filter(metadata, !!(sym(cohort_col)) %in% !!cohort)[[1]])
    cohortwise_sample_mat = as.matrix(raw_matrix[, cohort_samples])
    interm_df[['mean']] <- rowMeans(cohortwise_sample_mat, na.rm = TRUE)
    interm_df[['std']] <- matrixStats::rowSds(cohortwise_sample_mat, na.rm = TRUE)
    interm_df[['cv']] <- abs(interm_df[['std']]/interm_df[['mean']])*100
    cohortwise_cov_list[[cohort]] <- interm_df
  }
  cohortwise_cov_df <- dplyr::bind_rows(cohortwise_cov_list)
  
  message("Calculate Coefficient of Variation Completed...")
  
  return(cohortwise_cov_df)
}