#' create_multi_cohorts_metadata_for_limma
#'
#' It creates a metadata of two new conditions by pasting cohorts within cohort_a and cohort_b.
#' 
#' @param metadata The dataframe containing metadata information
#' @param cohort_col A metadata column where cohorts are present
#' @param cohort_a Vector of cohorts used as cohort_a
#' @param cohort_b Vector of cohorts used as cohort_b
#' @param new_colname New column where new conditions will be placed
#' @return dataframe, metadata with newly added column
#' @examples
#' create_multi_cohorts_metadata_for_limma(metadata, 'Cohort', c("Cohort1", "Cohort2"),  c("Cohort3", "Cohort4"), "multi_condition")
#' @export
create_multi_cohorts_metadata_for_limma <- function(metadata = NULL, cohort_col = NULL, 
                                                    cohort_a = NULL, cohort_b = NULL, 
                                                    new_colname = NULL){
  message("Create Multi Cohorts Metadata For Limma Started...")
  
  if (identical(metadata, NULL)){
    warning("The metadata is NULL")
    return (NULL)
  }
  
  if (identical(cohort_col, NULL)){
    warning("The cohort_col is NULL")
    return (NULL)  
  }
  
  if (identical(cohort_a, NULL)){
    warning("The cohort_a is NULL")
    return (NULL)  
  }
  
  if (identical(cohort_b, NULL)){
    warning("The cohort_b is NULL")
    return (NULL)  
  }
  
  if (identical(new_colname, NULL)){
    warning("The new_colname is NULL")
    return (NULL)  
  }
  
  common_a_and_b <- base::intersect(cohort_a, cohort_b)
  if (length(common_a_and_b) > 0){
    warning(paste0("The following cohorts are common in cohort_a and cohort_b: ", paste0(common_a_and_b, collapse = ", ")))
    return(NULL)
  }
  
  diff_cohort_a <- base::setdiff(cohort_a, metadata[[cohort_col]])
  if (length(diff_cohort_a) > 0){
    warning(paste0("The following cohorts in cohort_a are invalid: ", paste0(diff_cohort_a, collapse = ", ")))
    return(NULL)
  }
  
  diff_cohort_b <- base::setdiff(cohort_b, metadata[[cohort_col]])
  if (length(diff_cohort_b) > 0){
    warning(paste0("The following cohorts in cohort_b are invalid: ", paste0(diff_cohort_b, collapse = ", ")))
    return(NULL)
  }    
  
  if (identical(cohort_col, new_colname)){
    temp_cohort_col <- paste0(new_colname, "_old") 
    colnames(metadata)[which(names(metadata) == cohort_col)] <- temp_cohort_col
    cohort_col <- temp_cohort_col
  }
  
  new_colname <- as.character(new_colname)  
  metadata <- metadata[metadata[, cohort_col] %in% c(cohort_a, cohort_b), , drop = FALSE]  
  metadata[, new_colname] <- NA
  metadata[metadata[, cohort_col] %in% cohort_a, new_colname] <- paste(cohort_a, collapse = " - ")
  metadata[metadata[, cohort_col] %in% cohort_b, new_colname] <- paste(cohort_b, collapse = " - ")
  
  message("Create Multi Cohorts Metadata For Limma Completed...")
  
  return (metadata)
  
}