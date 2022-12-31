#' create_multi_cohorts_metadata_for_anova
#'
#' It adds the different combinations of columns after pasting them as new columns in metadata.
#' 
#' @param metadata The dataframe containing metadata information
#' @param cohort_col A vector of multiple columns from metadata used to make combinations of columns
#' @return dataframe, metadata with newly added columns
#' @examples
#' create_multi_cohorts_metadata_for_anova(metadata, cohort_col = c("Cohort", "Time"))
#' @export
create_multi_cohorts_metadata_for_anova <- function(metadata = NULL, cohort_col = NULL){
  message("Create Multi Cohorts Metadata For Anova Started...")  
  if (identical(metadata, NULL)){
    warning("The metadata is NULL")
    return (NULL)
  }
  
  if (!identical(class(metadata), "data.frame")){
    warning("The metadata is not a dataframe")
    return (NULL)
  }
  
  if (identical(cohort_col, NULL)){
    warning("The cohort_col is NULL")
    return (NULL)  
  }
  
  diff_cohort_col <- base::setdiff(cohort_col, colnames(metadata))  
  if (length(diff_cohort_col) > 0){
    warning(paste0("The following columns are not present in metadata: ", paste0(sQuote(diff_cohort_col), collapse = ", ")))
    return(NULL)  
  }  
  
  cohort_col <- unique(cohort_col) 
  if (length(cohort_col) < 2){
    warning("Please select more than two columns to add combinations of columns in metadata, returning actual metadata")  
  }
  else{
    tryCatch({   
      cohort_combinations <- list()
      for (cohort_ind in 2:length(cohort_col)){    
        combinations_df <- as.data.frame(gtools::combinations(length(cohort_col), cohort_ind, cohort_col, repeats.allowed = FALSE), stringsAsFactors = FALSE)
        for (comb_index in 1:nrow(combinations_df)){
          cohort_combinations <- c(cohort_combinations, list(as.character(combinations_df[comb_index, ])))
        }
      }
      
      for (cohort_comb in cohort_combinations){
        cohort_comb <- cohort_comb[order(match(cohort_comb, cohort_col))]
        metadata[, paste(cohort_comb, collapse = ":")] <- apply( metadata[ , cohort_comb] , 1 , paste , collapse = " - " ) 
      }
    },
    error = function(cond) {message(paste("\nCannot make combinations of cohort columns, caused an error: ", cond))}
    )
  }
  
  message("Create Multi Cohorts Metadata For Anova Completed...")
  
  return (metadata)  
  
}