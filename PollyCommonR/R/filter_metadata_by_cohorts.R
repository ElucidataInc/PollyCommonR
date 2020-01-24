#' filter_metadata_by_cohorts
#'
#' performs pca on the sample raw matrix
#'
#' @param metadata_df A dataframe with samples to cohort mapping.
#' @param condition A metadata column where cohorts are present.
#' @param selected_cohorts A vector of cohorts
#' @return A dataframe of filtered metadata.
#' @examples 
#' filter_metadata_by_cohorts(metadata_df, 'Cohort', c("Cohort1, Cohort2"))
#' @import dplyr
#' @export
filter_metadata_by_cohorts <- function(metadata, condition = 'Cohort', selected_cohorts = NULL){
  message("Filter Metadata Started...")
  require(dplyr)
  
  filtered_metadata_df <- metadata
  if (!identical(selected_cohorts, NULL)){
    common_cohorts <- intersect(selected_cohorts, metadata[,condition])
    if (!(length(common_cohorts) == 0)){
      filtered_metadata_df <- dplyr::filter(metadata, .data[[condition]] %in% common_cohorts)
    }else{
      message("No common cohorts, returning original metadata!")
    }
  }
  
  message("Filter Metadata Completed...")
  
  return (filtered_metadata_df)
}