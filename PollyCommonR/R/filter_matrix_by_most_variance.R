#' filter_matrix_by_most_variance
#'
#' Filter matrix by most variance
#'
#' @param sample_raw_mat A dataframe containing samples with raw values.
#' @param top_n Filter top n number of ids'.
#' @param algo An algorithm used to sort ids', options are anova and mad.
#' @param metadata_df A dataframe with samples to cohort mapping.
#' @param cohort_col A metadata column where cohorts are present.
#' @return The filtered sample_raw_mat.
#' @examples 
#' filter_matrix_by_most_variance(sample_raw_mat, top_n = NULL, algo = "anova", metadata_df, cohort_col = "Cohort")
#' @import dplyr stats
#' @export
filter_matrix_by_most_variance <- function(sample_raw_mat = NULL, top_n = NULL, algo = "mad", 
                                           metadata_df = NULL, cohort_col = NULL){
  message("Filter Matrix By Most Variance Started...")
  
  if (!identical(class(sample_raw_mat), "data.frame")){
    warning("The sample_raw_mat is not a dataframe, please provide valid sample_raw_mat")
    return(NULL) 
  }
  
  if (nrow(sample_raw_mat) < 1){
    warning("The sample_raw_mat is a blank dataframe, please provide valid sample_raw_mat")
    return(NULL) 
  }
  
  if (identical(top_n, NULL)){
    warning("The top_n parameter is NULL")
    return (NULL)
  }
  
  if (!(algo %in% c("anova", "mad"))){
    warning("Please choose correct algo from anova and mad")
    return (NULL)
  }
  
  filtered_data <- sample_raw_mat
  if (algo %in% "anova"){
    
    if (!identical(class(metadata_df), "data.frame")){
      warning("The metadata_df is not a dataframe, Please provide valid metadata_df")
      return(NULL) 
    }  
    
    if (nrow(metadata_df) < 1){
      warning("The metadata_df is a blank dataframe, please provide valid metadata_df")
      return(NULL) 
    }
    
    if (!(cohort_col %in% colnames(metadata_df))){
      warning(c(cohort_col, " column is not present in metadata"))
      return(NULL)
    }
    
    if (length(unique(metadata_df[[cohort_col]])) < 2) {
      warning("The number of cohorts should be greater than 2")
      return(NULL)
    }     
    
    score_results <- PollyCommonR::one_way_anova_on_matrix(sample_raw_mat, metadata_df, cohort_col)
    if (!identical(score_results, NULL)){
      score_results$Score_FP <- (score_results$F.Value) * log10(1/score_results$P.Value)
      score_results <- score_results[order(-score_results$Score_FP), , drop = FALSE]
      score_results <- head(score_results, top_n)
      filtered_data <- sample_raw_mat[score_results$id, , drop = FALSE]
    }    
  }
  else if (algo %in% "mad"){                             
    mad_score <- apply(sample_raw_mat, 1, function(x) {
      x <- x[is.finite(as.numeric(x))]
      mad_var <- stats::mad(x, na.rm = TRUE)
      return (mad_var)
    })
    mad_score <- mad_score[is.finite(mad_score)]
    if (length(mad_score) > 0){
      score_results <- data.frame(id = names(mad_score), mad_score, check.names = FALSE, stringsAsFactors = FALSE)                 
      score_results <- score_results[order(-score_results$mad_score), , drop = FALSE]        
      score_results <- head(score_results, top_n)
      filtered_data <- sample_raw_mat[score_results$id, , drop = FALSE]
    }
  }    
  
  message("Filter Matrix By Most Variance Completed...")
  
  return(filtered_data)
}