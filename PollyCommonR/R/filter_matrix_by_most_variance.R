#' filter_matrix_by_most_variance
#'
#' Filter matrix by most variance
#'
#' @param sample_raw_mat A dataframe containing samples with raw values.
#' @param top_n Filter top n number of ids'.
#' @param algo An algorithm used to sort ids', options are anova and mad.
#' @param metadata_df A dataframe with samples to cohort mapping.
#' @param cohort_col A vector of metadata columns used in anova test.
#' @param remove_duplicates remove duplicates using some column id from row descriptors.
#' @param row_desc A dataframe of row descriptors of the matrix.
#' @param id_col The id column of row descriptors dataframe used to remove duplicates.
#' @return The filtered sample_raw_mat.
#' @examples 
#' filter_matrix_by_most_variance(sample_raw_mat, top_n = NULL, algo = "anova", metadata_df, cohort_col = "Cohort")
#' @import dplyr stats
#' @export
filter_matrix_by_most_variance <- function(sample_raw_mat = NULL, top_n = NULL, algo = "mad", 
                                           metadata_df = NULL, cohort_col = NULL, 
                                           remove_duplicates = TRUE, row_desc = NULL, 
                                           id_col = NULL){
  message("Filter Matrix By Most Variance Started...")
  
  if (!identical(class(sample_raw_mat), "data.frame")){
    warning("The sample_raw_mat is not a dataframe, please provide valid sample_raw_mat")
    return(NULL) 
  }
  
  if (nrow(sample_raw_mat) < 1){
    warning("The sample_raw_mat is a blank dataframe, please provide valid sample_raw_mat")
    return(NULL) 
  }
  
  if (!(algo %in% c("anova", "mad"))){
    warning("Please choose correct algo from anova and mad")
    return (NULL)
  }
  
  if (identical(top_n, NULL) && identical(remove_duplicates, FALSE)){
    warning("The top_n parameter is NULL and remove_duplicates parameter is FALSE")
    return (NULL)
  }
  
  if (!identical(top_n, NULL) && is.na(as.numeric(top_n))){
    warning("The top_n should be a numberic value")
    return(NULL)
  }    
  
  if (remove_duplicates){
    if (!identical(class(row_desc), "data.frame")){
      warning("The row_desc is not a dataframe, please provide valid row_desc")
      return(NULL) 
    }   
    
    if (!all(row.names(sample_raw_mat) %in% row.names(row_desc))){
      warning("Not all rownames of sample_raw_mat are present in rownames of row_desc")
      return (NULL)  
    }
    
    if (identical(id_col, NULL)){
      warning("The id_col parameter is NULL")
      return (NULL)
    }      
    
    if (!(id_col %in% colnames(row_desc))){
      warning(c(id_col, " column is not present in row_desc"))
      return(NULL)
    }  
    
  }    
  
  score_results <- data.frame()  
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
    
    score_results <- PollyCommonR::compute_anova(sample_raw_mat, metadata_df, cohort_col)
    if (!identical(score_results, NULL)){
      score_results$Score_FP <- (score_results$F.Value) * log10(1/score_results$P.Value)
      score_results <- score_results[order(-score_results$Score_FP), , drop = FALSE]
      score_results <- score_results[!(duplicated(score_results$id)), , drop = FALSE]  
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
    }
  }
  
  if (identical(score_results, NULL) || nrow(score_results) < 1){      
    warning("Cannot filter data by most variance")
    return(NULL)
  }  
  
  if (remove_duplicates){ 
    score_results <- merge(score_results, row_desc, by.x = "id", by.y = 0, sort = FALSE)
    score_results <- score_results[!(duplicated(score_results[[id_col]])), , drop = FALSE]
  }
  
  if (!identical(top_n, NULL)){
    score_results <- head(score_results, top_n)   
  }
  
  filtered_data <- sample_raw_mat[unique(score_results$id), , drop = FALSE]
  
  message("Filter Matrix By Most Variance Completed...")
  
  return(filtered_data)
}