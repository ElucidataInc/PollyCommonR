#' compute_anova
#'
#' performs n way anova test on the sample raw matrix
#'
#' @param sample_raw_mat A dataframe containing samples with raw values.
#' @param metadata_df A dataframe with samples to cohort mapping.
#' @param cohort_col A vector of metadata columns used for n way anova.
#' @return A dataframe with calculated F.Value and P.Value.
#' @examples 
#' compute_anova(sample_raw_mat, metadata_df, cohort_col = "Cohort")
#' @import stats stringr
#' @export
compute_anova <- function(sample_raw_mat = NULL, metadata_df = NULL, cohort_col = "Cohort"){
  message("Compute Anova On Matrix Started...")
  require(stats)
  require(stringr)
  
  if (!identical(class(sample_raw_mat), "data.frame")){
    warning("The sample_raw_mat is not a dataframe, please provide valid sample_raw_mat")
    return(NULL) 
  }
  
  if (!identical(class(metadata_df), "data.frame")){
    warning("The metadata_df is not a dataframe, Please provide valid metadata_df")
    return(NULL) 
  }
  
  if (nrow(sample_raw_mat) < 1){
    warning("The sample_raw_mat is a blank dataframe, please provide valid sample_raw_mat")
    return(NULL) 
  }
  
  if (nrow(metadata_df) < 1){
    warning("The metadata_df is a blank dataframe, please provide valid metadata_df")
    return(NULL) 
  }
  
  metadata_sample <- metadata_df[, 1]
  raw_intensity_cols <- colnames(sample_raw_mat)
  sample_cols <- base::intersect(metadata_sample, raw_intensity_cols)
  
  if (length(sample_cols) == 0) {
    message("No common samples found, please provide valid data")
    return(NULL)
  }
  
  metadata_df <- metadata_df[metadata_df[, 1] %in% sample_cols, , drop = FALSE]
  diff_cohort_col <- base::setdiff(cohort_col, colnames(metadata_df))  
  if (length(diff_cohort_col) > 0){
    warning(paste0("The following columns are not present in metadata: ", paste0(diff_cohort_col, collapse = ", ")))
    return(NULL)
  }
  
  if (nrow(unique(metadata_df[, cohort_col, drop = FALSE])) < 2) {
    warning("The number of cohorts for common samples should be greater than or equal to 2")
    return(NULL)
  }
  
  identifier_cols <- raw_intensity_cols[!(raw_intensity_cols %in% sample_cols)]
  if (length(identifier_cols) > 0) {
    identifier_df <- data.frame(id = row.names(sample_raw_mat), sample_raw_mat[, identifier_cols, drop = FALSE], stringsAsFactors = FALSE, check.names = FALSE)
  }
  else {
    identifier_df <- data.frame()
  }
  
  sample_intensity_mat <- sample_raw_mat[, sample_cols, drop = FALSE]  
  sample_df <- data.frame(Sample = colnames(sample_intensity_mat), stringsAsFactors = FALSE, check.names = FALSE)
  anova_input_df <- merge(sample_df, metadata_df, by = 1, sort = FALSE)
  anova_results_df <- data.frame(stringsAsFactors = FALSE, check.names = FALSE)    
  for (row_name in row.names(sample_intensity_mat)){
    tryCatch({  
      anova_input_df$value <- as.numeric(sample_intensity_mat[row_name, anova_input_df$Sample])
      anova_input_df <- anova_input_df[apply(anova_input_df, 1, function(x) is.finite(as.numeric(x[['value']]))), , drop = FALSE]                                       
      frm <- paste("value", paste(cohort_col, collapse = " * "), sep = " ~ ")
      anv_lm <- stats::lm(stats::formula(frm), anova_input_df)
      aov_obj <- stats::aov(anv_lm)
      anova_r <- summary(aov_obj)[[1]]
      anova_r <- data.frame(id = row_name, interaction = stringr::str_trim(row.names(anova_r)), anova_r[, c("F value", "Pr(>F)")], stringsAsFactors = FALSE, check.names = FALSE)
      anova_r <- anova_r[!stringr::str_trim(row.names(anova_r)) %in% c("Residuals"), , drop = FALSE]
      colnames(anova_r) <- c("id", "interaction", "F.Value", "P.Value")
      row.names(anova_r) <- NULL                                                                    
      anova_results_df <- rbind(anova_results_df, anova_r)
    }, 
    error = function(cond) {message(paste("\nCannot run anova, caused an error: ", cond))}
    ) 
  }
  
  combined_anova_results_df <- data.frame()
  tryCatch({                                   
    if (nrow(anova_results_df) > 0){ 
      if (nrow(identifier_df) > 0){
        combined_anova_results_df <- merge(identifier_df, anova_results_df, by = 'id', sort = FALSE) 
      }
      else{
        combined_anova_results_df <- anova_results_df
      }
      if (all(c("F.Value", "P.Value") %in% colnames(combined_anova_results_df))){
        combined_anova_results_df <- combined_anova_results_df[rowSums(is.na(combined_anova_results_df[, c("F.Value", "P.Value")])) == 0, , drop = FALSE]
      }
    }
  }, 
  error = function(cond) {message(paste("\nCannot create anova result dataframe, caused an error: ", cond))}
  )
  
  if (nrow(combined_anova_results_df) < 1){
    warning("Unable to perform anova test on this dataset")
    return(NULL)
  }
  
  message("Compute Anova On Matrix Completed...")
  
  return(combined_anova_results_df)
}