#' one_way_anova_on_matrix
#'
#' performs pca on the sample raw matrix
#'
#' @param sample_raw_mat A dataframe containing samples with raw values.
#' @param metadata_df A dataframe with samples to cohort mapping.
#' @param cohort_col A metadata column where cohorts are present.
#' @return A dataframe with two extra columns (calculated F.Value and P.Value).
#' @examples 
#' one_way_anova_on_matrix(sample_raw_mat, metadata_df, cohort_col = "Cohort")
#' @import dplyr stats
#' @export
one_way_anova_on_matrix <- function(sample_raw_mat = NULL, metadata_df = NULL, cohort_col = "Cohort"){
  message("One Way Anova On Matrix Started...")
  
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
  
  if (!(cohort_col %in% colnames(metadata_df))){
    warning(c(cohort_col, " column is not present in metadata"))
    return(NULL)
  }
  if (length(unique(metadata_df[[cohort_col]])) < 2) {
    warning("The number of cohorts should be greater than 2")
    return(NULL)
  }
  metadata_sample <- metadata_df[, 1]
  raw_intensity_cols <- colnames(sample_raw_mat)
  sample_cols <- intersect(metadata_sample, raw_intensity_cols)
  if (length(sample_cols) == 0) {
    message("No common samples found, please provide valid data")
    return(NULL)
  }
  identifier_cols <- raw_intensity_cols[!(raw_intensity_cols %in% sample_cols)]
  if (length(identifier_cols) == 1) {
    identifier_df <- data.frame(sample_raw_mat[, identifier_cols, drop = FALSE])
    names(identifier_df) <- identifier_cols
  }
  else {
    identifier_df <- sample_raw_mat[, identifier_cols, drop = FALSE]
  }
  sample_intensity_mat <- sample_raw_mat[, sample_cols, drop = FALSE]
  sample_df <- data.frame(Sample = colnames(sample_intensity_mat))
  anova_input_df <- merge(sample_df, metadata_df, by = 1)
  anova_input_df[, cohort_col] <-  as.character(anova_input_df[[cohort_col]])
  anova_results <- apply(sample_intensity_mat, 1, function(x) {
    anova_input_df$value <- x
    anova_input_df <- anova_input_df[apply(anova_input_df, 1, function(x) is.finite(as.numeric(x[['value']]))), , drop = FALSE]
    F_val <- NA
    p_val <- NA
    if (nrow(anova_input_df) > 0){
      frm <- paste("value", cohort_col, sep = "~")
      anv1 <- stats::lm(stats::formula(frm), anova_input_df)
      a <- aov(anv1)
      if (all(c("F value", "Pr(>F)") %in% colnames(summary(a)[[1]]))) {
        F_val = summary(a)[[1]]["F value"][cohort_col, "F value"]
        p_val = summary(a)[[1]]["Pr(>F)"][[1]][1]
      }
    }
    return(list(F_statistic = F_val, p_value = p_val))
    
  })
  
  anova_results_df <- stats::setNames(
    data.frame(sapply(anova_results, function(x) x$F_statistic),
               sapply(anova_results, function(x) x$p_value)), 
    c("F.Value", "P.Value"))
  
  combined_anova_results_df <- dplyr::bind_cols(identifier_df, anova_results_df)   
  combined_anova_results_df <- data.frame(id = rownames(anova_results_df), combined_anova_results_df, stringsAsFactors = FALSE)
  combined_anova_results_df <- combined_anova_results_df[rowSums(is.na(combined_anova_results_df)) == 0, , drop = FALSE]
  
  if (nrow(combined_anova_results_df) < 1){
    warning("Each cohorts should have more than one sample")
    return(NULL)
  }
  
  message("One Way Anova On Matrix Completed...")
  return(combined_anova_results_df)
}