#' diff_exp_limma
#'
#' calls limma differential expression algorithms
#'
#' @param sample_raw_mat The sample matrix with samples in columns and ids as rownames.
#' @param metadata The dataframe containing metadata information
#' @param cohort_col A metadata column where cohorts are present
#' @param cohort_a Vector of cohorts used as cohort_a
#' @param cohort_b Vector of cohorts used as cohort_b
#' @param pval_adjust_method Provide pval adjust method ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). Also check p.adjust from stats.
#' @param log_flag Check if data is log transformed (TRUE) or not (FALSE). If not (FALSE) then do internally log2 transformation.
#' @return dataframe, If logFC >0, it implies abundance is greater in cohort_b, MaxExpr is calculated as the maximum of average of samples within cohorts.
#' @examples
#' diff_exp_limma(sample_raw_mat, metadata, 'Cohort', 'Cohort1', 'Cohort2')
#' @import limma dplyr
#' @export
diff_exp_limma <- function (sample_raw_mat = NULL, metadata = NULL, cohort_col = NULL, 
                            cohort_a = NULL, cohort_b = NULL, pval_adjust_method = "BH", log_flag = TRUE) {
  message("Calculate Differential Expression Limma Started...")
  require(limma)
  require(dplyr)
  
  if (identical(sample_raw_mat, NULL)){
    warning("The sample_raw_mat is NULL")
    return (NULL)  
  }
  
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
  
  if (identical(pval_adjust_method, NULL)){
    pval_adjust_method <- "BH"  
    warning("The pval_adjust_method is NULL, using 'BH' as default method.")  
  }    
  
  if (identical(log_flag, NULL)){
    log_flag <- TRUE  
    warning("The log_flag is NULL, using 'TRUE' as default method.")  
  }
  
  if (!(cohort_col %in% colnames(metadata))){
    warning(paste0(cohort_col, " is not present in metadata columns"))
    return (NULL)  
  }
  
  if (identical(cohort_col, "Comparison")){
    colnames(metadata)[which(names(metadata) == cohort_col)] <- "Cohort"  
    cohort_col <- "Cohort"
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
  
  metadata <- metadata[metadata[, cohort_col] %in% c(cohort_a, cohort_b), , drop = FALSE]
  metadata[, "Comparison"] <- NA
  metadata[metadata[, cohort_col] %in% cohort_a, "Comparison"] <- "A"
  metadata[metadata[, cohort_col] %in% cohort_b, "Comparison"] <- "B"
  
  common_samples <- base::intersect(metadata[, 1], colnames(sample_raw_mat))
  if (length(common_samples) < 1) {
    warning("No common samples in matrix and metadata")
    return(NULL)
  }
  
  diff_samples_metadata <- base::setdiff(metadata[, 1], colnames(sample_raw_mat))
  if (length(diff_samples_metadata) > 0) {
    warning(paste0("The following samples from metadata are not present in sample_raw_mat: ", paste0(diff_samples_metadata, collapse = ", ")))
  }    
  
  metadata <- metadata[!(metadata[, 1] %in% diff_samples_metadata), , drop = FALSE] 
  
  if (log_flag) {
    sample_raw_mat_log2 <- sample_raw_mat[, metadata[, 1], drop = FALSE]
  }
  else {
    sample_raw_mat_log2 <- log2(sample_raw_mat[, metadata[, 1], drop = FALSE])
  }
  
  cohort_a_samples <- metadata[metadata[, "Comparison"] == "A", ][, 1]
  cohort_b_samples <- metadata[metadata[, "Comparison"] == "B", ][, 1]
  if (length(cohort_a_samples) <= 1 | length(cohort_b_samples) <= 1) {
    warning("Since, your selected cohorts have no replicates in the data, you might want to change the cohorts or cohort condition and try again.")
    return(NULL)
  }
  
  limma_results_df <- NULL
  tryCatch(
    {
      condition <- metadata[, "Comparison"]
      design <- stats::model.matrix(~condition + 0)
      colnames(design) <- gsub("condition", "", colnames(design))
      contrast_matrix <- limma::makeContrasts(contrasts = c(paste("B", "A", sep = "-")), levels = design)
      fit <- limma::lmFit(sample_raw_mat_log2, design)
      fit <- limma::contrasts.fit(fit, contrast_matrix)
      fit <- limma::eBayes(fit)
      limma_results_df <- limma::topTable(fit, coef = paste("B", "A", sep = "-"), number = nrow(sample_raw_mat_log2), adjust.method = pval_adjust_method)
      limma_results_df <- limma_results_df[order(match(rownames(limma_results_df), rownames(sample_raw_mat_log2))), , drop = FALSE]
    },
    error = function(cond) {message(paste("\nCannot run limma, caused an error: ", cond))}
  )
  
  for (row_name in row.names(sample_raw_mat_log2)){
    tryCatch({  
      feature_df <- data.frame(metadata, stringsAsFactors = FALSE, check.names = FALSE)
      feature_df$value <- as.numeric(sample_raw_mat_log2[row_name, feature_df[, 1]])
      mean_df <- feature_df %>% dplyr::group_by(Comparison) %>% dplyr::summarise(mean = mean(value, na.rm = TRUE))
      limma_results_df[row_name, "MaxExpr"] <- max(mean_df$mean, na.rm = TRUE)
    },
    error = function(cond) {message(paste("Feature: ", row_name, "\nCannot calculate maximum of average of samples within cohorts, caused an error: ", cond))}    
    )
  }
  
  message("Calculate Differential Expression Limma Completed...")
  
  return(limma_results_df)    
}