#' compute_differential_expression
#'
#' It calculates the differential expression for cohort-B vs cohort-A or cohort-B/cohort-A comparison
#'
#' @param norm_data The dataframe/matrix with samples in columns and features in rows.
#' @param metadata The dataframe containing metadata information
#' @param cohort_col A metadata column where cohorts are present
#' @param cohort_a Vector of cohorts used as cohort_a
#' @param cohort_b Vector of cohorts used as cohort_b
#' @param algo Use limma or ttest
#' @param paired A logical indicating whether to use paired samples for ttest or fold change caculations.
#' @param equal_var A logical variable indicating whether to treat the two variances as being equal. 
#' If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used. 
#' It is only applicable with ttest algorithm.
#' @param nonpar A logical variable indicating whether to perform a non-parametric test i.e. wilcox.test (TRUE) or t.test. It is only applicable with ttest algorithm.
#' @param pval_adjust_method Provide pval adjust method ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). Also check p.adjust from stats.
#' @param fc_with_raw Use raw data (data without normalization) to calculate fold change.
#' @param raw_log_flag Check if raw data is log transformed (TRUE) or not (FALSE).
#' @param add_expr_stat Calculate MaxExpr and MinExpr (TRUE) where MaxExpr/MinExpr is calculated as the maximum/minimum of average of samples within cohorts. 
#' @param add_expr_data A logical variable indicating whether to add expression data to the differential expression results (TRUE)
#' @param expr_with_raw If expr_with_raw = TRUE then it will use raw data instead of norm data to add in the differential expression results and also calculate expression statistics
#' @param raw_data The dataframe/matrix with samples in columns and features in rows having same dimensions as norm data.
#' @return The differential expression results with state1 vs state2 (state1/state2) or cohort-B vs cohort-A (cohort-B/cohort-A) comparison, If logFC >0, it implies abundance is greater in cohort-B (state1).
#' @examples
#' compute_differential_expression(norm_data, metadata, 'Cohort', 'Cohort1', 'Cohort2')
#' @import limma dplyr reshape2
#' @export
compute_differential_expression <- function(norm_data = NULL, metadata = NULL, cohort_col = NULL, 
                                            cohort_a = NULL, cohort_b = NULL, algo = "limma",
                                            paired = FALSE, equal_var = TRUE, nonpar = FALSE, 
                                            pval_adjust_method = "BH", fc_with_raw = FALSE, raw_log_flag = FALSE, 
                                            add_expr_stat = TRUE, add_expr_data = FALSE, expr_with_raw = FALSE, 
                                            raw_data = NULL){
  
  message("Compute Differential Expression Started...")
  require(dplyr)
  require(reshape2)
  
  if (identical(norm_data, NULL)){
    warning("The norm_data is NULL")
    return (NULL)  
  }
  
  if (!identical(class(norm_data), "data.frame") && !identical(class(norm_data), "matrix")){
    warning("The norm_data is not a dataframe/matrix, please provide valid norm_data")
    return(NULL) 
  }    
  
  if (identical(fc_with_raw, TRUE) || (identical(add_expr_stat, TRUE) && identical(expr_with_raw, TRUE)) || (identical(add_expr_data, TRUE) && identical(expr_with_raw, TRUE))){
    if (identical(raw_data, NULL)){
      warning("The raw_data is NULL")
      return (NULL)
    }
    
    if (!identical(class(raw_data), "data.frame") && !identical(class(raw_data), "matrix")){
      warning("The raw_data is not a dataframe/matrix, please provide valid raw_data")
      return(NULL) 
    }
    
    if (!all(dim(norm_data) == dim(raw_data))){
      warning("The norm_data and raw_data have different dimensions (rows and columns)")
      return(NULL)    
    }
    
    if (!all(colnames(raw_data) %in% colnames(norm_data))){
      warning("Mismatch in columns of norm_data and raw_data")
      return(NULL)
    }
    
    if (!all(row.names(raw_data) %in% row.names(norm_data))){
      warning("Mismatch in rows of norm_data and raw_data")
      return(NULL)
    }
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
  
  norm_data <- as.data.frame(norm_data, stringsAsFactors = FALSE, check.names = FALSE)
  common_samples <- base::intersect(metadata[, 1], colnames(norm_data))
  if (length(common_samples) < 1) {
    warning("No common samples in matrix and metadata")
    return(NULL)
  }
  
  diff_samples_metadata <- base::setdiff(metadata[, 1], colnames(norm_data))
  if (length(diff_samples_metadata) > 0) {
    warning(paste0("The following samples from metadata are not present in norm_data: ", paste0(diff_samples_metadata, collapse = ", ")))
  }    
  
  metadata <- metadata[!(metadata[, 1] %in% diff_samples_metadata), , drop = FALSE]
  cohort_a_samples <- metadata[metadata[, "Comparison"] == "A", ][, 1]
  cohort_b_samples <- metadata[metadata[, "Comparison"] == "B", ][, 1]
  
  norm_data <- norm_data[, metadata[, 1], drop = FALSE]
  if (identical(fc_with_raw, TRUE)){
    raw_data <- as.data.frame(raw_data, stringsAsFactors = FALSE, check.names = FALSE)  
    raw_data <- raw_data[, metadata[, 1], drop = FALSE]
  }
  
  if (identical(algo, "ttest") || identical(fc_with_raw, TRUE)){
    if (identical(paired, TRUE)){
      if (length(cohort_a_samples) != length(cohort_b_samples)){
        warning("The cohort_a and cohort_b should have equal number of samples for paired comparison")
        return(NULL)    
      }  
    }
  }  
  
  diff_exp_data <- NULL
  tryCatch({    
    if (identical(algo, "ttest")){
      if (identical(fc_with_raw, TRUE)){
        fold_change_df <- PollyCommonR::compute_fold_change(data_mat = raw_data, samples_a = cohort_a_samples, samples_b = cohort_b_samples, paired = paired, log_flag = raw_log_flag)
      } else {
        fold_change_df <- PollyCommonR::compute_fold_change(data_mat = norm_data, samples_a = cohort_a_samples, samples_b = cohort_b_samples, paired = paired, log_flag = TRUE)  
      }
      ttest_res <- PollyCommonR::compute_t_test(data_mat = norm_data, samples_a = cohort_a_samples, samples_b = cohort_b_samples, pval_adjust_method = pval_adjust_method, paired = paired, equal_var = equal_var, nonpar = nonpar) 
      diff_exp_data <- base::merge(fold_change_df, ttest_res, by = "id", sort = FALSE)
    } else if (identical(algo, "limma")){
      limma_res <- PollyCommonR::diff_exp_limma(sample_raw_mat = norm_data, metadata = metadata, cohort_col = cohort_col, cohort_a = cohort_a, cohort_b = cohort_b, pval_adjust_method = pval_adjust_method, log_flag = TRUE)
      if (identical(fc_with_raw, TRUE)){
        fold_change_df <- PollyCommonR::compute_fold_change(data_mat = raw_data, samples_a = cohort_a_samples, samples_b = cohort_b_samples, paired = paired, log_flag = raw_log_flag)  
        limma_res = subset(limma_res, select = -c(logFC, AveExpr))  
        diff_exp_data <-  base::merge(fold_change_df, limma_res, by = "id", sort = FALSE)
      } else {
        if (identical(add_expr_stat, TRUE)){
          limma_res = subset(limma_res, select = -AveExpr)
        }
        diff_exp_data <- limma_res
      }
    } else {
      warning("Please select valid algo (ttest or limma)")  
    }
    
    comparison_df <- data.frame(id = diff_exp_data$id, state1 = paste(cohort_b, collapse = " - "), state2 = paste(cohort_a, collapse = " - "), stringsAsFactors = FALSE, check.names = FALSE)  
    diff_exp_data <-  base::merge(comparison_df, diff_exp_data, by = "id", sort = FALSE)  
    row.names(diff_exp_data) <- diff_exp_data$id
  },
  error = function(cond) {message(paste("\nCannot calculate differential expression, caused an error: ", cond))}
  )
  
  if (identical(add_expr_stat, TRUE)){
    if (identical(expr_with_raw, TRUE)){
      long_sample_raw_mat <- raw_data    
    } else { long_sample_raw_mat <- norm_data}
    expr_stat_df <- NULL  
    tryCatch({
      long_sample_raw_mat$id <- row.names(long_sample_raw_mat)
      long_sample_raw_mat <- reshape2::melt(long_sample_raw_mat, id.vars = "id")
      long_sample_raw_mat <- base::merge(long_sample_raw_mat, metadata, by.x = "variable", by.y = 1)
      sample_mat_mean_df <- long_sample_raw_mat %>% dplyr::group_by_at(c("id", cohort_col)) %>% dplyr::summarise(Mean = mean(value, na.rm = TRUE))
      average_expr_df <- long_sample_raw_mat %>% dplyr::group_by_at("id") %>% dplyr::summarise(AveExpr = mean(value, na.rm = TRUE))
      cohort_wise_expr_df <- sample_mat_mean_df %>% dplyr::group_by_at("id") %>% dplyr::summarise(MaxExpr = max(Mean, na.rm = TRUE), MinExpr = min(Mean, na.rm = TRUE))
      expr_stat_df <- base::merge(average_expr_df, cohort_wise_expr_df, by= "id", sort = FALSE)   
    },
    error = function(cond) {message(paste("\nCannot calculate maximum and minimum of average of samples within cohorts, caused an error: ", cond))}    
    )  
    
    if (!identical(diff_exp_data, NULL) && !identical(expr_stat_df, NULL)){
      tryCatch({
        diff_exp_data <- base::merge(diff_exp_data, expr_stat_df, by= "id", sort = FALSE)
        row.names(diff_exp_data) <- diff_exp_data$id
      },
      error = function(cond) {message(paste("\nCannot merge differential expression and expression stat dataframes, caused an error: ", cond))}    
      )
    }
  }
  
  if (identical(add_expr_data, TRUE)){
    if (identical(expr_with_raw, TRUE)){
      expr_data <- raw_data    
    } else { expr_data <- norm_data}
    
    tryCatch({
      expr_data <- expr_data[, c(cohort_b_samples, cohort_a_samples), drop = FALSE]
      expr_data$id <- row.names(expr_data)
      if (!identical(diff_exp_data, NULL) && !identical(expr_data, NULL)){
        diff_exp_data <- base::merge(diff_exp_data, expr_data, by= "id", sort = FALSE)
        row.names(diff_exp_data) <- diff_exp_data$id
      }
    },
    error = function(cond) {message(paste("\nCannot merge differential expression and expression data dataframes, caused an error: ", cond))}    
    )   
  }
  
  message("Compute Differential Expression Completed...")
  
  return(diff_exp_data)
}