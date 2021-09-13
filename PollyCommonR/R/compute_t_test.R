#' compute_t_test
#'
#' Perform t-test on dataframe/matrix where samples are in columns and features are in rows.
#'
#' @param data_mat A dataframe/matrix containing samples in columns and features in rows.
#' @param samples_a A vector of samples in cohort A.
#' @param samples_b A vector of samples in cohort B.
#' @param paired A logical indicating whether to perform paired t-test or not.
#' @param equal_var A logical variable indicating whether to treat the two variances as being equal. 
#' If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used.
#' @param nonpar A logical variable indicating whether to perform a non-parametric test i.e. wilcox.test (TRUE) or t.test.
#' @return A dataframe with calculated statistics.
#' @examples 
#' compute_t_test(data_mat, samples_a, samples_b)
#' @import stats
#' @export
compute_t_test <- function(data_mat = NULL, samples_a = NULL, samples_b = NULL, 
                           paired = FALSE, equal_var = TRUE, nonpar = FALSE, 
                           pval_adjust_method = "BH"){
  message("Compute T Test Started...")
  require(stats)
  
  if (identical(data_mat, NULL)){
    warning("The data_mat is NULL")
    return (NULL)  
  }
  
  if (!identical(class(data_mat), "data.frame") && !identical(class(data_mat), "matrix")){
    warning("The data_mat is not a dataframe/matrix, please provide valid norm_data")
    return(NULL) 
  }
  
  common_samples_a <- base::intersect(samples_a, colnames(data_mat))
  if (length(common_samples_a) < 1){
    warning("No common samples found between samples_a and data matrix")
    return(NULL)  
  }    
  
  common_samples_b <- base::intersect(samples_b, colnames(data_mat))  
  if (length(common_samples_b) < 1){
    warning("No common samples found between samples_b and data matrix")
    return(NULL)  
  }
  
  samples_diff <- base::setdiff(c(samples_a, samples_b), colnames(data_mat)) 
  if (length(samples_diff) > 0){
    warning(paste0("The following samples are not present in the data matrix :", paste0(samples_diff, collapse = ", ")))
  }
  
  if (identical(paired, TRUE)){
    if (!(length(common_samples_a) == length(common_samples_b))){
      warning("The samples_a and samples_b should have equal number of samples for paired comparison")
      return(NULL)    
    }
  }
  
  if (identical(nonpar, TRUE)){
    ttest <- function(x){ stats::wilcox.test(x = x[common_samples_a], y = x[common_samples_b], paired = paired, na.rm = TRUE)}
    test_stat_name <- "W"  
  }
  else {
    ttest <- function(x){ stats::t.test(x = x[common_samples_a], y = x[common_samples_b], paired = paired, var.equal = equal_var, na.rm = TRUE)}
    test_stat_name <- "t"  
  }
  
  ttest_res <- t(apply(data_mat, 1, function(x) {                                    
    test_temp <- try({ttest(x)})                                        
    if(class(test_temp) == "try-error") {
      return(c(NA, NA));
    }
    else{
      return(c(test_temp$statistic, test_temp$p.value));
    }
  }))
  
  ttest_res <- as.data.frame(ttest_res, stringsAsFactors = FALSE, check.names = FALSE)   
  colnames(ttest_res) <- c(test_stat_name, "P.Value")    
  ttest_res[, "adj.P.Val"] <- stats::p.adjust(ttest_res[, "P.Value"], method = pval_adjust_method)
  ttest_res <- data.frame(id = row.names(ttest_res), ttest_res, stringsAsFactors = FALSE, check.names = FALSE)  
  
  message("Compute T Test Completed...")
  
  return (ttest_res) 
}