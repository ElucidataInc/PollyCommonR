#' impute_missing_data
#'
#' Impute missing data using the following methods: 
#' "exclude": Exclude features having missing values; 
#' "replace_by_zero": Replace missing values by zero; 
#' "feature_lod": Replace missing values by limit of detection (LoD, 1/5 of the minimum positive value of each feature); 
#' "feature_min": Replace missing values by 1/2 of the minimum values of each feature row; 
#' "feature_mean": Replace missing values by the mean values of each feature row; 
#' "feature_median": Replace missing values by the median values of each feature row; 
#' "feature_knn": Estimate missing values by k-nearest neighbours based on similar features - KNN (feature-wise); 
#' "sample_knn": Estimate missing values by k-nearest neighbours based on similar samples - KNN (sample-wise); 
#' "bpca": Estimate missing values by Bayesian PCA (BPCA) method; 
#' "ppca": Estimate missing values by probabilistic PCA (PPCA) method; 
#' "svdImpute": Estimate missing values by Singular Value Decomposition (SVD) method;
#'
#' @param sample_raw_mat A dataframe/matrix containing samples in columns and features in rows
#' @param method The method used to impute missing data
#' @return The dataframe with imputated missing data
#' @examples 
#' impute_missing_data(sample_raw_mat, method)
#' @export
impute_missing_data <- function(sample_raw_mat = NULL, method = NULL){
  message("Impute Missing Data Started...")
  
  if (identical(sample_raw_mat, NULL)){
    warning("The sample_raw_mat is NULL")
    return(NULL) 
  }  
  
  if (!identical(class(sample_raw_mat), "data.frame") && !identical(class(sample_raw_mat), "matrix")){
    warning("The sample_raw_mat is not a dataframe/matrix, please provide valid sample_raw_mat")
    return(NULL) 
  }
  
  if (identical(method, NULL)){
    warning("The method is NULL")
    return(NULL) 
  }
  
  all_methods <- c("exclude", "replace_by_zero", "feature_lod", "feature_min", "feature_mean", 
                   "feature_median", "feature_knn", "sample_knn", "bpca", "ppca", "svdImpute")  
  if (!(method %in% all_methods)){
    warning(paste0("Please select valid method from : ", paste(sQuote(all_methods), collapse = ", ")))
  }    
  
  sample_raw_mat <- as.matrix(sample_raw_mat)  
  if(identical(method, "exclude")){
    correct_data_bool <- apply(is.na(sample_raw_mat), 1, sum) == 0
    sample_raw_mat <- sample_raw_mat[correct_data_bool, , drop = FALSE]
  }
  else if(identical(method, "replace_by_zero")){
    sample_raw_mat <- t(apply(sample_raw_mat, 1, function(x){
      if(sum(is.na(x)) > 0){ x[is.na(x)] <- 0}
      return (x)
    }))
  }
  else if(identical(method, "feature_lod")){
    sample_raw_mat <- t(apply(sample_raw_mat, 1, function(x) {
      lod <- min(x[x > 0], na.rm = TRUE)/5
      x[x == 0 | is.na(x)] <- lod
      return(x)
    }))  
  }
  else if(identical(method, "feature_min")){
    sample_raw_mat <- t(apply(sample_raw_mat, 1, function(x){
      if(sum(is.na(x)) > 0){ x[is.na(x)] <- min(x, na.rm = TRUE)/2}
      return (x)
    }))
  }    
  else if(identical(method, "feature_mean")){
    sample_raw_mat <- t(apply(sample_raw_mat, 1, function(x){
      if(sum(is.na(x)) > 0){ x[is.na(x)] <- mean(x, na.rm = TRUE)}
      return (x)
    }))
  }
  else if(identical(method, "feature_median")){
    sample_raw_mat <- t(apply(sample_raw_mat, 1, function(x){
      if(sum(is.na(x)) > 0){ x[is.na(x)] <- median(x, na.rm = TRUE)}
      return (x)
    }))
  }
  else if(identical(method, "feature_knn")){
    sample_raw_mat <- impute::impute.knn((sample_raw_mat))$data
  }
  else if(identical(method, "sample_knn")){
    sample_raw_mat <- t(impute::impute.knn(t(sample_raw_mat))$data)
  }
  else if(method %in% c("bpca", "ppca", "svdImpute")){
    sample_raw_mat <- t(pcaMethods::pca(t(sample_raw_mat), nPcs =5, method = method, center=T)@completeObs)
  }
  else {
    message("\nNo imputation of missing data performed")
  }    
  
  sample_raw_mat <- as.data.frame(sample_raw_mat, stringsAsFactors = FALSE, check.names = FALSE)
  message("Impute Missing Data Completed...")
  
  return (sample_raw_mat)  
}