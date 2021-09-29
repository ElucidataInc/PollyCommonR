#' compute_anova
#'
#' performs n way anova test on the sample raw matrix
#'
#' @param norm_data The dataframe/matrix with samples in columns and features in rows.
#' @param metadata The dataframe with samples to cohort mapping.
#' @param cohort_col A vector of metadata columns used for n way anova.
#' @param expr_stat Calculate MaxExpr and MinExpr (TRUE) where MaxExpr/MinExpr is calculated as the maximum/minimum of average of samples within cohorts.
#' @param expr_with_raw Use raw data (data without normalization) to calculate expression statistics.
#' @param raw_data The dataframe/matrix with samples in columns and features in rows having same dimensions as norm data.
#' @return A dataframe with calculated F.Value and P.Value.
#' @examples 
#' compute_anova(norm_data, metadata, cohort_col)
#' @import dplyr stats stringr
#' @export
compute_anova <- function(norm_data = NULL, metadata = NULL, cohort_col = "Cohort",
                          expr_stat = TRUE, expr_with_raw = FALSE, raw_data = NULL){
  message("Compute Anova Started...")
  require(dplyr)
  require(stats)
  require(stringr)
  
  if (identical(norm_data, NULL)){
    warning("The norm_data is NULL")
    return (NULL)  
  }
  
  if (!identical(class(norm_data), "data.frame") && !identical(class(norm_data), "matrix")){
    warning("The norm_data is not a dataframe/matrix, please provide valid norm_data")
    return(NULL) 
  } 
  
  if (nrow(norm_data) < 1){
    warning("The norm_data is a blank dataframe, please provide valid norm_data")
    return(NULL) 
  }
  
  if (identical(expr_stat, TRUE) && identical(expr_with_raw, TRUE)){
    if (identical(raw_data, NULL)){
      warning("The raw_data is NULL")
      return (NULL)
    }
    
    if (!identical(class(raw_data), "data.frame") && !identical(class(raw_data), "matrix")){
      warning("The raw_data is not a dataframe/matrix, please provide valid raw_data")
      return(NULL) 
    }
    
    if (!all(dim(norm_data) == dim(raw_data))){
      warning("The norm_data and raw_data have different dimentions (rows and columns)")
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
  
  if (!identical(class(metadata), "data.frame")){
    warning("The metadata is not a dataframe, Please provide valid metadata")
    return(NULL) 
  }
  
  if (nrow(metadata) < 1){
    warning("The metadata is a blank dataframe, please provide valid metadata")
    return(NULL) 
  }
  
  metadata_sample <- metadata[, 1]
  norm_data_cols <- colnames(norm_data)
  sample_cols <- base::intersect(metadata_sample, norm_data_cols)
  
  if (length(sample_cols) == 0) {
    message("No common samples found, please provide valid data")
    return(NULL)
  }
  
  metadata <- metadata[metadata[, 1] %in% sample_cols, , drop = FALSE]
  diff_cohort_col <- base::setdiff(cohort_col, colnames(metadata))  
  if (length(diff_cohort_col) > 0){
    warning(paste0("The following columns are not present in metadata: ", paste0(diff_cohort_col, collapse = ", ")))
    return(NULL)
  }
  
  if (nrow(unique(metadata[, cohort_col, drop = FALSE])) < 2) {
    warning("The number of cohorts for common samples should be greater than or equal to 2")
    return(NULL)
  }
  
  anova_cohort_cols <- cohort_col
  names(anova_cohort_cols) <- make.names(anova_cohort_cols)
  names(metadata)[match(unname(anova_cohort_cols), names(metadata))] <- names(anova_cohort_cols)
  
  identifier_cols <- norm_data_cols[!(norm_data_cols %in% sample_cols)]
  if (length(identifier_cols) > 0) {
    identifier_df <- data.frame(id = row.names(norm_data), norm_data[, identifier_cols, drop = FALSE], stringsAsFactors = FALSE, check.names = FALSE)
  }
  else {
    identifier_df <- data.frame()
  }
  
  get_anova_interactions <- function(cohort_vector){  
    cohort_comb_list <- list()
    for (cohort_ind in 1:length(cohort_vector)){ 
      cohort_comb_df <- as.data.frame(gtools::combinations(length(cohort_vector), cohort_ind, cohort_vector, repeats.allowed = FALSE), stringsAsFactors = FALSE)
      for (comb_index in 1:nrow(cohort_comb_df)){
        cohort_comb_list <- c(cohort_comb_list, list(as.character(cohort_comb_df[comb_index, ])))
      }
    }
    cohort_interactions <- vector()
    for (cohort_comb in cohort_comb_list){
      cohort_comb <- cohort_comb[order(match(cohort_comb, cohort_vector))]
      cohort_interactions <- c(cohort_interactions, paste(cohort_comb, collapse = ":"))
    }
    return (cohort_interactions)  
  }
  
  anova_interactions <- get_anova_interactions(unname(anova_cohort_cols))
  names(anova_interactions) <- get_anova_interactions(names(anova_cohort_cols))
  
  anova_data <- data.frame(stringsAsFactors = FALSE, check.names = FALSE)  
  run_anova <- function(anova_input_df = NULL, anova_cohort_cols = NULL, anova_interactions = NULL){ 
    row_anova <- NULL
    anova_r <- NULL
    tryCatch({  
      anova_input_df <- anova_input_df[apply(anova_input_df, 1, function(x) is.finite(as.numeric(x[['value']]))), , drop = FALSE]                                       
      frm <- paste("value", paste(names(anova_cohort_cols), collapse = " * "), sep = " ~ ")
      anv_lm <- stats::lm(stats::formula(frm), anova_input_df)
      aov_obj <- stats::aov(anv_lm)
      anova_r <- summary(aov_obj)[[1]]
      anova_r <- data.frame(interaction = stringr::str_trim(row.names(anova_r)), anova_r[, c("F value", "Pr(>F)")], stringsAsFactors = FALSE, check.names = FALSE)
      anova_r <- anova_r[!stringr::str_trim(row.names(anova_r)) %in% c("Residuals"), , drop = FALSE]
      colnames(anova_r) <- c("interaction", "F.Value", "P.Value")
      row.names(anova_r) <- NULL
      
      diff_interaction <- base::setdiff(names(anova_interactions), anova_r$interaction)                                                
      if (length(diff_interaction) > 0){
        interm_anova_df <- data.frame(interaction = diff_interaction, F.Value = NA, P.Value = NA, stringsAsFactors = FALSE, check.names = FALSE)
        anova_r <- rbind(anova_r, interm_anova_df)
      }
      row_anova <- anova_r                                     
    }, 
    error = function(cond) {message(paste("\nCannot run anova, caused an error: ", cond))}
    )
    
    if (identical(row_anova, NULL) | !(all(c("interaction", "F.Value", "P.Value") %in% colnames(row_anova)))){
      row_anova <- data.frame(interaction = names(anova_interactions), F.Value = NA, P.Value = NA, stringsAsFactors = FALSE, check.names = FALSE)
    }
    
    return (row_anova)                                   
  }
  
  convert_wide_to_long <- function(data_mat = NULL, metadata = NULL){                                          
    long_data_mat <- data_mat[, metadata[, 1], drop = FALSE]     
    long_data_mat$id <- row.names(long_data_mat)
    long_data_mat <- reshape2::melt(long_data_mat, id.vars = "id")
    long_data_mat <- base::merge(long_data_mat, metadata, by.x = "variable", by.y = 1)
    
    return (long_data_mat)  
  }
  
  if (identical(expr_stat, TRUE)){
    metadata <- suppressMessages(PollyCommonR::create_multi_cohorts_metadata_for_anova(metadata, cohort_col = anova_cohort_cols))
  }
  
  anova_data <- NULL
  tryCatch({                                     
    long_intensity_data <- convert_wide_to_long(data_mat = norm_data, metadata = metadata)                                                                           
    summarise_vars <- c(anova_cohort_cols, "value", "variable")
    anova_res <- long_intensity_data %>% dplyr::group_by_at(c("id")) %>% dplyr::summarise(run_anova(anova_input_df = across(all_of(summarise_vars)), anova_cohort_cols = anova_cohort_cols, anova_interactions = anova_interactions))
    anova_res <- as.data.frame(anova_res)
    anova_res <- anova_res[order(match(anova_res$id, row.names(norm_data))), , drop = FALSE]
    anova_data <- anova_res    
  },
  error = function(cond) {message(paste("\nCannot run anova on this dataset, caused an error: ", cond))}    
  )
  
  tryCatch({                                   
    if ((nrow(anova_data) > 0) && (nrow(identifier_df) > 0)){ 
      anova_data <- merge(identifier_df, anova_data, by = 'id', sort = FALSE)
    }
  }, 
  error = function(cond) {message(paste("\nCannot merge anova data and identifier data, caused an error: ", cond))}
  )                                       
  
  if (identical(expr_stat, TRUE)){
    expr_stat_df <- NULL
    tryCatch({  
      if (identical(expr_with_raw, TRUE)){
        long_intensity_data <- convert_wide_to_long(data_mat = raw_data, metadata = metadata)
      }
      
      expr_stat_res <- data.frame(stringsAsFactors = FALSE, check.names = FALSE)
      for (interaction_type in anova_interactions){                               
        sample_mat_mean_df <- long_intensity_data %>% dplyr::group_by_at(c("id", interaction_type)) %>% dplyr::summarise(Mean = mean(value, na.rm = TRUE))
        average_expr_df <- long_intensity_data %>% dplyr::group_by_at("id") %>% dplyr::summarise(AveExpr = mean(value, na.rm = TRUE))
        cohort_wise_expr_df <- sample_mat_mean_df %>% dplyr::group_by_at("id") %>% dplyr::summarise(MaxExpr = max(Mean, na.rm = TRUE), MinExpr = min(Mean, na.rm = TRUE))                      
        interm_stat_df <- base::merge(average_expr_df, cohort_wise_expr_df, by= "id", sort = FALSE)
        interm_stat_df$interaction <- interaction_type  
        expr_stat_res <- dplyr::bind_rows(expr_stat_res, interm_stat_df)    
      }
      
      expr_stat_df <- expr_stat_res  
    },
    error = function(cond) {message(paste("\nCannot run anova on this dataset, caused an error: ", cond))}    
    )
    
    if (!identical(anova_data, NULL) && !identical(expr_stat_df, NULL)){
      tryCatch({
        anova_data <- merge(anova_data, expr_stat_df, id = c("id", "interaction"), sort = FALSE)
      },
      error = function(cond) {message(paste("\nCannot merge anova data and expression stat dataframes, caused an error: ", cond))}    
      )
    }      
  }
  
  message("Compute Anova Completed...")
  
  return(anova_data)
}