#' compute_posthoc
#'
#' Performs pairwise t-tests and adjusts p-values using specified methods.
#'
#' @param norm_data The dataframe/matrix with samples in columns and features in rows.
#' @param metadata The dataframe with samples to cohort mapping.
#' @param cohort_col A vector of metadata columns used for grouping.
#' @param post_hoc_test The method used for adjusting p-values, options include "bonferroni", "fdr" (False Discovery Rate), and "holm".
#' @return A list containing:
#'   - posthoc_data: A dataframe with adjusted p-values for each pairwise comparison.
#'   - posthoc_data_long_df: A long-format dataframe containing detailed results of pairwise comparisons.
#'   - posthoc_data_heatmap: A custom wide-format dataframe suitable for visualization.
#' @examples 
#' compute_posthoc(norm_data, metadata, cohort_col, post_hoc_test = "bonferroni")
#' @import dplyr stats stringr tidyr
#' @export
compute_posthoc <- function(norm_data = NULL, metadata = NULL, cohort_col = "Cohort", post_hoc_test = "bonferroni"){
  message("Compute Post hoc Started...")
  require(dplyr)
  require(stats)
  require(stringr)
  require(tidyr)
  
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
  
  convert_wide_to_long <- function(data_mat = NULL, metadata = NULL){                                          
    long_data_mat <- data_mat[, metadata[, 1], drop = FALSE]     
    long_data_mat$id <- row.names(long_data_mat)
    long_data_mat <- reshape2::melt(long_data_mat, id.vars = "id")
    long_data_mat <- base::merge(long_data_mat, metadata, by.x = "variable", by.y = 1)
    
    return (long_data_mat)  
  }
  
  message("Converting data to long format...")
  long_intensity_data <- tryCatch({
    convert_wide_to_long(data_mat = norm_data, metadata = metadata)
  }, error = function(e) {
    message("Error occurred during data conversion to long format:", e)
    long_intensity_data <- NULL
    showNotification("Failed to convert data to long format. Please check the data and try again.", duration = 10, type = 'error', closeButton = TRUE)
    return(NULL)
  })
  message("Data conversion to long format completed.")
  
  subset_df <- long_intensity_data
  # subset_df <- subset_df[c("id", anova_cohort_cols, "value", "variable")]
  
  perform_pairwise_t_tests <- function(input_df, anova_cohort_cols, p_adjust_method = "bonferroni") {
    # Initialize an empty list to store the results
    results_list <- list()
    
    # Loop over each unique metabolite id
    for (id in unique(input_df$id)) {
      # Subset the input dataframe for the current metabolite id
      subset_df <- input_df[input_df$id == id, ]
      
      # Perform pairwise t-tests for each specified grouping variable
      for (col in anova_cohort_cols) {
        # Perform pairwise t-test
        pairwise_results <- pairwise.t.test(subset_df$value, subset_df[[col]],
                                            p.adjust.method = p_adjust_method)
        
        # Check if all values in pairwise results are NA or "-"
        if(all(is.na(pairwise_results$p.value) | pairwise_results$p.value == "-")) {
          # If all values are NA or "-", skip this iteration
          next
        }
        
        # Convert pairwise results to dataframe
        results_df <- as.data.frame(pairwise_results$p.value)
        
        # Remove comparisons between identical groups
        results_df[rownames(results_df) == col, ] <- NA
        
        # Convert row names to a new column
        results_df$Comparison_1 <- rownames(results_df)
        
        # Reshape the dataframe to long format
        long_results_df <- reshape2::melt(results_df, id.vars = "Comparison_1", 
                                          variable.name = "Comparison_2", 
                                          value.name = "adj.p")
        
        # Remove rows with NA values
        long_results_df <- long_results_df[!is.na(long_results_df$adj.p), ]
        
        # Combine Comparison_1 and Comparison_2 columns into a new column
        long_results_df$Comparison <- paste(long_results_df$Comparison_1, long_results_df$Comparison_2, sep = ":")
        
        # Remove the now redundant columns
        long_results_df <- long_results_df[c("Comparison", "adj.p")]
        
        # Add metabolite id and grouping variable name to the results
        long_results_df$ID <- id
        long_results_df$Interaction <- col
        long_results_df$P_Adjust_Method <- p_adjust_method
        
        # Reorder the columns as per the requirement
        long_results_df <- long_results_df[, c("ID", "Interaction", "Comparison", "adj.p", "P_Adjust_Method")]
        
        # Save the results in the list
        results_list[[length(results_list) + 1]] <- long_results_df
      }
    }
    
    # Combine all results into a single dataframe
    final_results_df <- do.call(rbind, results_list)
    
    # Reset rownames
    rownames(final_results_df) <- NULL
    
    # Return the final dataframe
    return(final_results_df)
  }
  
  # Define a function to generate pairwise column combinations
  generate_column_combinations <- function(cols) {
    if (length(cols) == 1) {
      message("Only one column selected.")
      return(cols)
    } else if (length(cols) == 2) {
      combination <- paste(cols[1], cols[2], sep = ":")
      if (combination %in% colnames(long_intensity_data)) {
        message("Pairwise column combination '", combination, "' exists in the data.")
        return(combination)
      } else {
        message("The specified column combination '", combination, "' does not exist in the data.")
        return(NULL)
      }
    } else {
      message("More than two columns selected. Pairwise comparisons will not be performed.")
      return(NULL)
    }
  }
  
  tryCatch({
    anova_cohort_cols <- as.vector(anova_cohort_cols)
    anova_cohort_cols <- unique(c(anova_cohort_cols, generate_column_combinations(anova_cohort_cols)))
  }, error = function(e) {
    message("An error occurred: ", conditionMessage(e))
    # Additional error handling code
    showNotification("Failed to generate column combinations. Please check the data and try again.", duration = 10, type = 'error', closeButton = TRUE)
    return(NULL)
  })
  
  message("Performing pairwise t-tests...")
  pairwise_results <- tryCatch({
    if(post_hoc_test == 'bonferroni'){
      message("Performing Bonferroni adjustment...")
      perform_pairwise_t_tests(long_intensity_data, anova_cohort_cols, p_adjust_method = "bonferroni")
    } else if(post_hoc_test == 'fdr'){
      message("Performing FDR adjustment...")
      perform_pairwise_t_tests(long_intensity_data, anova_cohort_cols, p_adjust_method = "fdr")
    } else if(post_hoc_test == 'holm'){
      message("Performing Holm adjustment...")
      perform_pairwise_t_tests(long_intensity_data, anova_cohort_cols, p_adjust_method = "holm")
    }
  }, error = function(e) {
    message("Error occurred during pairwise t-tests:", e)
    pairwise_results <- NULL
    showNotification("Failed to perform pairwise t-tests. Please check the data and try again.", duration = 10, type = 'error', closeButton = TRUE)
    return(NULL)
  })
  
  if (!is.null(pairwise_results)) {
    message("Pairwise t-tests completed successfully.")
  } else {
    message("Pairwise t-tests failed.")
  }
  
  
  # Convert Bonferroni or FDR or Holm results to wide format
  convert_to_wide_format <- function(df) {
    # Check if the necessary columns exist
    if (!all(c("ID", "Comparison", "adj.p") %in% colnames(df))) {
      stop("Required columns not found in the dataframe.")
    }
    
    # Check if adj.p is numeric
    if (!is.numeric(df$adj.p)) {
      stop("adj.p column should be numeric.")
    }
    
    df_filtered <- df[, c("ID", "Comparison", "adj.p")]

    # Remove duplicate rows
    df_filtered <- unique(df_filtered)

    # Pivot the data to wide format
    wide_df <- tidyr::pivot_wider(df_filtered, names_from = Comparison, values_from = adj.p)
    
    return(wide_df)
  }
  
  convert_to_wide_format_custom <- function(df) {
    # Check if the necessary columns exist
    if (!all(c("ID", "Comparison", "adj.p", "P_Adjust_Method") %in% colnames(df))) {
      stop("Required columns not found in the dataframe.")
    }
    
    # Check if adj.p is numeric
    if (!is.numeric(df$adj.p)) {
      stop("adj.p column should be numeric.")
    }
    
    df_filtered <- df[, c("ID", "Comparison", "adj.p", "P_Adjust_Method")]
    
    # Add adjustment method to the column headers
    df_filtered$Comparison <- paste(df_filtered$Comparison, "adj.p", df_filtered$P_Adjust_Method, sep = ", ")

    # Remove duplicate rows
    df_filtered <- unique(df_filtered)

    # Pivot the data to wide format
    wide_df <- tidyr::pivot_wider(df_filtered, names_from = Comparison, values_from = adj.p)
    
    # Remove the P_Adjust_Method column
    wide_df$P_Adjust_Method <- NULL
    
    return(wide_df)
  }
  
  message("Converting results to wide format...")
  wide_df <- tryCatch({
    convert_to_wide_format(pairwise_results)
  }, error = function(e) {
    message("Error occurred during conversion to wide format:", e)
    wide_df <- NULL
    showNotification("Failed to convert results to wide format. Please check the data and try again.", duration = 10, type = 'error', closeButton = TRUE)
    return(NULL)
  })
  message("Wide format conversion completed.")
  
  message("Converting results to custom wide format...")
  wide_df_custom <- tryCatch({
    convert_to_wide_format_custom(pairwise_results)
  }, error = function(e) {
    message("Error occurred during custom wide format conversion:", e)
    wide_df_custom <- NULL
    showNotification("Failed to convert results to custom wide format. Please check the data and try again.", duration = 10, type = 'error', closeButton = TRUE)
    return(NULL)
  })
  message("Custom wide format conversion completed.")
  
  
  message("Compute Post hoc Completed...")

  return(list(posthoc_data = wide_df, posthoc_data_long_df = pairwise_results, posthoc_data_heatmap = wide_df_custom))
}