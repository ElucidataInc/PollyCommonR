#' make_log2fc_matrix_from_intomix_input
#'
#' It creates the log2FC matrix from the intomix input format
#'
#' @param intomix_input The intomix input
#' @param selected_comparisons The vector of selected comparisons in "state1 vs state2" format 
#' @param pval_cutoff The pval cutoff to filter the data
#' @return A dataframe of comparisons in columns and ids in rows
#' @examples
#' make_log2fc_matrix_from_intomix_input(intomix_input, selected_comparisons, pval_cutoff = 0.05)
#' @import dplyr
#' @export
make_log2fc_matrix_from_intomix_input <- function(intomix_input = NULL, selected_comparisons = NULL, pval_cutoff = 0.05){
  message("Make Log2FC Matrix From Intomix Input Started...")
  
  if (identical(intomix_input, NULL)){
    warning("The intomix_input is NULL")
    return (NULL)  
  }
  
  if (!identical(as.character(class(intomix_input)), "data.frame")){
    warning(" The intomix_input is not a dataframe")
    return (NULL) 
  }
  
  required_cols <- c("uniqueId", "ID", "state1", "state2", "pval", "log2FC")
  if (!all(required_cols %in% colnames(intomix_input))){
    warning(c("The intomix_input dataframe should have the following columns: ", paste0(required_cols, collapse = ", ")))
    return (NULL)      
  }
  
  if (identical(pval_cutoff, NULL)){
    warning("The pval_cutoff is NULL, using 0.05 as default value.")
    pval_cutoff <- 0.05
  }
  
  intomix_input$comparison <- paste(intomix_input$state1, intomix_input$state2, sep = " vs ")
  total_comparisons <- unique(intomix_input$comparison)
  common_comparisons <- total_comparisons
  
  if (length(selected_comparisons) != 0) {
    common_comparisons <- base::intersect(selected_comparisons, total_comparisons)
  }
  
  if (length(common_comparisons) == 0) {
    warning("Invalid selected comparisons, the comparisons should be in 'state1 vs state2' format.")
    return(NULL)
  }  
  
  comparisonwise_data_list <- list()
  for ( single_comp in common_comparisons){
    single_comp_data <- dplyr::filter(intomix_input, comparison == single_comp, pval <= pval_cutoff)
    if (nrow(single_comp_data) >= 1){  
      single_comp_data[[single_comp]] <- single_comp_data$log2FC
      comparisonwise_data_list[[single_comp]] <- single_comp_data[, c("uniqueId", single_comp), drop = FALSE]
    }
  }
  
  log2fc_matrix_df <-  Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c("uniqueId"), all = TRUE), comparisonwise_data_list)
  
  if (identical(log2fc_matrix_df, NULL)){
    warning("Unable to create log2fc matrix, please change parameters")
    return(NULL)      
  }
  
  row.names(log2fc_matrix_df) <- log2fc_matrix_df$uniqueId
  log2fc_matrix_df$uniqueId <- NULL
  
  message("Make Log2FC Matrix From Intomix Input Completed...")
  
  return (log2fc_matrix_df)                   
}