#' make_intomix_input
#'
#' Create an intomix input file
#'
#' @param norm_data A sample intenisty data with samples in columns.
#' @param metadata dataframe containing metadata information
#' @param cohorts A vector of cohorts used for cohorts comparison
#' @param cohorts_compare_data A dataframe having predefined cohorts comparison
#' @param cohort_col A metadata column where cohorts are present
#' @param rownames_col A raw_intensity_df column which is to be used for assigning rownames.
#' @param remove_duplicates Remove duplicate ids rows based on pval
#' @param id_col A column id which will be used to remove duplicates
#' @param p_val_correct_methods character vector containing p-val correction methods ('Bonferroni', 'BH(FDR)') to be applied.
#' @param log_flag Check if data is log transformed (TRUE) or not (FALSE). If not (FALSE) then do internally log2 transformation.
#'
#' @return A dataframe of differential expression for all cohorts comparison.
#' @examples
#' make_intomix_input(norm_data = NULL, metadata = NULL, cohorts = NULL,  
#'                    cohorts_compare_data = NULL, cohort_col = 'Cohort', 
#'                    rownames_col = 'uniqueId', remove_duplicates = TRUE, 
#'                    id_col = 'compoundId', p_val_correct_methods = 'BH',
#'                    log_flag = TRUE )
#' @import stringr
#' @export
make_intomix_input <- function(norm_data = NULL, metadata = NULL,
                               cohorts = NULL,  cohorts_compare_data = NULL, cohort_col = 'Cohort', 
                               rownames_col = 'uniqueId', remove_duplicates = TRUE, 
                               id_col = 'compoundId', p_val_correct_methods = 'BH',
                               log_flag = TRUE ){
  message("Make Intomix Input Started...")
  
  if (identical(rownames_col, NULL)){
    warning("no 'rownames_col' parameter specified which is a norm_data column")
    return(NULL)
  }
  if (!(rownames_col %in% colnames(norm_data))){
    warning(c(rownames_col, " is not present in norm_data columns"))
    return(NULL)
  }
  
  if (identical(id_col, NULL)){
    warning("no 'id_col' parameter specified which is a norm_data column")
    return(NULL)
  }
  if (!(id_col %in% colnames(norm_data))){
    warning(c(id_col, " is not present in norm_data columns"))
    return(NULL)
  }
  
  if (identical(cohorts_compare_data, NULL)){  
    
    if (identical(cohorts, NULL)){
      warning("No 'cohorts' parameter specified")
      return(NULL)
    }
    
    metadata[[cohort_col]] <- stringr::str_trim(metadata[[cohort_col]])
    cohorts <- stringr::str_trim(cohorts)
    common_cohorts <- intersect(unique(metadata[[cohort_col]]), cohorts)
    
    if (length(common_cohorts) < 2) {
      warning("Please provide more than two common cohorts present in cohorts vector and metadata")
      return(NULL)
    }
    
    cohorts_comb_df <- expand.grid(state1 = common_cohorts, state2 = common_cohorts, stringsAsFactors = F)
    cohorts_comb_df <- cohorts_comb_df[apply(cohorts_comb_df,1 , function(x) x["state2"] != x["state1"]),]
    cohorts_comb_df <- cohorts_comb_df[order(cohorts_comb_df$state1), ]
    
  } else {
    if (!identical(class(cohorts_compare_data), "data.frame")){
      warning("Please provide a valid dataframe. The 'cohorts_compare_data' is a dataframe with 'state1' and 'state2' columns")
      return(NULL) 
    } else {
      if (!all((c("state1", "state2") %in% colnames(cohorts_compare_data)))){
        warning("The 'cohorts_compare_data' dataframe should have 'state1' and 'state2' columns")
        return(NULL)  
      } else {
        cohorts_comb_df <- cohorts_compare_data
      }
    }   
  }   
  
  metadata[[cohort_col]] <- stringr::str_trim(metadata[[cohort_col]])
  norm_data$uniqueId <- make.unique(stringr::str_trim(as.character(norm_data[, rownames_col])))
  sample_raw_mat <- PollyCommonR::sample_intensity_matrix(norm_data, metadata, rownames_col = "uniqueId")
  
  overall_diff_exp <- data.frame()
  for (each_comb in 1:nrow(cohorts_comb_df)){
    cohort_a <- stringr::str_trim(as.character(cohorts_comb_df[each_comb,]$state2))
    cohort_b <- stringr::str_trim(as.character(cohorts_comb_df[each_comb,]$state1)) 
    if (cohort_a == cohort_b){
      next
    }
    
    compare_cohorts <- c(cohort_a, cohort_b)
    
    if (!all(compare_cohorts %in% metadata[[cohort_col]])){
      next
    }
    
    diff_exp <- NULL
    try(diff_exp <- PollyCommonR::diff_exp_limma(prot_norm_mat = sample_raw_mat, 
                                                 metadata = metadata, cohort_col,
                                                 compare_cohorts[1], compare_cohorts[2],
                                                 p_val_correct_methods, log_flag), silent = TRUE)
    
    if (identical(diff_exp, NULL)){
      next
    }
    
    diff_exp_update <- data.frame(id = rownames(diff_exp), state1 = compare_cohorts[2], 
                                  state2 = compare_cohorts[1], diff_exp, stringsAsFactors = FALSE)
    merged_diff_exp <- merge(norm_data[, c("uniqueId", id_col)], diff_exp_update, by.x = "uniqueId", by.y = "id")
    
    if (remove_duplicates){
      merged_diff_exp_sorted <- merged_diff_exp[order(merged_diff_exp$P.Value, decreasing = FALSE), ]
      merged_diff_exp_filtered <- merged_diff_exp_sorted[!(duplicated(merged_diff_exp_sorted[[id_col]])), ]
      overall_diff_exp <- rbind(overall_diff_exp, merged_diff_exp_filtered)
    } else {
      overall_diff_exp <- rbind(overall_diff_exp, merged_diff_exp) 
    }    
    
  }
  
  if (nrow(overall_diff_exp) < 1){
    warning("The 'cohorts_compare_data' does not contain valid cohorts for comparison")
    return(NULL) 
  }
  
  overall_diff_exp$ID <- overall_diff_exp[[id_col]]
  intomix_input_df <- overall_diff_exp[,c('uniqueId', 'ID','state1','state2','P.Value','logFC')]
  intomix_input_df_filtered <- intomix_input_df
  colnames(intomix_input_df_filtered) <- c('uniqueId','ID','state1','state2','pval','log2FC')
  
  message("Make Intomix Input Started...")
  
  return (intomix_input_df_filtered)
}