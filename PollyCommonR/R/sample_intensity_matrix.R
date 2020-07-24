#' sample_intensity_matrix
#'
#' create a sample matrix with samples in columns and ids in rownames.
#'
#' @param raw_intensity_df A dataframe containing samples with raw values.
#' @param metadata_df A dataframe with samples to cohort mapping.
#' @param rownames_col A raw_intensity_df column which is to be used for assigning rownames.
#' @param drop_na Drop rows having NA's.
#' @param replace_na_with_zero Replace all NA's with zero.
#' @return A marix with samples in columns and ids in rows with raw values.
#' @examples 
#' sample_intensity_matrix(raw_intensity_df, metadata_df, rownames_col = 'Id')
#' @import stringr
#' @export
sample_intensity_matrix <- function(raw_intensity_df = NULL, metadata_df = NULL, 
                                    rownames_col = NULL, drop_na = FALSE, replace_na_with_zero = FALSE){
  message("Make Sample Intensity Matrix Started...")
  
  metadata_sample <- metadata_df[, 1]
  raw_intensity_cols <- colnames(raw_intensity_df)
  sample_cols <- intersect(metadata_sample, raw_intensity_cols)
  
  if (length(sample_cols) == 0){
    message("No common samples found, please provide valid data")
    return(NULL)
  }
  sample_intensity_mat <- raw_intensity_df[ , sample_cols, drop = FALSE]
  
  if (identical(rownames_col, NULL)==FALSE){
    if (rownames_col %in% colnames(raw_intensity_df)){
      make_unique_names <- make.unique(stringr::str_trim(as.character(raw_intensity_df[ ,rownames_col])))
      rownames(sample_intensity_mat) <- make_unique_names
    } else{
      message(c(rownames_col, " is not present in raw_intensity_df columns, returning sample matrix with default rownames"))
    }
  } else{
    message("no rownames_col parameter specified, returning sample matrix with default rownames")
  }
  
  sample_intensity_mat[, colnames(sample_intensity_mat)] <- apply(sample_intensity_mat[, colnames(sample_intensity_mat), drop = FALSE], 2, function(x) as.numeric(x))
  
  if (drop_na){
    sample_intensity_mat <- sample_intensity_mat[rowSums(is.na(sample_intensity_mat)) == 0, , drop = FALSE]
  }
  
  if (replace_na_with_zero){
    sample_intensity_mat[is.na(sample_intensity_mat)] <- 0
  }
  
  message("Make Sample Intensity Matrix Completed...")
  
  return(sample_intensity_mat)
}