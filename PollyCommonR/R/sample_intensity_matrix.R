#' sample_intensity_matrix
#'
#' performs pca on the sample raw matrix
#'
#' @param raw_intensity_df A dataframe containing samples with raw values.
#' @param metadata_df A dataframe with samples to cohort mapping.
#' @param rownames_col A raw_intensity_df column which is to be used for assigning rownames.
#' @return A marix with samples in columns and ids in rows with raw values.
#' @examples 
#' sample_intensity_matrix(raw_intensity_df, metadata_df, rownames_col = 'Id')
#' @export
sample_intensity_matrix <- function(raw_intensity_df = NULL, metadata_df = NULL, rownames_col = NULL){
  message("Make Sample Intensity Matrix Started...")
  
  metadata_sample <- metadata_df[,1]
  raw_intensity_cols <- colnames(raw_intensity_df)
  sample_cols <- intersect(metadata_sample, raw_intensity_cols)
  if (length(sample_cols) == 0){
    message("No common samples found, please provide valid data")
    return(NULL)
  }
  sample_intensity_mat <- raw_intensity_df[,sample_cols]
  
  if (identical(rownames_col, NULL)==FALSE){
    if (rownames_col %in% colnames(raw_intensity_df)){
      make_unique_names <- make.unique(as.character(raw_intensity_df[,rownames_col]))
      rownames(sample_intensity_mat) <- make_unique_names
    } else{
      message(c(rownames_col, " is not present in raw_intensity_df columns, returning sample matrix with default rownames"))
    }
  } else{
    message("no rownames_col parameter specified, returning sample matrix with default rownames")
  }
  
  sample_intensity_mat <- sample_intensity_mat[rowSums(is.na(sample_intensity_mat)) == 0, ]
  
  message("Make Sample Intensity Matrix Completed...")
  
  return(sample_intensity_mat)
}

