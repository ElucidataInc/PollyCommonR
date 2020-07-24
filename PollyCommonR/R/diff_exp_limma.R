#' diff_exp_limma
#'
#' calls limma differential expression algorithms
#'
#' @param prot_norm_mat sample matrix with samples in columns and ids as rownames.
#' @param metadata dataframe containing metadata information
#' @param cohort_condition A metadata column where cohorts are present
#' @param cohort_a string cohort_a
#' @param cohort_b string cohort_b
#' @param p_val_correct_methods character vector containing p-val correction methods ('Bonferroni', 'BH(FDR)') to be applied.
#'
#' @param log_flag Check if data is log transformed (TRUE) or not (FALSE). If not (FALSE) then do internally log2 transformation.
#'
#' @return dataframe, If logFC >0, it implies abundance is greater in cohort_b.
#' @examples
#' diff_exp_limma(prot_norm_mat, metadata, 'Cohort', 'Cohort1', 'Cohort2', p_val_correct_methods = 'BH(FDR)',log_flag = FALSE)
#' @import limma dplyr
#' @export
diff_exp_limma <- function (prot_norm_mat, metadata, cohort_condition, cohort_a, 
                            cohort_b, p_val_correct_methods = "BH", log_flag = TRUE) {
  message("Calculate Differential Expression Limma Started...")
  require(limma)
  require(dplyr)
  
  if (!(cohort_a %in% metadata[[cohort_condition]] && cohort_b %in% metadata[[cohort_condition]])) {
    warning(paste0("Invalid ", cohort_a, " or ", cohort_b))
    return(NULL)
  }
  cohort_a <- make.names(gsub("condition", "", cohort_a))
  cohort_b <- make.names(gsub("condition", "", cohort_b))
  metadata[[cohort_condition]] <- make.names(gsub("condition", "", metadata[[cohort_condition]]))
  metadata[, 1] <- make.names(metadata[, 1])
  colnames(prot_norm_mat) <- make.names(colnames(prot_norm_mat))
  if ((cohort_a == "") | (cohort_b == "")) {
    warning("One of your cohorts is empty, try changing the cohort condition.")
    return(NULL)
  }
  prot_norm_mat_samples <- colnames(prot_norm_mat)
  metadata_samples <- as.character(metadata[, 1])
  common_samples <- intersect(metadata_samples, prot_norm_mat_samples)
  if (length(common_samples) == 0) {
    warning("No common samples in matrix and metadata")
    return(NULL)
  }
  metadata <- dplyr::filter(metadata, !!(sym(colnames(metadata)[1])) %in% !!common_samples)
  metadata <- dplyr::filter(metadata, !!(sym(cohort_condition)) %in% c(cohort_a, cohort_b))
  if (log_flag) {
    prot_norm_mat_log2 <- prot_norm_mat[, metadata[, 1]]
  }
  else {
    prot_norm_mat_log2 <- log2(prot_norm_mat[, metadata[, 1]])
  }
  condition <- metadata[[cohort_condition]]
  design <- model.matrix(~condition + 0)
  colnames(design) <- gsub("condition", "", colnames(design))
  contrast_matrix <- limma::makeContrasts(contrasts = c(paste(cohort_b, cohort_a, sep = "-")), levels = design)
  cohort_a_columns <- metadata[metadata[cohort_condition] == cohort_a, ][, 1]
  cohort_b_columns <- metadata[metadata[cohort_condition] == cohort_b, ][, 1]
  if (length(cohort_a_columns) <= 1 | length(cohort_b_columns) <= 1) {
    warning("Since, your selected cohorts have no replicates in the data, you might want to change the cohorts or cohort condition and try again.")
    return(NULL)
  }
  else {
    fit <- limma::lmFit(prot_norm_mat_log2, design)
    fit <- limma::contrasts.fit(fit, contrast_matrix)
    fit <- limma::eBayes(fit)
    limma_results_df <- limma::topTable(fit, coef = paste(cohort_b, cohort_a, sep = "-"), number = nrow(prot_norm_mat_log2))
    for (p_val_correct in p_val_correct_methods) {
      if (p_val_correct == "Bonferroni") {
        p_val_correct_method <- "bonferroni"
      }
      else {
        p_val_correct_method <- "BH"
      }
      limma_results_df[, p_val_correct] <- p.adjust(limma_results_df$P.Value, 
                                                    method = p_val_correct_method)
    }

    limma_results_df <- limma_results_df[order(match(rownames(limma_results_df), rownames(prot_norm_mat_log2))), , drop = FALSE]

    message("Calculate Differential Expression Limma Completed...")
    return(limma_results_df)
  }
}