#' create_cohortwise_cov_barplot
#'
#' Makes boxplot for cov values across cohorts
#'
#' @param calculated_cov_df A dataframe with calculated cov
#' @param cohorts_order The order of cohorts
#' @return A plotly object
#' @examples
#' create_cohortwise_cov_barplot(calculated_cov_df, id_order = c('cmpd1','cmpd2'), id_col = 'id', cohorts_order = c('CohortA','CohortB'))
#' @import dplyr stringr plotly
#' @export
create_cohortwise_cov_barplot <- function(calculated_cov_df, id_order = NULL, id_col = 'id', cohorts_order = NULL){
  
  message("Create Coefficient of Variation Boxplot Started...")
  require(dplyr)
  require(stringr)
  require(plotly)
  
  if (identical(id_order, NULL)){
    warning("Please provide valid ids")
    return(NULL)
  }
  
  if (identical(id_col, NULL)){
    warning("Please provide valid id column name")
    return(NULL)
  }
  
  if (!(id_col %in% colnames(calculated_cov_df))){
    warning(c(id_col, " column is not present in the input dataframe"))
    return(NULL)
  }
  
  calculate_cov_df[[id_col]] <- stringr::str_trim(calculate_cov_df[[id_col]])
  ids_vec <- unique(calculated_cov_df[[id_col]])
  common_ids <- c()
  if (length(id_order) != 0){
    id_order <- stringr::str_trim(id_order)
    common_ids <- intersect(id_order, ids_vec)
  }
  if (length(common_ids) == 0){
    warning("Invalid selected ids")
    return(NULL)
  }
  diff_ids <- setdiff(id_order, common_ids)
  if (length(diff_ids) != 0){
    message(c("The following are not valid ids : ", diff_ids))
  }
  
  calculate_cov_df[['cohort']] <- stringr::str_trim(calculate_cov_df[['cohort']])
  cohorts_vec <- unique(calculated_cov_df[['cohort']])
  filtered_cohorts_vec <- cohorts_vec
  common_cohorts <- c()
  if (length(cohorts_order) != 0){
    cohorts_order <- stringr::str_trim(cohorts_order)
    common_cohorts <- intersect(cohorts_order, cohorts_vec)
  }
  if (length(common_cohorts) != 0){
    filtered_cohorts_vec = common_cohorts
  }
  calculated_cov_df_filtered <-  calculated_cov_df %>% dplyr::filter(!!(sym(id_col)) %in% common_ids,  cohort %in% filtered_cohorts_vec)
  p <- plot_ly()
  for(id_val in common_ids){
    cmpd_df <- calculated_cov_df_filtered %>% dplyr::filter(!!(sym(id_col)) %in% id_val)
    p <- add_trace(p, x = cmpd_df[,'cohort'], y = cmpd_df[,'cv'], type = 'bar', name = id_val)
  }
  p <- p %>% layout(
    title = "CV Distribution across Cohorts",
    xaxis = list(
      title = "Cohorts",
      titlefont = list(size=16),
      categoryorder = "array",
      categoryarray = filtered_cohorts_vec
    ),
    yaxis = list(title = "Coefficient of Variation (%)", titlefont = list(size=16)),
    showlegend = TRUE
  ) %>% plotly::config(displaylogo = FALSE,
                       modeBarButtons = list(list("zoomIn2d"), 
                                             list("zoomOut2d"), 
                                             list('toImage')),
                       mathjax = 'cdn')
  
  message("Create Coefficient of Variation Boxplot Completed...")
  
  return (p)
  
}