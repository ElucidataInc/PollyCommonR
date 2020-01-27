#' create_cohortwise_cov_boxplot
#'
#' Makes boxplot for cov values across cohorts
#'
#' @param calculated_cov_df A dataframe with calculated cov
#' @param cohorts_order The order of cohorts
#' @return A plotly object
#' @examples
#' create_cohortwise_cov_boxplot(calculated_cov_df, cohorts_order = c('CohortA','CohortB'))
#' @import dplyr stringr plotly
#' @export
create_cohortwise_cov_boxplot <- function(calculated_cov_df, cohorts_order = NULL){
  
  message("Create Coefficient of Variation Boxplot Started...")
  require(dplyr)
  require(stringr)
  require(plotly)
  
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
  calculated_cov_df_filtered <-  calculated_cov_df %>% dplyr::filter(cohort %in% filtered_cohorts_vec)
  
  p <- plot_ly(calculated_cov_df_filtered, y = ~cv, color = ~cohort, type = "box") %>%
    layout(
      title = "CV Distribution across Cohorts",
      xaxis = list(
        title = "Cohorts",
        titlefont = list(size=16),
        categoryorder = "array",
        categoryarray = filtered_cohorts_vec
      ),
      yaxis = list(title = "Coefficient of Variation (%)", titlefont = list(size=16)),
      showlegend = FALSE
    ) %>% plotly::config(displaylogo = FALSE,
                         modeBarButtons = list(list("zoomIn2d"), 
                                               list("zoomOut2d"), 
                                               list('toImage')),
                         mathjax = 'cdn')
  
  message("Create Coefficient of Variation Boxplot Completed...")
  
  return (p)
  
}