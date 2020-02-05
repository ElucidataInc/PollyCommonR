#' create_cohortwise_cov_barplot
#'
#' Makes boxplot for cov values across cohorts
#'
#' @param calculated_cov_df A dataframe with calculated cov
#' @param id_order A vector of ids to plot
#' @param id_col The column name where ids are present
#' @param cohorts_order The order of cohorts
#' @param interactive make plot interactive if set TRUE
#' @return A plotly object
#' @examples
#' create_cohortwise_cov_barplot(calculated_cov_df, id_order = c('cmpd1','cmpd2'), id_col = 'id', cohorts_order = c('CohortA','CohortB'))
#' @import dplyr stringr plotly
#' @export
create_cohortwise_cov_barplot <- function(calculated_cov_df, id_order = NULL, id_col = 'id', 
                                          cohorts_order = NULL, interactive = FALSE){
  
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
  
  calculated_cov_df[[id_col]] <- stringr::str_trim(calculated_cov_df[[id_col]])
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
  
  calculated_cov_df[['cohort']] <- stringr::str_trim(calculated_cov_df[['cohort']])
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
  all_cv <- calculated_cov_df_filtered$cv
  all_cv <- all_cv[!is.na(all_cv) & !is.infinite(all_cv)]
  
  if (interactive == TRUE){
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
  } else {
    p<-ggplot(calculated_cov_df_filtered, aes(x=cohort, y=cv, fill=!!(sym(id_col)))) +
      geom_bar(stat="identity", position = "dodge")+
      ggtitle("CV Distribution across Cohorts")+
      labs(x = "Cohorts",
           y = "Coefficient of Variation (%)")+
      scale_x_discrete(limits = filtered_cohorts_vec, expand = c(0.07,0))+
      scale_y_continuous(breaks=seq(0, max(all_cv), 25))+
      ggsci::scale_color_aaas() + # filling the point colors
      theme(axis.line = element_line(size = 1, colour = "black"), # axis line of size 1 inch in black color
            panel.grid.major = element_blank(), # major grids included
            panel.grid.minor = element_blank(), # no minor grids
            panel.border = element_blank(), panel.background = element_blank(), # no borders and background color
            plot.title = element_text(colour="black", size = 18, face = "plain", hjust=0.5),
            axis.title = element_text(colour="black", size = 14, face = "plain"), # axis title
            axis.text.x = element_text(colour="black", size = 10, angle = 90,
                                       hjust = 1, margin=unit(c(0.5,0.5,0.1,0.1), "cm"),
                                       face = "plain"), # x-axis text in fontsize 10
            axis.text.y = element_text(colour="black", size = 10,
                                       margin=unit(c(0.5,0.5,0.1,0.1), "cm"), 
                                       face = "plain"), # y-axis text in fontsize 10
            axis.ticks.length = unit(-0.25, "cm"))
  }
  
  message("Create Coefficient of Variation Boxplot Completed...")
  
  return (p)
  
}