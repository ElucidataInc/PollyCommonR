#' create_cohortwise_cov_boxplot
#'
#' Makes boxplot for cov values across cohorts
#'
#' @param calculated_cov_df A dataframe with calculated cov
#' @param cohorts_order The order of cohorts
#' @param interactive make plot interactive if set TRUE
#' @return A plotly object
#' @examples
#' create_cohortwise_cov_boxplot(calculated_cov_df, cohorts_order = c('CohortA','CohortB'))
#' @import dplyr stringr ggplot2 plotly
#' @export
create_cohortwise_cov_boxplot <- function(calculated_cov_df, cohorts_order = NULL, interactive = FALSE){
  
  message("Create Coefficient of Variation Boxplot Started...")
  require(dplyr)
  require(stringr)
  require(ggplot2)
  require(plotly)
  
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
  calculated_cov_df_filtered <-  calculated_cov_df %>% dplyr::filter(cohort %in% filtered_cohorts_vec)
  all_cv <- calculated_cov_df_filtered$cv
  all_cv <- all_cv[!is.na(all_cv) & !is.infinite(all_cv)]
  if (interactive == TRUE){
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
                                                   list("zoomOut2d")),
                             mathjax = 'cdn')
  } else{
    p <- ggplot(calculated_cov_df_filtered, aes(x = cohort, y = cv, fill=cohort))+
      geom_boxplot(show.legend = FALSE)+
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