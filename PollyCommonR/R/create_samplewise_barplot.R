#' create_samplewise_barplot
#'
#' Makes samplewise barplot
#'
#' @param sample_raw_mat sample_raw_mat matrix/dataframe containing raw values
#' @param metadata_df A dataframe with samples to cohort mapping
#' @param id_name An id for which barplot to be made
#' @param cohorts_order The order of cohorts
#' @param cohort_col A cohort column present in metadata
#' @param interactive Make plot interactive using plotly
#' @param x_label Label x-axis
#' @param y_label Label y-axis
#' @param title_label Title of the plot
#' @return ggplot object or plotly object
#' @examples
#' create_samplewise_barplot(sample_raw_mat = NULL, metadata_df = NULL, id_name = NULL, 
#'                           cohorts_order = NULL, cohort_col='Cohort', 
#'                           x_label = "", y_label = "", title_label = "")
#' @import dplyr stringr ggplot2 plotly
#' @export
create_samplewise_barplot <- function(sample_raw_mat = NULL, metadata_df = NULL,
                                      id_name = NULL, cohorts_order = NULL, 
                                      cohort_col='Cohort', interactive = FALSE, 
                                      x_label = "Sample", y_label = "Intensity", title_label = ""){
  message("Create Samplewise Barplot Started...")
  require(dplyr)
  require(stringr)
  require(ggplot2)
  require(plotly)  
  
  ids_vec <- unique(rownames(sample_raw_mat))
  if (!(id_name %in% ids_vec)){
    warning(c(id_name, " is not a valid id_name, please choose from sample_raw_mat rownames"))
    return (NULL)
  }
  
  if (!(cohort_col %in% colnames(metadata_df))){
    warning(c(cohort_col, " is not a valid cohort_col, please choose from metadata_df colnames"))
    return (NULL)    
  }
  
  metadata_df[[cohort_col]] <- stringr::str_trim(metadata_df[[cohort_col]])
  cohorts_vec <- unique(metadata_df[[cohort_col]])
  filtered_cohorts_vec <- cohorts_vec
  common_cohorts <- c()
  if (length(cohorts_order) != 0){
    cohorts_order <- stringr::str_trim(cohorts_order)
    common_cohorts <- base::intersect(cohorts_order, cohorts_vec)
  }
  if (length(common_cohorts) != 0){
    filtered_cohorts_vec = common_cohorts
  }
  transposed_mat <- as.data.frame(t(sample_raw_mat))
  transposed_mat$Sample <- rownames(transposed_mat)
  metadata_df$Sample <- metadata_df[,1]
  mat_with_metadata <- merge(metadata_df, transposed_mat, by="Sample")
  filtered_mat <- data.frame()
  for (cohort in filtered_cohorts_vec){
    interm_mat <- mat_with_metadata %>% dplyr::filter(!!(sym(cohort_col)) %in% cohort)
    interm_mat <- dplyr::arrange(interm_mat, Sample)
    filtered_mat <- rbind(filtered_mat, interm_mat)  
  }    
  
  if (interactive == TRUE){
    p <- plot_ly(x = filtered_mat[["Sample"]], y = filtered_mat[[id_name]], 
                 color = filtered_mat[[cohort_col]], type = "bar") %>%
      layout(
        title = title_label,
        xaxis = list(
          title = x_label,
          titlefont = list(size=16),
          categoryorder = "array",
          categoryarray = filtered_mat[["Sample"]]
        ),
        yaxis = list(title = y_label, titlefont = list(size=16)),
        showlegend = TRUE
      ) %>% plotly::config(displaylogo = FALSE,
                           modeBarButtons = list(list("zoomIn2d"), 
                                                 list("zoomOut2d")),
                           mathjax = 'cdn')   
  } else{
    p <- ggplot(filtered_mat, aes(x=Sample, y=!!(sym(id_name)), fill=!!(sym(cohort_col))))+
      geom_bar(stat="identity", position = "dodge")+
      ggtitle(title_label)+
      labs(x = x_label, y = y_label)+
      scale_x_discrete(limits = filtered_mat$Sample, expand = c(0.09,0.09))+
      ggsci::scale_color_aaas() + # filling the point colors
      theme(axis.line = element_line(size = 1, colour = "black"), # axis line of size 1 inch in black color
            panel.grid.major = element_blank(), # major grids included
            panel.grid.minor = element_blank(), # no minor grids
            panel.border = element_blank(), panel.background = element_blank(), # no borders and background color
            plot.title = element_text(colour="black", size = 18, face = "plain", hjust=0.5),
            axis.title = element_text(colour="black", size = 14, face = "plain"), # axis title
            axis.text.x = element_text(colour="black", size = 10, angle = 90,
                                       margin=unit(c(0.2,0.2,0.1,0.1), "cm"),
                                       face = "plain"), # x-axis text in fontsize 10
            axis.text.y = element_text(colour="black", size = 10,
                                       margin=unit(c(0.2,0.2,0.1,0.1), "cm"), 
                                       face = "plain"), # y-axis text in fontsize 10
            axis.ticks.length = unit(0.25, "cm"))
  }
  
  message("Create Samplewise Barplot Started...")
  
  return (p)
}