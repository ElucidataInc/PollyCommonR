#' create_cohortwise_identifier_barplot
#'
#' Makes cohortwise barplot for selected identifiers
#'
#' @param sample_raw_mat sample_raw_mat matrix/dataframe containing raw values
#' @param metadata_df A dataframe with samples to cohort mapping
#' @param id_order A vector of id's for which barplot to be made
#' @param cohorts_order The order of cohorts
#' @param cohort_col A cohort column present in metadata
#' @param x_label Label x-axis
#' @param y_label Label y-axis
#' @param title_label Title of the plot
#' @param x_text_wrap_n A number to wrap x axis text to next line after every nth position
#' @param x_text_angle The x-axis text angle
#' @param y_text_angle The y-axis text angle
#' @param interactive Make plot interactive using plotly
#' @return ggplot object or plotly object
#' @examples
#' create_cohortwise_identifier_barplot(sample_raw_mat = NULL, metadata_df = NULL, id_order = NULL, 
#'                           cohorts_order = NULL, cohort_col='Cohort', 
#'                           x_label = "", y_label = "", title_label = "")
#' @import dplyr stringr ggplot2 plotly
#' @export
create_cohortwise_identifier_barplot <- function (sample_raw_mat = NULL, metadata_df = NULL, 
                                                  id_order = NULL, cohorts_order = NULL, 
                                                  cohort_col = "Cohort", x_label = NULL, y_label = NULL,
                                                  title_label = NULL, x_text_wrap_n = NULL, x_text_angle = NULL,
                                                  y_text_angle = NULL, interactive = FALSE) {
  message("Create Cohortwise Identifier Barplot Started...")
  require(dplyr)
  require(stringr)
  require(ggplot2)
  require(plotly)
  
  if (identical(id_order, NULL)) {
    warning("The id_order is null, please choose from sample_raw_mat rownames")
    return(NULL)
  }
  
  if (!(cohort_col %in% colnames(metadata_df))) {
    warning(c(cohort_col, " is not a valid cohort_col, please choose from metadata_df colnames"))
    return(NULL)
  }
  
  rownames(sample_raw_mat) <- stringr::str_trim(rownames(sample_raw_mat))
  ids_vec <- unique(rownames(sample_raw_mat))
  common_ids <- c()
  if (length(id_order) != 0) {
    id_order <- stringr::str_trim(id_order)
    common_ids <- base::intersect(id_order, ids_vec)
  }
  if (length(common_ids) == 0) {
    warning("Invalid selected ids")
    return(NULL)
  }
  diff_ids <- setdiff(id_order, common_ids)
  if (length(diff_ids) != 0) {
    message(c("The following are not valid ids : ",  paste0(diff_ids, collapse = ", ")))
  }
  
  if (identical(x_label, NULL)) {
    x_label <- ""
  }
  if (identical(y_label, NULL)) {
    y_label <- ""
  } 
  
  if (identical(title_label, NULL)) {
    title_label <- ""
  }
  
  if(!identical(x_text_wrap_n, NULL)){
    x_text_wrap_n <- as.numeric(x_text_wrap_n)
    if (is.na(x_text_wrap_n)){
      warning("The x_text_wrap_n is not a numeric value") 
      return (NULL)
    }  
  }
  
  if(!identical(x_text_angle, NULL)){
    x_text_angle <- as.numeric(x_text_angle)
    if (is.na(x_text_angle)){
      warning("The x_text_angle is not a numeric value") 
      return (NULL)
    }  
  }
  
  if(!identical(y_text_angle, NULL)){
    y_text_angle <- as.numeric(y_text_angle)
    if (is.na(y_text_angle)){
      warning("The y_text_angle is not a numeric value") 
      return (NULL)
    }  
  }    
  
  metadata_df[[cohort_col]] <- stringr::str_trim(metadata_df[[cohort_col]])
  cohorts_vec <- unique(metadata_df[[cohort_col]])
  filtered_cohorts_vec <- cohorts_vec
  common_cohorts <- c()
  if (length(cohorts_order) != 0) {
    cohorts_order <- stringr::str_trim(cohorts_order)
    common_cohorts <- base::intersect(cohorts_order, cohorts_vec)
  }
  if (length(common_cohorts) != 0) {
    filtered_cohorts_vec <- common_cohorts
  }
  
  sample_raw_mat <- sample_raw_mat[common_ids, , drop = FALSE]
  
  if (!identical(x_text_wrap_n, NULL)){
    row.names(sample_raw_mat) <- suppressMessages(sapply(row.names(sample_raw_mat), function(x) paste(PollyCommonR::split_string_every_nth_char(x, x_text_wrap_n), collapse = "-\n")))
    common_ids <- suppressMessages(sapply(common_ids, function(x) paste(PollyCommonR::split_string_every_nth_char(x, x_text_wrap_n), collapse = "-\n")))
  }
  
  transposed_mat <- as.data.frame(t(sample_raw_mat))
  
  transposed_mat$Sample <- rownames(transposed_mat)
  metadata_df$Sample <- metadata_df[, 1]
  metadata_df <- metadata_df[c("Sample", cohort_col)]  
  mat_with_metadata <- merge(metadata_df, transposed_mat, by = "Sample")
  mat_with_metadata <- mat_with_metadata %>% reshape2::melt(id.vars=c("Sample", cohort_col)) %>% 
    dplyr::rename_at(vars(c("variable", "value")), ~c("Id", "Value"))
  mat_with_metadata_mean_sd <- mat_with_metadata %>% dplyr::group_by_at(c(cohort_col, "Id")) %>%
    dplyr::summarize(mean = mean(Value, na.rm = TRUE), sd = sd(Value, na.rm = TRUE))
  
  filtered_mat_mean_sd <- mat_with_metadata_mean_sd %>% 
    dplyr::filter(Id %in% common_ids, !!(sym(cohort_col)) %in% filtered_cohorts_vec)
  filtered_mat_mean_sd[[cohort_col]] <- factor(filtered_mat_mean_sd[[cohort_col]], levels = filtered_cohorts_vec)
  filtered_mat_mean_sd[["Id"]] <- factor(filtered_mat_mean_sd[["Id"]], levels = common_ids)
  filtered_mat_mean_sd <- filtered_mat_mean_sd[order(filtered_mat_mean_sd[["Id"]]),]
  
  if (interactive == TRUE) {
    p <- plot_ly()
    for (selected_cohort_name in filtered_cohorts_vec) {
      mat_cohort <- filtered_mat_mean_sd %>% dplyr::filter(!!(sym(cohort_col)) %in% selected_cohort_name)
      p <- add_trace(p, x = as.character(mat_cohort[['Id']]), y = mat_cohort[["mean"]], 
                     error_y = list(array = mat_cohort[["sd"]], color = '#000000'),
                     type = "bar", name = selected_cohort_name)
    }
    
    p <- p %>% layout(title = title_label, margin=list(b=100), 
                      xaxis = list(title = x_label, titlefont = list(size = 16), 
                                   categoryorder = "array",
                                   categoryarray = filtered_cohorts_vec), 
                      yaxis = list(title = y_label, titlefont = list(size = 16)),
                      showlegend = TRUE) %>%
      plotly::config(displaylogo = FALSE, modeBarButtons = list(list("zoomIn2d"), list("zoomOut2d"), 
                                                                list("toImage")), mathjax = "cdn")
  }
  else{
    if (identical(x_text_angle, NULL)){
      if (length(common_ids) == 1){ x_text_angle <- 0}
      else { x_text_angle <- 90} 
    }
    
    p <- ggplot(filtered_mat_mean_sd, aes(x = Id, y = mean, fill = !!(sym(cohort_col)))) + 
      geom_bar(stat = "identity", position = "dodge") + 
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
      labs(x = x_label, y = y_label, title = title_label) + 
      scale_x_discrete(limits = common_ids) +
      theme(axis.line = element_line(size = 1, colour = "black"), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.border = element_blank(), panel.background = element_blank(), 
            plot.title = element_text(colour = "black", size = 18, 
                                      face = "plain", hjust=0.5), 
            axis.title = element_text(colour = "black",
                                      size = 14, face = "plain"), 
            axis.text.x = element_text(colour = "black", size = 10, angle = x_text_angle, 
                                       margin = unit(c(0.2, 0.2, 0.1, 0.1), "cm"), face = "plain"),
            axis.text.y = element_text(colour = "black", size = 10,  angle = y_text_angle,
                                       margin = unit(c(0.2, 0.2, 0.1, 0.1), "cm"), face = "plain"), 
            axis.ticks.length = unit(0.25, "cm"))
  }
  
  message("Create Cohortwise Identifier Barplot Started...")
  
  return(p)
}