#' create_cohortwise_identifier_boxplot
#'
#' Makes cohortwise boxplot for selected identifiers
#'
#' @param sample_raw_mat sample_raw_mat matrix/dataframe containing raw values
#' @param metadata_df A dataframe with samples to cohort mapping
#' @param id_name A vector of id's for which barplot to be made
#' @param cohorts_order The order of cohorts
#' @param cohort_col A cohort column present in metadata
#' @param x_label Label x-axis
#' @param y_label Label y-axis
#' @param title_label Title of the plot along with the metabolite name
#' @param compare_cohorts Logical input whether to compare cohortts and draw comparison lines and  p values accross specific user defined criteria
#' @param cohort_comparisons A list of cohort comparisons defined by users
#' 
#' 
#' @return ggplot object or plotly object
#' @examples
#' create_cohortwise_identifier_barplot(sample_raw_mat = NULL, metadata_df = NULL, id_name = NULL, 
#'                           cohorts_order = NULL, cohort_col='Cohort', 
#'                           x_label = "", y_label = "", title_label = "")
#' @import dplyr stringr ggplot2 plotly
#' @export
create_cohortwise_identifier_boxplot <- function (sample_raw_mat = NULL, 
                                                      metadata_df = NULL,
                                                      id_name = NULL, 
                                                      cohorts_order = NULL, 
                                                      cohort_col = "Cohort",
                                                      x_label = "", 
                                                      y_label = "", 
                                                      title_label = "",
                                                      interactive = TRUE,
                                                      compare_cohorts = TRUE,
                                                      cohort_comparisons = list(c("Cohort_1", "Cohort_2"), 
                                                                                c("Cohort_5", "Cohort_3"))) {
  message("Create Cohortwise Identifier Boxplot Started...")
  x = c("dplyr",
  "stringr",
  "ggplot2",
  "plotly",
  "ggpubr")
  suppressPackageStartupMessages(lapply(x, require, character.only = TRUE)
)
  
  if (interactive== TRUE & compare_cohorts ==TRUE) {
    stop("\n The function geom_GeomSignif()(i.e used to comapre and draw p value significance asterix on top of plots) has yet to be implemented in ``plotly'', the package that is used to create interactive plots.\n")
  }
  
  if ( identical(title_label, NULL)) {
    title_label <- paste(id_name,collapse=" | ")
  }
  
  if (identical(id_name, NULL)) {
    warning("The id_name is null, please choose a single metabolite name from sample_raw_mat rownames")
    return(NULL)
  }
  
  if (compare_cohorts ==TRUE & length(id_name)>1) {
    warning("\nThis function draws boxplot for only SINGLE Metabolite at a time, Please provide single input in id_name argument\n")
    return(NULL) 
  }
 
  if (!(id_name %in% rownames(sample_raw_mat))) {
    warning(c(cohort_col, " \nNot a valid Metabolite Id, please choose from input matrix rownames"))
    return(NULL)
  }
  
  
  if (!(cohort_col %in% colnames(metadata_df))) {
    warning(c(cohort_col, " is not a valid cohort_col, please choose from metadata_df colnames"))
    return(NULL)
  }
  
  rownames(sample_raw_mat) <- stringr::str_trim(rownames(sample_raw_mat))
  ids_vec <- unique(rownames(sample_raw_mat))
  common_ids <- c()
  if (length(id_name) != 0) {
    id_name <- stringr::str_trim(id_name)
    common_ids <- base::intersect(id_name, ids_vec)
  }
  if (length(common_ids) == 0) {
    warning("Invalid selected ids")
    return(NULL)
  }
  diff_ids <- setdiff(id_name, common_ids)
  if (length(diff_ids) != 0) {
    message(c("The following are not valid ids : ",  paste0(diff_ids, collapse = ", ")))
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
  transposed_mat <- as.data.frame(t(sample_raw_mat))
  
  transposed_mat$Sample <- rownames(transposed_mat)
  metadata_df$Sample <- metadata_df[, 1]
  metadata_df <- metadata_df[c("Sample", cohort_col)]  
  mat_with_metadata <- merge(metadata_df, transposed_mat, by = "Sample")
  mat_with_metadata <- mat_with_metadata %>% reshape2::melt(id.vars=c("Sample", cohort_col)) %>% 
    dplyr::rename_at(vars(c("variable", "value")), ~c("Id", "Value")) 
  mat_with_metadata <- mat_with_metadata  %>% 
    dplyr::filter(Id %in% common_ids, !!(sym(cohort_col)) %in% filtered_cohorts_vec)
  mat_with_metadata[[cohort_col]] <- factor(mat_with_metadata[[cohort_col]], levels = filtered_cohorts_vec)
  mat_with_metadata[["Id"]] <- factor(mat_with_metadata[["Id"]], levels = common_ids)
  mat_with_metadata <- mat_with_metadata[order(mat_with_metadata[["Id"]]),]  
  
  
  if (compare_cohorts == TRUE) {
    
    if (!(id_name %in% rownames(sample_raw_mat))) {
      warning(c(cohort_col, " \nNot a valid Metabolite Id, please choose from input matrix rownames\n"))
      return(NULL)
    }
    
    p <- ggplot(mat_with_metadata, aes(x = !!(sym(cohort_col)), y = Value, fill = !!(sym(cohort_col)))) + 
      #facet_wrap(~Id, scales = "free")+
      geom_boxplot(width = 0.6) +
      stat_boxplot(geom = 'errorbar', width = 0.6) +
      geom_point(aes(fill = !!(sym(cohort_col))), size =1.5, shape = 21, position = position_jitterdodge()) +
      xlab("")+
      ylab("")+
      theme(legend.title = element_blank())+
      labs(x = x_label, y = y_label) + 
      #scale_x_discrete(limits = common_ids) +
      ggsci::scale_color_aaas() + 
      theme(axis.line = element_line(size = 1, colour = "black"), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.border = element_blank(), panel.background = element_blank(), 
            plot.title = element_text(colour = "black", size = 16, 
                                      face = "bold", hjust=0.5), 
            axis.title = element_text(colour = "black",
                                      size = 14, face = "plain"), 
            axis.text.x = element_text(colour = "black", size = 10, angle = 90, 
                                       margin = unit(c(0.2, 0.2, 0.1, 0.1), "cm"), face = "plain"),
            axis.text.y = element_text(colour = "black", size = 10, 
                                       margin = unit(c(0.2, 0.2, 0.1, 0.1), "cm"), face = "plain"), 
            axis.ticks.length = unit(0.25, "cm"))+
      theme(legend.position="none")
    #}    
    
    
    p = p+ggtitle(title_label)
    
    #if (compare_cohorts == TRUE) {
    message("..........\nAlways make sure the cohort comparison contains the cohorts that you supply \n as input matrix to this function...\n")
    
    ymax = max(mat_with_metadata$Value)
    
    
    p = p +  
      annotate("text", x = 5.3, y = ymax/0.95, label = "ns  P > 0.05")+
      annotate("text", x = 5.3, y = ymax/0.90, label = "*  P ≤ 0.05")+
      annotate("text", x = 5.3, y = ymax/0.85, label = "**  P ≤ 0.01")+
      annotate("text", x = 5.3, y = ymax/0.80, label = "***  P ≤ 0.001")+
      stat_compare_means(comparisons = cohort_comparisons , label = "p.signif", method = "t.test")
    
  }else {
    if (interactive == TRUE) {
      p <- plot_ly()
      for (selected_cohort_name in filtered_cohorts_vec) {
        mat_cohort <- mat_with_metadata %>% dplyr::filter(!!(sym(cohort_col)) %in% selected_cohort_name)
        p <- add_trace(p, x = as.character(mat_cohort[['Id']]), y = mat_cohort[["Value"]],
                       type = "box", text = mat_cohort[['Sample']], name = selected_cohort_name, boxpoints = "all", jitter = 0.3, pointpos = 0)
      }
      
      p <- p %>% layout(title = "", margin=list(b=100), 
                        xaxis = list(title = x_label, titlefont = list(size = 16), 
                                     categoryorder = "array", 
                                     categoryarray = filtered_cohorts_vec), 
                        yaxis = list(title = y_label, titlefont = list(size = 16)),
                        boxmode = 'group', showlegend = TRUE) %>%
        plotly::config(displaylogo = FALSE, modeBarButtons = list(list("zoomIn2d"), list("zoomOut2d"), 
                                                                  list("toImage")), mathjax = "cdn")
    }
    else{
      p <- ggplot(mat_with_metadata, aes(x = Id, y = Value, fill = !!(sym(cohort_col)))) + 
        geom_boxplot(width = 0.6) +
        stat_boxplot(geom = 'errorbar', width = 0.6) +
        geom_point(aes(fill = !!(sym(cohort_col))), size =1.5, shape = 21, position = position_jitterdodge()) +
        # ggtitle(title_label) + 
        labs(x = x_label, y = y_label) + 
        scale_x_discrete(limits = common_ids) +
        ggsci::scale_color_aaas() + 
        theme(axis.line = element_line(size = 1, colour = "black"), 
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_blank(), panel.background = element_blank(), 
              plot.title = element_text(colour = "black", size = 18, 
                                        face = "plain", hjust=0.5), 
              axis.title = element_text(colour = "black",
                                        size = 14, face = "plain"), 
              axis.text.x = element_text(colour = "black", size = 10, angle = 90, 
                                         margin = unit(c(0.2, 0.2, 0.1, 0.1), "cm"), face = "plain"),
              axis.text.y = element_text(colour = "black", size = 10, 
                                         margin = unit(c(0.2, 0.2, 0.1, 0.1), "cm"), face = "plain"), 
              axis.ticks.length = unit(0.25, "cm"))
    }
  }
  

  message("plotting of Boxplot Finished...")
  
  return(p)
}
