#' plot_anova
#'
#' makes a anova plot
#'
#' @param anova_data The anova data from compoute anova function
#' @param p_val_cutoff The pval cutoff
#' @param interaction_type Type of interaction from interaction column in anova data
#' @param annotate_id A vector of ids (ids from id column of anova data) to be annotated
#' @param row_desc A dataframe of row descriptors of the matrix
#' @param annotate_col A row descriptor column which is used to show names for annotated ids on plot
#' @param text_hover_col A row descriptor column which is used to hover text on plot
#' @param category_col A row descriptor column which is used to add shapes to different categories
#' @param x_label Label x-axis
#' @param y_label Label y-axis
#' @param title_label Title of the plot
#' @param marker_size Size of marker point
#' @param plot_id source id for the plotly plot
#' @param interactive make plot interactive (default is TRUE)
#' @return plotly or ggplot object
#' @examples
#' plot_anova(diff_exp, p_val_cutoff = 0.05)
#' @import dplyr ggplot2 plotly ggsci ggrepel latex2exp
#' @export
plot_anova <- function(anova_data = NULL, p_val_cutoff = NULL, interaction_type = NULL,
                       annotate_id = NULL, row_desc = NULL, annotate_col = NULL, 
                       text_hover_col = NULL, category_col = NULL, x_label = NULL, 
                       y_label = NULL, title_label = NULL, marker_size = 8, plot_id = NULL, 
                       interactive = TRUE) {
  message("Make Volcano Plot Started...")
  require(dplyr)  
  require(ggplot2)
  require(plotly)
  require(ggsci)
  require(ggrepel)
  require(latex2exp)
  
  if (identical(anova_data, NULL)){
    warning("The anova_data is NULL")
    return (NULL)
  }
  
  if (!identical(class(anova_data), "data.frame")){
    warning("The anova_data is not a dataframe, please provide valid anova_data")
    return(NULL) 
  }     
  
  if (!("id" %in% colnames(anova_data))){
    warning("The anova_data should contain 'id' column")
    return (NULL)
  }
  
  if (!("P.Value" %in% colnames(anova_data))){
    warning("The anova_data should contain 'P.Value' column")
    return (NULL)
  }  
  
  if (!identical(row_desc, NULL)){
    if (!identical(class(row_desc), "data.frame")){
      warning("The row_desc is not a dataframe, please provide valid row_desc")
      return(NULL) 
    }   
    
    if (!identical(annotate_col, NULL)){
      if (!(annotate_col %in% colnames(row_desc))){
        warning(c(annotate_col, " column is not present in row_desc"))
        return (NULL)
      } 
    }      
    
    if (!identical(text_hover_col, NULL)){
      if (!(text_hover_col %in% colnames(row_desc))){
        warning(c(text_hover_col, " column is not present in row_desc"))
        return (NULL)
      }      
    }
    
    if (!identical(category_col, NULL)){
      if (!(category_col %in% colnames(row_desc))){
        warning(c(category_col, " column is not present in row_desc"))
        return (NULL)
      }
    }
    
    anova_data <- base::merge(anova_data, row_desc, by.x = "id", by.y=0, sort = FALSE)
    if (nrow(anova_data) < 1){
      warning("No common data between anova_data ('id' column) and row_desc ('rownames')")
      return (NULL)
    }   
  }
  else {
    annotate_col <- "id"
    text_hover_col <- NULL
    category_col <- NULL    
  }
  
  if (!identical(annotate_id, NULL)){
    common_annotate_id <- unique(base::intersect(annotate_id, anova_data$id))
    if (length(common_annotate_id) < 1){
      warning("The annotate ids are not matched with id column of anova_data")    
    }
    else {
      diff_annotate_id <- unique(base::setdiff(annotate_id, anova_data$id))
      if (length(diff_annotate_id) >= 1){
        warning(paste("The following annotate ids are not matched with id column of anova_data :", paste(sQuote(diff_annotate_id), collapse = ", "), collapse = " ")) 
      }  
    }       
  }
  
  if ("interaction" %in% colnames(anova_data)){
    all_interaction <- unique(anova_data[, "interaction"])
    if (!identical(interaction_type, NULL)){
      diff_interaction <- base::setdiff(interaction_type, all_interaction)  
      if (length(diff_interaction) > 0){
        warning(paste("The following interaction types are not present in interaction column of anova data :", paste(sQuote(diff_interaction), collapse = ", "), collapse = " "))    
        return (NULL)
      }
      anova_data <- anova_data[anova_data[, "interaction"] %in% interaction_type, ] 
    }
  }
  
  id_index_df <- data.frame(index = 1:length(unique(anova_data$id)), id = unique(anova_data$id), stringsAsFactors = FALSE, check.names = FALSE)
  anova_data <- merge(anova_data, id_index_df, by = "id", sort = FALSE)
  anova_data$threshold <- "not significant"
  
  anova_data <- dplyr::mutate(anova_data, threshold = dplyr::case_when(P.Value <= p_val_cutoff ~ "significant", TRUE ~ threshold))
  
  significance_table <- base::table(anova_data$threshold)
  print (significance_table)
  
  anova_data <- anova_data[!is.infinite(anova_data[, "P.Value"]), ]
  anova_data <- anova_data[!sapply(anova_data[, "P.Value"], anyNA), ]
  
  if (nrow(anova_data) < 1){
    warning("The anova data has zero valid rows")
    return (NULL)
  }
  
  if (identical(row_desc, NULL) | identical(category_col, NULL)){
    for (sig_type in unique(anova_data$threshold)){
      anova_data[anova_data$threshold == sig_type, "threshold"] <- paste0(sig_type, " (", significance_table[[sig_type]], ")") 
    }
  }
  
  if ("interaction" %in% colnames(anova_data)){
    anova_data$text_hover <- paste0(anova_data$id, "\n", "interaction: ", anova_data$interaction)
  } else { anova_data$text_hover <- anova_data$id}
  anova_data$category_sym <- NULL  
  if (!identical(row_desc, NULL)){
    if (!identical(text_hover_col, NULL)){
      anova_data$text_hover <- anova_data[[text_hover_col]]
    } 
    if (!identical(category_col, NULL)){
      anova_data[[category_col]][unlist(lapply(anova_data[[category_col]], function(x) x %in% c(NA, "NA", "")))] <- "NA"
      anova_data$category_sym <- anova_data[[category_col]]                                             
    }
  }
  
  significance_color <- c("grey", "red")  
  xaxis_lab_gg <- latex2exp::TeX("$Features$")
  xaxis_lab_pl <- plotly::TeX("\\text{Features}")
  yaxis_lab_gg <- latex2exp::TeX("$-\\log_{10}(p \\, value)$")
  yaxis_lab_pl <- plotly::TeX("-\\log_{10}(\\text{p value})")
  
  x_col = "index"
  y_col = "-log10(P.Value)" 
  x_val <- anova_data$index
  y_val <- -log10(anova_data[, "P.Value"])                                              
  
  if (interactive == TRUE){
    if (identical(x_label, NULL)){ x_label <- xaxis_lab_pl }
    if (identical(y_label, NULL)){ y_label <- yaxis_lab_pl }
    
    filtered_anova_data <- anova_data[anova_data$id %in% annotate_id, ]
    if(nrow(filtered_anova_data) == 0){
      a <- NULL
    } else {
      annotate_text <- filtered_anova_data$id  
      if (!identical(row_desc, NULL)){  
        if (!identical(annotate_col, NULL)){ annotate_text <- filtered_anova_data[[annotate_col]] }
      }
      
      a <- list(
        x = filtered_anova_data$index,
        y = -log10(filtered_anova_data[, "P.Value"]),
        text = annotate_text,
        xref = "x",
        yref = "y",
        showarrow = TRUE,
        arrowhead = 3.5,
        ax = 20,
        ay = -40
      )}
    p <- plotly::plot_ly(source = plot_id) %>%
      add_trace(
        x = x_val, y = y_val, customdata = anova_data$id,
        type = "scattergl", mode = "markers",
        marker = list(size = marker_size),
        color = anova_data$threshold,
        colors = significance_color,
        symbol = anova_data$category_sym,
        text = anova_data$text_hover 
      ) %>%
      
      layout(
        title = list(text = title_label, xref = "paper", yref = "paper"),
        yaxis = list(title = y_label),
        xaxis = list(title = x_label),
        annotations = a,
        showlegend = TRUE
      ) %>%
      add_annotations(text="Significance", xref="paper", yref="paper",
                      x=1.04, xanchor="left",
                      y=0.98, yanchor="bottom",
                      font = list(size = (23.91034/18)*13),
                      legendtitle=TRUE, showarrow=FALSE ) %>% 
      plotly::config(displaylogo = FALSE,
                     modeBarButtons = list(list("zoom2d"),
                                           list("select2d"),
                                           list("lasso2d"),
                                           list("autoScale2d"),
                                           list("resetScale2d"),
                                           list("pan2d"),
                                           list("zoomIn2d"), 
                                           list("zoomOut2d"),
                                           list("hoverClosestCartesian"),
                                           list('toImage')), 
                     mathjax = 'cdn')
  } else {
    if (identical(annotate_col, NULL)){ annotate_col <- "id" }
    if (identical(x_label, NULL)){ x_label <- xaxis_lab_gg }
    if (identical(y_label, NULL)){ y_label <- yaxis_lab_gg }      
    p <- ggplot(anova_data, aes_string(x = x_col, y = y_col, color = "threshold", fill = "threshold", shape = category_col), text = id) + 
      geom_point( size = marker_size/2, alpha = 0.7) + # scatter plot function with shape of points defined as 21 scale.   
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
      ggtitle(title_label) +       
      labs(x = x_label, y = y_label, color = "Significance", fill = "Significance", shape = "Category") + # x and y axis labels
      theme(legend.position = "right", legend.direction = "vertical", # legend positioned at the bottom, horizantal direction,
            axis.line = element_line(size=1, colour = "black"),	# axis line of size 1 inch in black color
            panel.grid.major = element_blank(),	# major grids included
            panel.grid.minor = element_blank(),	# no minor grids
            panel.border = element_blank(), panel.background = element_blank(), # no borders and background color
            plot.title = element_text(colour="black", size = 18, face = "plain", hjust=0.5),
            axis.title = element_text(colour="black", size = 15, face = "bold"), # axis title 
            axis.text.x = element_text(colour="black", size = 10, margin=unit(c(0.2,0.2,0.1,0.1), "cm"), face = "bold"), # x-axis text in fontsize 10
            axis.text.y = element_text(colour="black", size = 10, margin=unit(c(0.2,0.2,0.1,0.1), "cm"), face = "bold"), # y-axis text in fontsize 10
            legend.text = element_text(size = 10, face = "bold"),
            legend.title = element_text(colour="black", size=12, face="bold"),
            axis.ticks.length = unit(0.25, "cm")) # ticks facing inward with 0.25cm length
    p <- p + ggrepel::geom_text_repel(data = subset(anova_data, id %in% annotate_id), aes_string(label = annotate_col), color = "black", size = 3,
                                      box.padding = unit(0.5, 'lines'),
                                      point.padding = unit(1.6, 'lines'),
                                      segment.size = 0.5,
                                      arrow = arrow(length = unit(0.01, 'npc')),
                                      force = 0.5,
                                      max.iter = 3e3)
  }
  
  message("Make Volcano Plot Completed...")
  
  return (p)
}