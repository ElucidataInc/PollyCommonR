#' plot_anova
#'
#' makes the anova plot
#'
#' @param anova_data The anova data from compoute anova function
#' @param p_val_cutoff The pval cutoff
#' @param interaction_type Type of interaction from interaction column in anova data
#' @param annotate_id A vector of ids (ids from id column of anova data) to be annotated
#' @param row_desc A dataframe of row descriptors of the matrix
#' @param annotate_col A row descriptor column which is used to show names for annotated ids on plot
#' @param text_hover_col A row descriptor column which is used to hover text on plot
#' @param category_col A row descriptor column which is used to add shapes to different categories
#' @param marker_size_by_expr Show size of markers by expression values, TRUE/FALSE
#' @param marker_expr_col A numeric expression column present in anova_data
#' @param marker_size_range A numeric vector of minimum and maximum values of size of markers
#' @param marker_size Size of marker point
#' @param marker_opacity The opacity of the markers
#' @param x_label Label x-axis
#' @param y_label Label y-axis
#' @param title_label Title of the plot
#' @param plot_id source id for the plotly plot
#' @param plotly_highlight Highlight points on the plotly plot
#' @param highlight_on Highlight points when clicking the points (plotly_click), hovering on points (plotly_hover) or selecting multiple points (plotly_selected)
#' @param highlight_off Switch off highlighting points when plotly_doubleclick, plotly_deselect or plotly_relayout
#' @param highlight_persistent should selections persist (i.e., accumulate), default is TRUE
#' @param highlight_color The color of the highlighted points
#' @param highlight_opacitydim A number between 0 and 1 used to reduce the opacity of non-selected traces (by multiplying with the existing opacity)
#' @param highlight_debounce The amount of time to wait before firing an event (in milliseconds). The default of 0 means do not debounce at all. The color of the highlighted points
#' @param interactive make plot interactive (default is TRUE)
#' @return plotly or ggplot object
#' @examples
#' plot_anova(anova_data, p_val_cutoff = 0.05)
#' @import dplyr ggplot2 plotly ggsci ggrepel latex2exp
#' @export
plot_anova <- function(anova_data = NULL, p_val_cutoff = 0.05, interaction_type = NULL,
                       annotate_id = NULL, row_desc = NULL, annotate_col = NULL, 
                       text_hover_col = NULL, category_col = NULL, marker_size_by_expr = TRUE, 
                       marker_expr_col = "MaxExpr", marker_size_range = c(5, 25), marker_size = 8,  
                       marker_opacity = 0.5, x_label = NULL, y_label = NULL, title_label = NULL, 
                       plot_id = NULL, plotly_highlight = FALSE, highlight_on = "plotly_click",
                       highlight_off = "plotly_doubleclick", highlight_persistent = FALSE,
                       highlight_color = "blue", highlight_opacitydim = 0.8, highlight_debounce = 0, 
                       interactive = TRUE) {
  message("Plot Anova Started...")
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
  
  if (identical(p_val_cutoff, NULL)){
    warning("The p_val_cutoff is NULL")
    return (NULL)
  }
  else {
    p_val_cutoff <- as.numeric(p_val_cutoff)
    if (is.na(p_val_cutoff)){
      warning("The p_val_cutoff is not a numeric value") 
      return (NULL)
    }  
  }     
  
  if (identical(marker_size_by_expr, TRUE)){
    if (identical(marker_expr_col, NULL)){ 
      warning("The marker_expr_col is NULL")
      return (NULL)  
    } 
    
    if (!(marker_expr_col %in% colnames(anova_data))){
      warning(paste0("The ", marker_expr_col, " column is not present in anova_data"))  
      return (NULL)
    }
    
    if (!is.numeric(anova_data[, marker_expr_col])){
      warning(paste0("The ", marker_expr_col, " is not a numeric column of anova_data"))  
      return (NULL)  
    }  
    
    if(identical(marker_size_range, NULL)){
      warning("The marker_size_range is NULL")
      return (NULL)
    }      
    
    if (!is.numeric(marker_size_range)){
      warning("The marker_size_range is not a numeric vector") 
      return (NULL)
    }
    
    if (length(marker_size_range) != 2){
      warning("The marker_size_range is not a numeric vector of two elements") 
      return (NULL)
    }     
  }
  else {
    if(identical(marker_size, NULL)){
      warning("The marker_size is NULL")
      return (NULL)
    }
    else {  
      marker_size <- as.numeric(marker_size)
      if (is.na(marker_size)){
        warning("The marker_size is not a numeric value") 
        return (NULL)
      }  
    }  
  }
  
  if(identical(marker_opacity, NULL)){
    warning("The marker_opacity is NULL")
    return (NULL)
  }  
  else {
    marker_opacity <- as.numeric(marker_opacity)
    if (is.na(marker_opacity)){
      warning("The marker_opacity is not a numeric value") 
      return (NULL)
    }  
  }
  
  if (identical(interactive, TRUE)){
    if (identical(plotly_highlight, TRUE)){
      highlight_on_options <- c('plotly_click', 'plotly_hover', 'plotly_selected') 
      if (identical(highlight_on, NULL) || !(highlight_on %in% highlight_on_options)){
        warning(paste0("Select highlight_on from :", paste(sQuote(highlight_on_options), collapse = ", "), "\nUsing 'plotly_click' as default")) 
        highlight_on <- "plotly_click"    
      } 
      
      highlight_off_options <- c('plotly_doubleclick', 'plotly_deselect', 'plotly_relayout')
      if (identical(highlight_off, NULL) || !(highlight_off %in% highlight_off_options)){
        warning(paste0("Select highlight_off from :", paste(sQuote(highlight_off_options), collapse = ", "), "\nUsing 'plotly_doubleclick' as default")) 
        highlight_off <- "plotly_doubleclick"    
      }
      
      if (identical(highlight_persistent, NULL) || !is.logical(highlight_persistent)){
        warning("The highlight_persistent should be be TRUE/FALSE, using TRUE as default")
        highlight_persistent <- TRUE
      }
      
      color_name <- NULL  
      try(color_name <- grDevices::col2rgb(highlight_color, alpha =TRUE), silent = T)  
      if (identical(highlight_color, NULL) || identical(color_name, NULL)){
        warning("Select the valid highlight_color, using 'blue' as default")
        highlight_color <- "blue"
      }
      
      if (identical(highlight_debounce, NULL) || !is.numeric(highlight_debounce)){
        warning("Select the valid highlight_debounce, using 0 as default")
        highlight_debounce <- 0
      }
      
      if (identical(highlight_opacitydim, NULL) || !is.numeric(highlight_opacitydim)){
        warning("Select the valid highlight_opacitydim, using 0.8 as default")
        highlight_opacitydim <- 0.8
      }        
    } 
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
  
  ave_expr <- as.numeric(anova_data[, marker_expr_col])
  ave_expr <- ave_expr[is.finite(ave_expr)]
  ave_expr_range <- c(min(ave_expr), max(ave_expr))
  if (identical(marker_size_by_expr, TRUE)){
    slope_m <- (marker_size_range[2] - marker_size_range[1]) / (ave_expr_range[2] - ave_expr_range[1])
    eq_constant <- marker_size_range[2] - (slope_m * ave_expr_range[2])                                           
    anova_data$marker_size <- sapply(anova_data[, marker_expr_col], function(x) {
      y <- (slope_m * x) + eq_constant
      if (!is.finite(y)){ y <- marker_size_range[1]}
      return(y)
    })                                       
  }
  else { anova_data$marker_size <- marker_size}                                                                                     
  
  if (marker_expr_col %in% colnames(anova_data)){ 
    anova_data$text_hover<-  paste0(paste0("id: ", anova_data$id),
                                    "<br>", paste0("-log10(P.Value)): ", - log10(anova_data[["P.Value"]])),
                                    "<br>", paste0("P.Value: ", anova_data[["P.Value"]]),
                                    "<br>", paste0(marker_expr_col, ": ", anova_data[, marker_expr_col]))
  }
  else {
    anova_data$text_hover<-  paste0(paste0("id: ", anova_data$id),
                                    "<br>", paste0("-log10(P.Value)): ", - log10(anova_data[["P.Value"]])),
                                    "<br>", paste0("P.Value: ", anova_data[["P.Value"]]))
  }
  
  if ("interaction" %in% colnames(anova_data)){
    anova_data$text_hover <- paste0(anova_data$text_hover, "<br>", paste0("interaction: ", anova_data$interaction))
  }
  
  anova_data$category_sym <- 1  
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
  y_col <- "-log10(P.Value)"                                              
  
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
    
    if (identical(plotly_highlight, TRUE)){ anova_data <- plotly::highlight_key(anova_data, ~id)}
    
    if (identical(marker_size_by_expr, TRUE)){
      p <- plotly::plot_ly(data = anova_data, source = plot_id,
                           x = stats::as.formula(paste0("~", x_col)), y = stats::as.formula(paste0("~", y_col)),
                           customdata = ~id, type = "scatter", mode = "markers", size = ~marker_size, 
                           fill = ~'', sizes = marker_size_range,
                           marker = list(sizemode = 'diameter', opacity = marker_opacity),
                           color = ~threshold, colors = significance_color,
                           symbol = ~category_sym, text = ~text_hover)
    }
    else {
      p <- plotly::plot_ly(data = anova_data, source = plot_id,
                           x = stats::as.formula(paste0("~", x_col)), y = stats::as.formula(paste0("~", y_col)),
                           customdata = ~id, type = "scatter", mode = "markers",
                           marker = list(size = marker_size, sizemode = 'diameter', opacity = marker_opacity),
                           color = ~threshold, colors = significance_color,
                           symbol = ~category_sym, text = ~text_hover)        
    }
    
    p <- p %>% 
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
    
    if (identical(plotly_highlight, TRUE)){
      p <- plotly::highlight(p, on = highlight_on, off = highlight_off, persistent = highlight_persistent,
                             color = highlight_color, opacityDim = highlight_opacitydim, selected = plotly::attrs_selected(showlegend = FALSE), 
                             debounce = highlight_debounce)
    }
    
  } else {
    if (identical(annotate_col, NULL)){ annotate_col <- "id" }
    if (identical(x_label, NULL)){ x_label <- xaxis_lab_gg }
    if (identical(y_label, NULL)){ y_label <- yaxis_lab_gg }      
    anova_data[!is.finite(anova_data[, marker_expr_col]), marker_expr_col] <- ave_expr_range[1]
    p <- ggplot(anova_data, aes_string(x = x_col, y = y_col, color = "threshold", fill = "threshold", shape = category_col), text = id)
    if (identical(marker_size_by_expr, TRUE)){
      p <- p + geom_point(aes_string(size = marker_expr_col), alpha = marker_opacity) + 
        ggplot2::scale_size(range = marker_size_range/3)    
    }
    else {
      p <- p + geom_point(size = marker_size/2, alpha = marker_opacity) +
        ggplot2::guides(size = FALSE)
    }         
    p <- p+  
      labs(x = x_label, y = y_label, title = title_label, size = "Expression", color = "Significance", fill = "Significance", shape = "Category") + # x and y axis labels
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
  
  message("Plot Anova Completed...")
  
  return (p)
}