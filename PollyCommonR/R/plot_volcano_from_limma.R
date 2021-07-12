#' plot_volcano_from_limma
#'
#' makes a volcano plot
#'
#' @param diff_exp The differential expression dataframe from limma
#' @param log2fc_range The absolute log2FC value to set threshold
#' @param p_val_cutoff The pval cutoff
#' @param p_val_type The pval type (P.Value or adj.P.Val)
#' @param annotate_id A vector of ids (rownames of row_desc) to be annotated
#' @param row_desc A dataframe of row descriptors of the matrix
#' @param annotate_col A row descriptor column which is used to show names for annotated ids on plot
#' @param text_hover_col A row descriptor column which is used to hover text on plot
#' @param category_col A row descriptor column which is used to add shapes to different categories
#' @param marker_size_by_expr Show size of markers by average expression values (AveExpr), TRUE/FALSE
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
#' plot_volcano_from_limma(diff_exp, log2fc_range = 0.5, p_val_cutoff = 0.05, p_val_type = "P.Value")
#' @import dplyr ggplot2 plotly ggsci ggrepel latex2exp
#' @export
plot_volcano_from_limma <- function(diff_exp = NULL, log2fc_range = 1, p_val_cutoff = 0.05, 
                                    p_val_type = "P.Value", annotate_id = NULL, row_desc = NULL, 
                                    annotate_col = NULL, text_hover_col = NULL, category_col = NULL,
                                    marker_size_by_expr = TRUE, marker_size_range = c(5, 25),
                                    marker_size = 8,  marker_opacity = 0.5, x_label = NULL, y_label = NULL,
                                    title_label = NULL, plot_id = NULL, plotly_highlight = FALSE, 
                                    highlight_on = "plotly_click", highlight_off = "plotly_doubleclick", 
                                    highlight_persistent = FALSE, highlight_color = "blue", 
                                    highlight_opacitydim = 0.8, highlight_debounce = 0, interactive = TRUE) {
  message("Make Volcano Plot Started...")
  require(dplyr)
  require(ggplot2)
  require(plotly)
  require(ggsci)
  require(ggrepel)
  require(latex2exp)
  
  if (identical(diff_exp, NULL)){
    warning("The diff_exp is NULL")
    return (NULL)
  }
  
  if (!identical(class(diff_exp), "data.frame")){
    warning("The diff_exp is not a dataframe, please provide valid diff_exp")
    return(NULL) 
  }
  
  if (identical(log2fc_range, NULL)){
    warning("The log2fc_range is NULL")
    return (NULL)
  }
  else {
    log2fc_range <- as.numeric(log2fc_range)
    if (is.na(log2fc_range)){
      warning("The log2fc_range is not a numeric value") 
      return (NULL)
    }  
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
  
  pval_type <- c("P.Value", "adj.P.Val")
  common_pval <- base::intersect(pval_type, colnames(diff_exp))
  if (length(common_pval) < 1){
    warning(paste0("The diff_exp should have at least one or both pval type from ", paste0(pval_type, collapse = " and ")))
    return (NULL)
  }  
  
  if (!("logFC" %in% colnames(diff_exp))){
    warning("The diff_exp should contain 'logFC' column")
    return (NULL)
  }
  
  required_cols <- c("logFC", common_pval)  
  
  if (identical(p_val_type, NULL)){
    p_val_type <- pval_type[1]  
    warning("The p_val_type is NULL, using 'P.Value' as default")
  }
  
  if (!(p_val_type %in% pval_type)){
    warning("Invalid p_val_type, select 'P.Value' or 'adj.P.Val'")
    return(NULL)  
  }
  
  if (!(p_val_type %in% colnames(diff_exp))){
    warning(paste0("The ", p_val_type, " is not present in diff_exp"))
    return (NULL)
  }
  
  if (identical(marker_size_by_expr, TRUE)){
    if (!("AveExpr" %in% colnames(diff_exp))){
      warning("The average expression (AveExpr) column is not present in diff_exp")
      return (NULL)
    }
    
    if (identical(marker_size_range, NULL)){ 
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
    
    if (!all(row.names(diff_exp) %in% row.names(row_desc))){
      warning("Not all rownames of diff_exp are present in rownames of row_desc")
      return (NULL)  
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
    
    diff_exp <- base::transform(base::merge(diff_exp,row_desc,by=0, sort = FALSE), row.names=Row.names, Row.names=NULL) 
  }
  else {
    annotate_col <- "id"
    text_hover_col <- NULL
    category_col <- NULL    
  }  
  
  if (!identical(annotate_id, NULL)){
    common_annotate_id <- base::intersect(annotate_id, row.names(diff_exp))
    if (length(common_annotate_id) < 1){
      warning("The annotate ids are not matched with rownames of diff_exp")    
    }
    else {
      diff_annotate_id <- base::setdiff(annotate_id, row.names(diff_exp))
      if (length(diff_annotate_id) >= 1){
        warning(paste("The following annotate ids are not matched with rownames of diff_exp :", paste(sQuote(diff_annotate_id), collapse = ", "), collapse = " ")) 
      }  
    }       
  }
  
  diff_exp$id <- row.names(diff_exp)
  diff_exp$threshold <- "not significant"
  
  diff_exp <- dplyr::mutate(diff_exp, threshold = dplyr::case_when((abs(logFC) >= log2fc_range) & (!!(dplyr::sym(p_val_type)) <= p_val_cutoff) ~ "significant", TRUE ~ threshold))
  row.names(diff_exp) <- diff_exp$id
  significance_table <- base::table(diff_exp$threshold)
  print (significance_table)
  
  diff_exp <- diff_exp[!is.infinite(rowSums(diff_exp[, required_cols])), ]
  diff_exp <- diff_exp[!apply(diff_exp[, required_cols], 1, anyNA), ]
  
  if (nrow(diff_exp) < 1){
    warning("Differential Expression dataframe has zero valid rows")
    return (NULL)
  }
  
  if (identical(row_desc, NULL) | identical(category_col, NULL)){
    for (sig_type in unique(diff_exp$threshold)){
      diff_exp[diff_exp$threshold == sig_type, "threshold"] <- paste0(sig_type, " (", significance_table[[sig_type]], ")") 
    }
  }
  
  ave_expr <- diff_exp$AveExpr
  ave_expr <- ave_expr[is.finite(ave_expr)]
  ave_expr_range <- c(min(ave_expr), max(ave_expr))
  if (identical(marker_size_by_expr, TRUE)){
    slope_m <- (marker_size_range[2] - marker_size_range[1]) / (ave_expr_range[2] - ave_expr_range[1])
    eq_constant <- marker_size_range[2] - (slope_m * ave_expr_range[2])                                           
    diff_exp$marker_size <- sapply(diff_exp$AveExpr, function(x) {
      y <- (slope_m * x) + eq_constant
      if (!is.finite(y)){ y <- marker_size_range[1]}
      return(y)
    })                                       
  }
  else { diff_exp$marker_size <- marker_size}                                                                                     
  
  if ("AveExpr" %in% colnames(diff_exp)){ 
    diff_exp$text_hover<-  paste0(paste0("id: ", diff_exp$id),
                                  "<br>", paste0("logFC: ", diff_exp$logFC),
                                  "<br>", paste0(gsub("pval", p_val_type, "-log10(pval)"), ": ", - log10(diff_exp[[p_val_type]])),
                                  "<br>", paste0(p_val_type, ": ", diff_exp[[p_val_type]]),
                                  "<br>", paste0("Average Expression (AveExpr): ", diff_exp$AveExpr))
  }
  else {
    diff_exp$text_hover<-  paste0(paste0("id: ", diff_exp$id),
                                  "<br>", paste0("logFC: ", diff_exp$logFC),
                                  "<br>", paste0(gsub("pval", p_val_type, "-log10(pval)"), ": ", - log10(diff_exp[[p_val_type]])),
                                  "<br>", paste0(p_val_type, ": ", diff_exp[[p_val_type]]))
  } 
  
  diff_exp$category_sym <- 1
  if (!identical(row_desc, NULL)){
    if (!identical(text_hover_col, NULL)){
      diff_exp$text_hover <- diff_exp[[text_hover_col]]
    } 
    if (!identical(category_col, NULL)){
      diff_exp[[category_col]][unlist(lapply(diff_exp[[category_col]], function(x) x %in% c(NA, "NA", "")))] <- "NA"
      diff_exp$category_sym <- diff_exp[[category_col]]                                             
    }
  }
  
  significance_color <- c("grey", "red")  
  xaxis_lab_gg <- latex2exp::TeX("$\\log_{2}(fold \\, change)$")
  xaxis_lab_pl <- plotly::TeX("\\log_{2}(\\text{fold change})")                                         
  if(identical(p_val_type, "P.Value")){
    yaxis_lab_gg <- latex2exp::TeX("$-\\log_{10}(p \\, value)$")
    yaxis_lab_pl <- plotly::TeX("-\\log_{10}(\\text{p value})")
  }
  else if(identical(p_val_type, "adj.P.Val")){
    yaxis_lab_gg <- latex2exp::TeX("$-\\log_{10}(adjusted \\,p \\, value)$")
    yaxis_lab_pl <- plotly::TeX("-\\log_{10}(\\text{adjusted p value})")  
  }
  
  x_col = "logFC"
  y_col <- paste0("-log10(", p_val_type,")")                                          
  
  if (interactive == TRUE){
    if (identical(x_label, NULL)){ x_label <- xaxis_lab_pl }
    if (identical(y_label, NULL)){ y_label <- yaxis_lab_pl }
    
    filtered_diff_exp <- diff_exp[diff_exp$id %in% annotate_id, ]
    if(nrow(filtered_diff_exp) == 0){
      a <- NULL
    } else {
      annotate_text <- filtered_diff_exp$id  
      if (!identical(row_desc, NULL)){  
        if (!identical(annotate_col, NULL)){ annotate_text <- filtered_diff_exp[[annotate_col]] }
      }
      
      a <- list(
        x = filtered_diff_exp$logFC,
        y = -log10(filtered_diff_exp[, p_val_type]),
        text = annotate_text,
        xref = "x",
        yref = "y",
        showarrow = TRUE,
        arrowhead = 3.5,
        ax = 20,
        ay = -40
      )}
    
    if (identical(plotly_highlight, TRUE)){ diff_exp <- plotly::highlight_key(diff_exp, ~id)}  
    
    if (identical(marker_size_by_expr, TRUE)){
      p <- plotly::plot_ly(data = diff_exp, source = plot_id,
                           x = stats::as.formula(paste0("~", x_col)), y = stats::as.formula(paste0("~", y_col)),
                           customdata = ~id, type = "scatter", mode = "markers", size = ~marker_size, 
                           fill = ~'', sizes = marker_size_range,
                           marker = list(sizemode = 'diameter', opacity = marker_opacity),
                           color = ~threshold, colors = significance_color,
                           symbol = ~category_sym, text = ~text_hover)
    }
    else {
      p <- plotly::plot_ly(data = diff_exp, source = plot_id,
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
    
    diff_exp[!is.finite(diff_exp$AveExpr), "AveExpr"] <- ave_expr_range[1]
    p <- ggplot(diff_exp, aes_string(x = x_col, y = y_col, color = "threshold", fill = "threshold", shape = category_col), text = id)
    if (identical(marker_size_by_expr, TRUE)){
      p <- p + geom_point(aes_string(size = "AveExpr"), alpha = marker_opacity) + 
        ggplot2::scale_size(range = marker_size_range/3)    
    }
    else {
      p <- p + geom_point(size = marker_size/2, alpha = marker_opacity) +
        ggplot2::guides(size = FALSE )
    }  
    p <- p +
      labs(x = x_label, y = y_label, title = title_label, size = "Average Expression", color = "Significance", fill = "Significance", shape = "Category") + # x and y axis labels
      theme(legend.position = "right", legend.direction = "vertical", # legend positioned at the bottom, horizantal direction,
            axis.line = element_line(size=1, colour = "black"),	# axis line of size 1 inch in black color
            panel.grid.major = element_blank(), # major grids included
            panel.grid.minor = element_blank(), # no minor grids
            panel.border = element_blank(), panel.background = element_blank(), # no borders and background color
            plot.title = element_text(colour="black", size = 18, face = "plain", hjust=0.5),
            axis.title = element_text(colour="black", size = 15, face = "bold"), # axis title 
            axis.text.x = element_text(colour="black", size = 10, margin=unit(c(0.2,0.2,0.1,0.1), "cm"), face = "bold"), # x-axis text in fontsize 10
            axis.text.y = element_text(colour="black", size = 10, margin=unit(c(0.2,0.2,0.1,0.1), "cm"), face = "bold"), # y-axis text in fontsize 10
            legend.text = element_text(size = 10, face = "bold"),
            legend.title = element_text(colour="black", size=12, face="bold"),
            axis.ticks.length = unit(0.25, "cm")) # ticks facing inward with 0.25cm length
    p <- p + ggrepel::geom_text_repel(data = subset(diff_exp, id %in% annotate_id), aes_string(label = annotate_col), color = "black", size = 3,
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