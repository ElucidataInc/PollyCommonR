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
#' @param x_label Label x-axis
#' @param y_label Label y-axis
#' @param title_label Title of the plot
#' @param marker_size Size of marker point
#' @param plot_id source id for the plotly plot
#' @param interactive make plot interactive (default is TRUE)
#' @return plotly or ggplot object
#' @examples
#' plot_volcano_from_limma(diff_exp, log2fc_range = 0.5, p_val_cutoff = 0.05, p_val_type = "P.Value")
#' @import ggplot2 plotly ggsci ggrepel latex2exp
#' @export
plot_volcano_from_limma <- function(diff_exp = NULL, log2fc_range = NULL, p_val_cutoff = NULL, 
                                    p_val_type = "P.Value", annotate_id = NULL, row_desc = NULL, 
                                    annotate_col = NULL, text_hover_col = NULL, category_col = NULL,
                                    x_label = NULL, y_label = NULL, title_label = NULL, marker_size = 8, 
                                    plot_id = NULL, interactive = TRUE) {
  message("Make Volcano Plot Started...")
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
  
  diff_exp <- diff_exp[!is.infinite(rowSums(diff_exp[, required_cols])), ]
  diff_exp <- diff_exp[!apply(diff_exp[, required_cols], 1, anyNA), ]
  
  if (nrow(diff_exp) < 1){
    warning("Differential Expression dataframe has zero valid rows")
    return (NULL)
  }
  
  diff_exp <- data.frame(id = row.names(diff_exp), diff_exp, stringsAsFactors = FALSE, check.names = FALSE)
  diff_exp$threshold <- "not significant"  
  diff_exp[(abs(diff_exp$logFC) >= log2fc_range) & (diff_exp[, p_val_type] <= p_val_cutoff), "threshold"] <- "significant"
  significance_table <- base::table(diff_exp$threshold)
  print (significance_table)
  
  if (identical(row_desc, NULL) | identical(category_col, NULL)){
    for (sig_type in unique(diff_exp$threshold)){
      diff_exp[diff_exp$threshold == sig_type, "threshold"] <- paste0(sig_type, " (", significance_table[[sig_type]], ")") 
    }
  }
  
  diff_exp$text_hover <- row.names(diff_exp)
  diff_exp$category_sym <- NULL  
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
  y_col = gsub("pval", p_val_type, "-log10(pval)")  
  x_val <- diff_exp$logFC
  y_val <- -log10(diff_exp[, p_val_type])                                               
  
  if (interactive == TRUE){
    if (identical(x_label, NULL)){ x_label <- xaxis_lab_pl }
    if (identical(y_label, NULL)){ y_label <- yaxis_lab_pl }

    filtered_diff_exp <- diff_exp[row.names(diff_exp) %in% annotate_id, ]
    if(nrow(filtered_diff_exp) == 0){
      a <- NULL
    } else {
      annotate_text <- row.names(filtered_diff_exp)  
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
    p <- plotly::plot_ly(source = plot_id) %>%
      add_trace(
        x = x_val, y = y_val, customdata = diff_exp$id,
        type = "scattergl", mode = "markers",
        marker = list(size = marker_size),
        color = diff_exp$threshold,
        colors = significance_color,
        symbol = diff_exp$category_sym,
        text = diff_exp$text_hover  
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
                     modeBarButtons = list(list("zoomIn2d"), 
                                           list("zoomOut2d"), 
                                           list('toImage')), 
                     mathjax = 'cdn')
  } else {
    if (identical(annotate_col, NULL)){ annotate_col <- "id" }
    if (identical(x_label, NULL)){ x_label <- xaxis_lab_gg }
    if (identical(y_label, NULL)){ y_label <- yaxis_lab_gg }      
    
    p <- ggplot(diff_exp, aes_string(x = x_col, y = y_col, color = "threshold", fill = "threshold", shape = category_col), text = id) + 
      geom_point( size = marker_size/2, alpha = 0.7) + # scatter plot function with shape of points defined as 21 scale.   
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