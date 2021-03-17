#' volcano_plot
#'
#' makes a plotly volcano plot
#'
#' @param diff_exp_rdesc The differential expression dataframe from limma
#' @param log2fc_range Vector of log2FC range
#' @param p_val_cutoff The pval cutoff
#' @param fdr_cutoff The FDR cutoff
#' @param annotate_id A vector of ids (rownames of row_desc) to be annotated
#' @param row_desc A dataframe of row descriptors of the matrix
#' @param annotate_col A row descriptor column which is used to show names for annotated ids on plot
#' @param text_hover_col A row descriptor column which is used to hover text on plot
#' @param category_col A row descriptor column which is used to add shapes to different categories
#' @param x_label Label x-axis
#' @param y_label Label y-axis
#' @param title_label Title of the plot
#' @param marker_size Size of marker point
#' @param interactive make plot interactive (default is TRUE)
#' @return plotly object
#' @examples
#' plot_volcano_from_limma(diff_exp_rdesc, log2fc_range = 0.5, p_val_cutoff = 0.05, fdr_cutoff = NULL,
#'  annotate_id = c('a','b'), interactive = TRUE)
#' @import ggplot2 plotly ggsci ggrepel latex2exp
#' @export
plot_volcano_from_limma <- function(diff_exp_rdesc = NULL, log2fc_range = NULL, p_val_cutoff = NULL, 
                                    fdr_cutoff = NULL, annotate_id = NULL, row_desc = NULL, 
                                    annotate_col = NULL, text_hover_col = NULL, category_col = NULL,
                                    x_label = NULL, y_label = NULL, title_label = NULL, marker_size = 8, 
                                    interactive = TRUE) {
  message("Make Volcano Plot Started...")
  require(ggplot2)
  require(plotly)
  require(ggsci)
  require(ggrepel)
  require(latex2exp)
  
  if (identical(diff_exp_rdesc, NULL)){
    warning("Differential Expression dataframe is NULL")
    return (NULL)
  }
  
  if (!identical(row_desc, NULL)){
    if (!identical(class(row_desc), "data.frame")){
      warning("The row_desc is not a dataframe, please provide valid row_desc")
      return(NULL) 
    }   
    
    if (!all(row.names(diff_exp_rdesc) %in% row.names(row_desc))){
      warning("Not all rownames of diff_exp_rdesc are present in rownames of row_desc")
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
    
    diff_exp_rdesc <- base::transform(base::merge(diff_exp_rdesc,row_desc,by=0), row.names=Row.names, Row.names=NULL) 
  }
  else {
    annotate_col <- "id"
    text_hover_col <- NULL
    category_col <- NULL    
  }  
  
  if (!identical(annotate_id, NULL)){
    common_annotate_id <- base::intersect(annotate_id, row.names(diff_exp_rdesc))
    if (length(common_annotate_id) < 1){
      warning("The annotate ids are not matched with rownames of diff_exp_rdesc")    
    }
    else {
      diff_annotate_id <- base::setdiff(annotate_id, row.names(diff_exp_rdesc))
      if (length(diff_annotate_id) >= 1){
        warning(paste("The following annotate ids are not matched with rownames of diff_exp_rdesc :", paste(sQuote(diff_annotate_id), collapse = ", "), collapse = " ")) 
      }  
    }       
  }
  
  diff_exp_rdesc <- diff_exp_rdesc[!is.infinite(rowSums(diff_exp_rdesc[, c("logFC", "P.Value", "adj.P.Val")])), ]
  diff_exp_rdesc <- diff_exp_rdesc[!apply(diff_exp_rdesc[, c("logFC","P.Value", "adj.P.Val")], 1, anyNA), ]
  
  if (nrow(diff_exp_rdesc) < 1){
    warning("Differential Expression dataframe has zero valid rows")
    return (NULL)
  }
  
  diff_exp_rdesc$threshold <- "not significant"
  diff_exp_rdesc$id <- row.names(diff_exp_rdesc)
  pal <- c("grey", "red")
  
  if(identical(fdr_cutoff, NULL)){
    diff_exp_rdesc[(abs(diff_exp_rdesc$logFC) >= log2fc_range) &
                     (diff_exp_rdesc[, "P.Value"] <= p_val_cutoff), "threshold"] <- "significant"
    print(table(diff_exp_rdesc$threshold))
    
    xaxis_lab_gg <- latex2exp::TeX("$\\log_{2}(fold \\, change)$")
    yaxis_lab_gg <- latex2exp::TeX("$-\\log_{10}(p \\, value)$")
    
    xaxis_lab_pl <- plotly::TeX("\\log_{2}(\\text{fold change})")
    yaxis_lab_pl <- plotly::TeX("-\\log_{10}(\\text{p value})")
    
    x_col = "logFC"
    y_col = "-log10(P.Value)"
    
    x_val <- diff_exp_rdesc$logFC
    y_val <- -log10(diff_exp_rdesc$P.Value)
    
    pval_type <- "P.Value"
    
  }else{
    diff_exp_rdesc[(abs(diff_exp_rdesc$logFC) >= log2fc_range) &
                     (diff_exp_rdesc[, "adj.P.Val"] <= fdr_cutoff), "threshold"] <- "significant"
    print(table(diff_exp_rdesc$threshold))
    
    xaxis_lab_gg <- latex2exp::TeX("$\\log_{2}(fold \\, change)$")
    yaxis_lab_gg <- latex2exp::TeX("$-\\log_{10}(adjusted \\,p \\, value)$")
    
    xaxis_lab_pl <- plotly::TeX("\\log_{2}(\\text{fold change})")
    yaxis_lab_pl <- plotly::TeX("-\\log_{10}(\\text{adjusted p value})")
    
    x_col <- "logFC"
    y_col <- "-log10(adj.P.Val)"
    
    x_val <- diff_exp_rdesc$logFC
    y_val <- -log10(diff_exp_rdesc$adj.P.Val)
    
    pval_type <- "adj.P.Val"
    
  }
  
  diff_exp_rdesc$text_hover <- row.names(diff_exp_rdesc)
  diff_exp_rdesc$category_sym <- NULL  
  if (!identical(row_desc, NULL)){     
    if (!identical(text_hover_col, NULL)){
      diff_exp_rdesc$text_hover <- diff_exp_rdesc[[text_hover_col]]
    } 
    if (!identical(category_col, NULL)){
      diff_exp_rdesc[[category_col]][unlist(lapply(diff_exp_rdesc[[category_col]], function(x) x %in% c(NA, "NA", "")))] <- "NA"
      diff_exp_rdesc$category_sym <- diff_exp_rdesc[[category_col]]                                             
    }                                            
  }
  
  if (interactive == TRUE){
    if (identical(x_label, NULL)){ x_label <- xaxis_lab_pl }
    if (identical(y_label, NULL)){ y_label <- yaxis_lab_pl }
    diff_exp_rdesc$text_hover <- row.names(diff_exp_rdesc)
    diff_exp_rdesc$category_sym <- NULL  
    if (!identical(row_desc, NULL)){     
      if (!identical(text_hover_col, NULL)){
        diff_exp_rdesc$text_hover <- diff_exp_rdesc[[text_hover_col]]
      } 
      if (!identical(category_col, NULL)){
        diff_exp_rdesc[[category_col]][unlist(lapply(diff_exp_rdesc[[category_col]], function(x) x %in% c(NA, "NA", "")))] <- "NA"
        diff_exp_rdesc$category_sym <- diff_exp_rdesc[[category_col]]                                             
      }                                            
    }      
    
    filtered_diff_exp <- diff_exp_rdesc[row.names(diff_exp_rdesc) %in% annotate_id, ]
    if(nrow(filtered_diff_exp) == 0){
      a <- NULL
    } else {
      annotate_text <- row.names(filtered_diff_exp)  
      if (!identical(row_desc, NULL)){  
        if (!identical(annotate_col, NULL)){ annotate_text <- filtered_diff_exp[[annotate_col]] }
      }
      
      a <- list(
        x = filtered_diff_exp$logFC,
        y = -log10(filtered_diff_exp[[pval_type]]),
        text = annotate_text,
        xref = "x",
        yref = "y",
        showarrow = TRUE,
        arrowhead = 3.5,
        ax = 20,
        ay = -40
      )}
    p <- plotly::plot_ly() %>%
      add_trace(
        x = x_val, y = y_val,
        type = "scattergl", mode = "markers",
        marker = list(size = marker_size),
        color = diff_exp_rdesc$threshold,
        colors = pal,
        symbol = diff_exp_rdesc$category_sym,
        text = diff_exp_rdesc$text_hover  
      ) %>%
      
      layout(
        title = title_label,
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
    
    p <- ggplot(diff_exp_rdesc, aes_string(x = x_col, y = y_col, color = "threshold", fill = "threshold", shape = category_col), text = id) + 
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
    p <- p + ggrepel::geom_text_repel(data = subset(diff_exp_rdesc, id %in% annotate_id), aes_string(label = annotate_col), color = "black", size = 3,
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