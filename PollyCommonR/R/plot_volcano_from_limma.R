#' volcano_plot
#'
#' makes a plotly volcano plot
#'
#' @param diff_exp_rdesc The differential expression dataframe from limma
#' @param log2fc_range Vector of log2FC range
#' @param p_val_cutoff The pval cutoff
#' @param fdr_cutoff The FDR cutoff
#' @param annotate_id A vector of genes/metabolites to be annotated
#' @param title_label Title of the plot
#' @param interactive make plot interactive (default is TRUE)
#' @return plotly object
#' @examples
#' plot_volcano_from_limma(diff_exp_rdesc, log2fc_range = 0.5, p_val_cutoff = 0.05, fdr_cutoff = NULL,
#'  annotate_id = c('a','b'), interactive = TRUE)
#' @import ggplot2 plotly ggsci ggrepel latex2exp
#' @export
plot_volcano_from_limma <- function(diff_exp_rdesc = NULL, log2fc_range = NULL, p_val_cutoff = NULL, 
                                    fdr_cutoff = NULL, annotate_id = NULL, 
                                    title_label = "", interactive = TRUE) {
  message("Make Volcano Plot Started...")
  require(ggplot2)
  require(plotly)
  require(ggsci)
  require(ggrepel)
  require(latex2exp)
  
  diff_exp_rdesc <- diff_exp_rdesc[!is.infinite(rowSums(diff_exp_rdesc)),]
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
  
  if (interactive == TRUE){
    filtered_diff_exp <- diff_exp_rdesc[row.names(diff_exp_rdesc) %in% annotate_id, ]
    if(nrow(filtered_diff_exp) == 0){
      a <- NULL
    }else{
      a <- list(
        x = filtered_diff_exp$logFC,
        y = -log10(filtered_diff_exp[[pval_type]]),
        text = rownames(filtered_diff_exp),
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
        color = diff_exp_rdesc$threshold,
        colors = pal,
        text = rownames(diff_exp_rdesc)
      ) %>%
      layout(
        title = title_label,
        yaxis = list(title = yaxis_lab_pl),
        xaxis = list(title = xaxis_lab_pl),
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
  }else{
    p <- ggplot(diff_exp_rdesc, aes_string(x = x_col, y = y_col, fill = 'threshold')
                , text = id) + # calls the ggplot function with dose on the x-axis and len on the y-axis 
      geom_point(shape = 21, size = 4, alpha = 0.7) + # scatter plot function with shape of points defined as 21 scale.
      ggtitle(title_label)+
      labs(x = xaxis_lab_gg, y = yaxis_lab_gg,  fill = "Significance") + # x and y axis labels
      ggsci::scale_color_aaas() + # filling the point colors
      theme(legend.position = "right", legend.direction = "vertical", # legend positioned at the bottom, horizantal direction,
            axis.line = element_line(size=1, colour = "black"),	# axis line of size 1 inch in black color
            panel.grid.major = element_blank(),	# major grids included
            panel.grid.minor = element_blank(),	# no minor grids
            panel.border = element_blank(), panel.background = element_blank(), # no borders and background color
            plot.title = element_text(colour="black", size = 18, face = "plain", hjust=0.5),
            axis.title = element_text(colour="black", size = 15, face = "bold"), # axis title 
            axis.text.x = element_text(colour="black", size = 10, margin=unit(c(0.5,0.5,0.1,0.1), "cm"), face = "bold"), # x-axis text in fontsize 10
            axis.text.y = element_text(colour="black", size = 10, margin=unit(c(0.5,0.5,0.1,0.1), "cm"), face = "bold"), # y-axis text in fontsize 10
            legend.text = element_text(size = 10, face = "bold"),
            legend.title = element_text(colour="black", size=12, face="bold"),
            axis.ticks.length = unit(-0.25, "cm")) # ticks facing inward with 0.25cm length
    p <- p + ggrepel::geom_text_repel(data = subset(diff_exp_rdesc, id %in% annotate_id), aes(label = id),size = 5,
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
