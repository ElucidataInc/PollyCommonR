#' volcano_plot
#'
#' makes a plotly volcano plot
#'
#' @param diff_exp_rdesc The differential expression dataframe from limma
#' @param log2fc_range Vector of log2FC range
#' @param p_val_cutoff The pval cutoff
#' @param fdr_cutoff The FDR cutoff
#' @return plotly object
#' @examples
#' plot_volcano_from_limma(diff_exp_rdesc, log2fc_range = 0.5, p_val_cutoff = 0.05, fdr_cutoff = NULL)
#' @export
plot_volcano_from_limma <- function(diff_exp_rdesc = NULL, log2fc_range = NULL, p_val_cutoff = NULL, fdr_cutoff = NULL) {

  message("Make Volcano Plot Started...")
  require(plotly)

  if(identical(fdr_cutoff, NULL)){
    diff_exp_rdesc <- diff_exp_rdesc[!is.infinite(rowSums(diff_exp_rdesc)),]
    diff_exp_rdesc$threshold <- "not significant"
    diff_exp_rdesc[(abs(diff_exp_rdesc$logFC) >= log2fc_range) &
      (diff_exp_rdesc[, "P.Value"] <= p_val_cutoff), "threshold"] <- "significant"
    print(table(diff_exp_rdesc$threshold))
    pal <- c("grey", "red")
    p <- plot_ly() %>%
      add_trace(
        x = diff_exp_rdesc$logFC, y = -log10(diff_exp_rdesc$P.Value),
        type = "scattergl", mode = "markers",
        color = diff_exp_rdesc$threshold,
        colors = pal,
        text = rownames(diff_exp_rdesc)
      ) %>%
      layout(
        title = "Volcano plot",
        yaxis = list(title = "Negative log10 p-value"),
        xaxis = list(title = "log2 Fold Change")
      ) %>% plotly::config(displaylogo = FALSE,
        modeBarButtons = list(list("toImage"), list("zoomIn2d"), list("zoomOut2d")))

    message("Make Volcano Plot Completed...")

    return(p)
  } else{
    diff_exp_rdesc <- diff_exp_rdesc[!is.infinite(rowSums(diff_exp_rdesc)),]
    diff_exp_rdesc$threshold <- "not significant"
    diff_exp_rdesc[(abs(diff_exp_rdesc$logFC) >= log2fc_range) &
      (diff_exp_rdesc[, "adj.P.Val"] <= fdr_cutoff), "threshold"] <- "significant"
    print(table(diff_exp_rdesc$threshold))
    pal <- c("grey", "red")
    p <- plot_ly() %>%
      add_trace(
        x = diff_exp_rdesc$logFC, y = -log10(diff_exp_rdesc$adj.P.Val),
        type = "scattergl", mode = "markers",
        color = diff_exp_rdesc$threshold,
        colors = pal,
        text = rownames(diff_exp_rdesc)
      ) %>%
      layout(
        title = "Volcano plot",
        yaxis = list(title = "Negative log10 FDR"),
        xaxis = list(title = "log2 Fold Change")
      ) %>% plotly::config(displaylogo = FALSE,
        modeBarButtons = list(list("toImage"), list("zoomIn2d"), list("zoomOut2d")))

    message("Make Volcano Plot Completed...")

    return(p)
  }

}
