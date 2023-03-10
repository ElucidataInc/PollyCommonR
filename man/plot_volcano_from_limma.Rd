% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_volcano_from_limma.R
\name{plot_volcano_from_limma}
\alias{plot_volcano_from_limma}
\title{plot_volcano_from_limma}
\usage{
plot_volcano_from_limma(
  diff_exp = NULL,
  log2fc_range = 1,
  p_val_cutoff = 0.05,
  p_val_type = "P.Value",
  annotate_id = NULL,
  row_desc = NULL,
  annotate_col = NULL,
  text_hover_col = NULL,
  category_col = NULL,
  marker_size_by_expr = FALSE,
  marker_expr_col = NULL,
  marker_size_range = c(5, 25),
  marker_size = 8,
  marker_opacity = 0.5,
  x_label = NULL,
  y_label = NULL,
  title_label = NULL,
  plot_id = NULL,
  plotly_highlight = FALSE,
  highlight_on = "plotly_click",
  highlight_off = "plotly_doubleclick",
  highlight_persistent = FALSE,
  highlight_color = "blue",
  highlight_opacitydim = 0.8,
  highlight_debounce = 0,
  interactive = TRUE
)
}
\arguments{
\item{diff_exp}{The differential expression dataframe from limma}

\item{log2fc_range}{The absolute log2FC value to set threshold}

\item{p_val_cutoff}{The pval cutoff}

\item{p_val_type}{The pval type (P.Value or adj.P.Val)}

\item{annotate_id}{A vector of ids (rownames of row_desc) to be annotated}

\item{row_desc}{A dataframe of row descriptors of the matrix}

\item{annotate_col}{A row descriptor column which is used to show names for annotated ids on plot}

\item{text_hover_col}{A row descriptor column which is used to hover text on plot}

\item{category_col}{A row descriptor column which is used to add shapes to different categories}

\item{marker_size_by_expr}{Show size of markers by expression values, TRUE/FALSE}

\item{marker_expr_col}{A numeric expression column present in diff_exp}

\item{marker_size_range}{A numeric vector of minimum and maximum values of size of markers}

\item{marker_size}{Size of marker point}

\item{marker_opacity}{The opacity of the markers}

\item{x_label}{Label x-axis}

\item{y_label}{Label y-axis}

\item{title_label}{Title of the plot}

\item{plot_id}{source id for the plotly plot}

\item{plotly_highlight}{Highlight points on the plotly plot}

\item{highlight_on}{Highlight points when clicking the points (plotly_click), hovering on points (plotly_hover) or selecting multiple points (plotly_selected)}

\item{highlight_off}{Switch off highlighting points when plotly_doubleclick, plotly_deselect or plotly_relayout}

\item{highlight_persistent}{should selections persist (i.e., accumulate), default is TRUE}

\item{highlight_color}{The color of the highlighted points}

\item{highlight_opacitydim}{A number between 0 and 1 used to reduce the opacity of non-selected traces (by multiplying with the existing opacity)}

\item{highlight_debounce}{The amount of time to wait before firing an event (in milliseconds). The default of 0 means do not debounce at all. The color of the highlighted points}

\item{interactive}{make plot interactive (default is TRUE)}
}
\value{
plotly or ggplot object
}
\description{
makes a volcano plot
}
\examples{
plot_volcano_from_limma(diff_exp, log2fc_range = 0.5, p_val_cutoff = 0.05, p_val_type = "P.Value")
}
