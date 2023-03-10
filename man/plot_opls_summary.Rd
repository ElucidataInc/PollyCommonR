% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_opls_summary.R
\name{plot_opls_summary}
\alias{plot_opls_summary}
\title{plot_opls_summary}
\usage{
plot_opls_summary(
  opl,
  title_label = "Model Summary Diagnostics",
  title_label_size = 18,
  axis_title_size = 14,
  opls_cohort_text_format = "bold",
  opls_legend_text_align = "right",
  opls_legend_title_size = 18,
  opls_plot_axis_format = "bold",
  opls_plot_axis_text_size = 14
)
}
\arguments{
\item{opl}{An opls object}

\item{title_label}{Title of the plot}

\item{title_label_size}{Size of title label}

\item{axis_title_size}{Size of axis title}

\item{opls_cohort_text_format}{set text format of legends}

\item{opls_legend_text_align}{align legends}

\item{opls_legend_title_size}{set font size of cohort title}

\item{opls_plot_axis_format}{set axis format}

\item{opls_plot_axis_text_size}{set axis text size}
}
\value{
a ggplot2 object
}
\description{
Plots the summary statistics of the opls object.
}
\examples{
plot_opls_summary(OPLS_object,title_label="My Graph")
}
