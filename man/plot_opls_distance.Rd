% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_opls_distance.R
\name{plot_opls_distance}
\alias{plot_opls_distance}
\title{plot_opls_distance}
\usage{
plot_opls_distance(
  dist_data,
  metadata = NULL,
  condition = NULL,
  significance = FALSE,
  significance_data = NULL,
  interactive = TRUE,
  title_label = "",
  x_title = "Score Component Distance",
  y_title = "Orthogonal Component Distance",
  marker_size = 6,
  title_label_size = 18,
  axis_title_size = 14,
  opls_cohort_text_format = "bold",
  opls_cohort_text_align = "right",
  opls_cohort_title_size = 18,
  opls_cohort_sample_size = 15,
  opls_plot_axis_format = "bold",
  opls_plot_axis_text_size = 14
)
}
\arguments{
\item{dist_data}{A dataframe containing the distance data.}

\item{metadata}{Dataframe containing samples to cohort mapping information}

\item{condition}{A metadata column where cohorts are present}

\item{significance}{A boolean indicating whether to plot significance level lines.}

\item{significance_data}{A numeric of length 2 containing the significance level distances of both axes.}

\item{interactive}{make plot interactive (default is TRUE)}

\item{title_label}{Title of the plot}

\item{x_title}{The title for x-axis}

\item{y_title}{The title for y-axis}

\item{marker_size}{The size of marker point}

\item{title_label_size}{Size of title label}

\item{axis_title_size}{Size of axis title}

\item{opls_cohort_text_format}{set text format of cohort legends}

\item{opls_cohort_text_align}{align cohort legends}

\item{opls_cohort_title_size}{set font size of cohort title}

\item{opls_cohort_sample_size}{set font size of cohorts}

\item{opls_plot_axis_format}{set axis format}

\item{opls_plot_axis_text_size}{set axis text size}
}
\value{
a plotly object
}
\description{
Plot the distance outliers plot.
}
\examples{
plot_opls_distance(dist_data,Metadata,"condition",significance=TRUE, significance_data)
}
