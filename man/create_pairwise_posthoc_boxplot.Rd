% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_pairwise_posthoc_boxplot.R
\name{create_pairwise_posthoc_boxplot}
\alias{create_pairwise_posthoc_boxplot}
\title{create_pairwise_posthoc_boxplot}
\usage{
create_pairwise_posthoc_boxplot(
  posthoc_data_long_df = NULL,
  intensity_data = NULL,
  selected_metabolite = NULL,
  selected_interaction = NULL,
  selected_comparisons = NULL,
  filter_pvalues = FALSE,
  plot_axis_rotation = FALSE
)
}
\arguments{
\item{posthoc_data_long_df}{A dataframe containing long format post-hoc data.}

\item{intensity_data}{A dataframe containing intensity data.}

\item{selected_metabolite}{The selected metabolite for analysis.}

\item{selected_interaction}{The selected interaction for analysis.}

\item{selected_comparisons}{The selected comparisons for analysis.}

\item{filter_pvalues}{Logical indicating whether to filter comparisons with adj. p-value >= 0.05.}

\item{plot_axis_rotation}{Logical indicating whether to rotate axis labels and wrap text.}
}
\value{
A ggplot2 object representing the pairwise post-hoc boxplot.
}
\description{
Creates a pairwise post-hoc boxplot for selected metabolite and interaction.
}
