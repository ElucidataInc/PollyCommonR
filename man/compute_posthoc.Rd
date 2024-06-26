% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_posthoc.R
\name{compute_posthoc}
\alias{compute_posthoc}
\title{compute_posthoc}
\usage{
compute_posthoc(
  norm_data = NULL,
  metadata = NULL,
  cohort_col = "Cohort",
  post_hoc_test = "bonferroni",
  significant_ids = NULL
)
}
\arguments{
\item{norm_data}{The dataframe/matrix with samples in columns and features in rows.}

\item{metadata}{The dataframe with samples to cohort mapping.}

\item{cohort_col}{A vector of metadata columns used for grouping.}

\item{post_hoc_test}{The method used for adjusting p-values, options include "bonferroni", "fdr" (False Discovery Rate), and "holm".}
}
\value{
A list containing:
  - posthoc_data: A dataframe with adjusted p-values for each pairwise comparison.
  - posthoc_data_long_df: A long-format dataframe containing detailed results of pairwise comparisons.
  - posthoc_data_heatmap: A custom wide-format dataframe suitable for visualization.
}
\description{
Performs pairwise t-tests and adjusts p-values using specified methods.
}
\examples{
compute_posthoc(norm_data, metadata, cohort_col, post_hoc_test = "bonferroni")
}
