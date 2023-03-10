% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_samplewise_barplot.R
\name{create_samplewise_barplot}
\alias{create_samplewise_barplot}
\title{create_samplewise_barplot}
\usage{
create_samplewise_barplot(
  sample_raw_mat = NULL,
  metadata_df = NULL,
  id_name = NULL,
  cohorts_order = NULL,
  cohort_col = "Cohort",
  interactive = FALSE,
  x_label = "Sample",
  y_label = "Intensity",
  title_label = ""
)
}
\arguments{
\item{sample_raw_mat}{sample_raw_mat matrix/dataframe containing raw values}

\item{metadata_df}{A dataframe with samples to cohort mapping}

\item{id_name}{An id for which barplot to be made}

\item{cohorts_order}{The order of cohorts}

\item{cohort_col}{A cohort column present in metadata}

\item{interactive}{Make plot interactive using plotly}

\item{x_label}{Label x-axis}

\item{y_label}{Label y-axis}

\item{title_label}{Title of the plot}
}
\value{
ggplot object or plotly object
}
\description{
Makes samplewise barplot
}
\examples{
create_samplewise_barplot(sample_raw_mat = NULL, metadata_df = NULL, id_name = NULL, 
                          cohorts_order = NULL, cohort_col='Cohort', 
                          x_label = "", y_label = "", title_label = "")
}
