% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_intomix_input.R
\name{make_intomix_input}
\alias{make_intomix_input}
\title{make_intomix_input}
\usage{
make_intomix_input(
  norm_data = NULL,
  metadata = NULL,
  cohorts = NULL,
  cohorts_compare_data = NULL,
  cohort_col = "Cohort",
  rownames_col = "uniqueId",
  remove_duplicates = TRUE,
  id_col = "compoundId",
  pval_type = "P.Value",
  pval_adjust_method = "BH",
  log_flag = TRUE,
  drop_na = FALSE,
  replace_na_with_zero = FALSE
)
}
\arguments{
\item{norm_data}{A sample intenisty data with samples in columns.}

\item{metadata}{dataframe containing metadata information}

\item{cohorts}{A vector of cohorts used for cohorts comparisons}

\item{cohorts_compare_data}{A dataframe having predefined cohorts comparisons}

\item{cohort_col}{A metadata column where cohorts are present}

\item{rownames_col}{A raw_intensity_df column which is to be used for assigning rownames.}

\item{remove_duplicates}{Remove duplicate ids rows based on pval}

\item{id_col}{A column id which will be used to remove duplicates}

\item{pval_type}{Select p value type from P.Value or adj.P.Val}

\item{pval_adjust_method}{Provide pval adjust method ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). Also check p.adjust from stats.}

\item{log_flag}{Check if data is log transformed (TRUE) or not (FALSE). If not (FALSE) then do internally log2 transformation.}

\item{drop_na}{Drop rows having NA's used in sample_intensity_matrix function.}

\item{replace_na_with_zero}{Replace all NA's with zero used in sample_intensity_matrix function.}
}
\value{
A dataframe of differential expression for all cohorts comparison.
}
\description{
Create an intomix input file
}
\examples{
make_intomix_input(norm_data = NULL, metadata = NULL, cohorts = NULL,  
                   cohorts_compare_data = NULL, cohort_col = 'Cohort', 
                   rownames_col = 'uniqueId', remove_duplicates = TRUE, 
                   id_col = 'compoundId', pval_adjust_method = 'BH',
                   log_flag = TRUE )
}
