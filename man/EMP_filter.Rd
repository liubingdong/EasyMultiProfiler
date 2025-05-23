% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Support_EMP_filter.R
\name{EMP_filter}
\alias{EMP_filter}
\title{Filer experssion or abundance data that match a condition}
\usage{
EMP_filter(
  obj,
  sample_condition,
  feature_condition,
  filterSample = NULL,
  filterFeature = NULL,
  experiment = NULL,
  show_info = NULL,
  action = "select",
  keep_result = FALSE,
  use_cached = TRUE
)
}
\arguments{
\item{obj}{EMPT object.}

\item{sample_condition}{Expressions that return a logical value, and are defined in terms of the variables in coldata. If multiple expressions are included, they are combined with the &，| operator.}

\item{feature_condition}{Expressions that return a logical value, and are defined in terms of the variables in rowdata. If multiple expressions are included, they are combined with the &，| operator.}

\item{filterSample}{A series of character strings. Select samples in the data exactly.}

\item{filterFeature}{A series of character strings. Select samples in the data exactly.}

\item{experiment}{A character string. Experiment name in the MultiAssayExperiment object.}

\item{show_info}{A character string. Set the class of EMPT to show properly.}

\item{action}{A character string. You can use the filterSample and filterFeature parameters in conjunction with this. The choice is whether to keep filterSample and filterFeature (select), or simply exclude them (kick).}

\item{keep_result}{If the input is TRUE, it means to keep all analysis results,regardless of how samples and features change. If the input is a name, it means to keep the corresponding analysis results.}

\item{use_cached}{A boolean. Whether the function use the results in cache or re-compute.}
}
\value{
EMPT object
}
\description{
Filer experssion or abundance data that match a condition
}
\examples{
data(MAE)
## from MultiAssayExperiment
MAE |>
  EMP_filter(sample_condition = BMI>20 & Sex == 'M') |>
  EMP_summary()

## from EMPT
MAE |>
  EMP_assay_extract(experiment = 'host_gene') |>
  EMP_identify_assay(method = 'edgeR',estimate_group = 'Group') |>
  EMP_diff_analysis(method='DESeq2',.formula = ~Group) |>
  EMP_filter(feature_condition = fdr < 0.05)
  
MAE |>
  EMP_assay_extract(experiment = 'taxonomy') |>
  EMP_alpha_analysis() |>
  EMP_filter(sample_condition = shannon >3 | invsimpson >19)

## Precise selection
MAE |>
  EMP_assay_extract(experiment = 'taxonomy') |>
  EMP_alpha_analysis() |>
  EMP_filter(sample_condition = shannon >3 | invsimpson >19,
             filterSample = c('P11774','P31579'),action = 'kick') # Accurately kick samples based on satisfying sample_condition
}
