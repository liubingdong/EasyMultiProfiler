% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Support_EMP_import.R
\name{EMP_taxonomy_import}
\alias{EMP_taxonomy_import}
\title{Import microbial data into SummariseExperiment}
\usage{
EMP_taxonomy_import(
  file = NULL,
  data = NULL,
  humann_format = FALSE,
  file_format = NULL,
  start_level = "Kindom",
  assay_name = NULL,
  duplicate_feature = NULL,
  sep = if (humann_format == "TRUE") "|" else ";"
)
}
\arguments{
\item{file}{A file path.}

\item{data}{A dataframe.The row must be the feature and the column is the sample.}

\item{humann_format}{A boolean. Whether the function import the data in the Metaphlan or Humann format.}

\item{file_format}{A string including biom and gzv.}

\item{start_level}{A string. Specific the start level of input data from Domain,Kindom,Phylum,Class,Order,Family,Genus,Species,Strain.(default:Kindom)}

\item{assay_name}{A character string. Indicate what kind of result the data belongs to, such as counts, relative abundance, TPM, etc.}

\item{duplicate_feature}{A boolean. Whether the feature exist the dupicated name.}

\item{sep}{The field separator character. Values on feature column of the file are separated by this character. (defacult:';',when humann_format=T defacult:'|')}
}
\value{
SummmariseExperiment object
}
\description{
Import microbial data into SummariseExperiment
}
\details{
Paramter file_format and humann_format help the function import data properly. Data in humann format is usually is generated from Metaphlan and Humann. Data in biom format is usually is generated from Qiime1. Data in qzv format is usually is generated from Qiime2.
}
\examples{
# More examples and tutorial could be found on website: 
# https://liubingdong.github.io/EasyMultiProfiler_tutorial/
}
