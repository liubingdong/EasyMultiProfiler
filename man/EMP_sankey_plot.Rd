% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02.EMP.utilis.R, R/Plot_EMP_sankey.R
\name{EMP_sankey_plot}
\alias{EMP_sankey_plot}
\alias{EMP_sankey_plot,EMP_cor_analysis-method}
\alias{EMP_sankey_plot.EMP_cor_analysis}
\title{EMP_sankey_plot}
\usage{
EMP_sankey_plot(obj, ...)

\S4method{EMP_sankey_plot}{EMP_cor_analysis}(obj, ...)

EMP_sankey_plot.EMP_cor_analysis(
  obj,
  positive_col = "#CC79A7",
  negtive_col = "steelblue",
  fontSize = 7,
  nodeWidth = 15,
  nodePadding = 10,
  height = NULL,
  width = NULL,
  palette = c("#009E73", "#F0E442", "#CB563f", "#6c4f9d", "#4DBBD5FF", "#00A087FF",
    "#b66e1e", "#8491B4FF"),
  ...
)
}
\arguments{
\item{obj}{EMPT}

\item{...}{Additional parameters, see also \code{\link[networkD3]{sankeyNetwork}}}

\item{positive_col}{Edge colour for positive relation. [Default:#CC79A7]}

\item{negtive_col}{Edge colour for negitive relation. [Default:steelblue]}

\item{fontSize}{Numeric font size in pixels for the node text labels.}

\item{nodeWidth}{Numeric width of each node.}

\item{nodePadding}{Numeric essentially influences the width height.}

\item{height}{Numeric height for the network graph's frame area in pixels.}

\item{width}{Numeric width for the network graph's frame area in pixels.}

\item{palette}{Colour palette for nodes.}
}
\value{
EMP object
}
\description{
EMP_sankey_plot
}
\examples{
data(MAE)
k1 <- MAE |>
  EMP_assay_extract('taxonomy') |>
  EMP_collapse(estimate_group = 'Genus',collapse_by = 'row') |>
  EMP_diff_analysis(method='DESeq2', .formula = ~Group) |>
  EMP_filter(feature_condition = pvalue<0.05)

k2 <- MAE |>
  EMP_collapse(experiment = 'untarget_metabol',na_string=c('NA','null','','-'),
               estimate_group = 'MS2kegg',method = 'sum',collapse_by = 'row') |>
  EMP_diff_analysis(method='DESeq2', .formula = ~Group) |>
  EMP_filter(feature_condition = pvalue<0.05 & abs(fold_change) > 1.5)

k3 <- MAE |>
  EMP_assay_extract('geno_ec') |>
  EMP_diff_analysis(method='DESeq2', .formula = ~Group) |>
  EMP_filter(feature_condition = pvalue<0.05 & abs(fold_change) > 2)

(k1 + k3 + k2) |> EMP_cor_analysis() |>
  EMP_sankey_plot(height=600,width=700,fontSize=10,nodeWidth=15,nodePadding=5) # more parameters
}
