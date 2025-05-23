% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02.EMP.utilis.R, R/Plot_EMP_heatmap.R
\name{EMP_heatmap_plot}
\alias{EMP_heatmap_plot}
\alias{EMP_heatmap_plot,EMP_assay_heatmap_union-method}
\alias{EMP_heatmap_plot,EMP_cor_analysis-method}
\alias{EMP_heatmap_plot,EMP_WGCNA_cor_analysis-method}
\alias{EMP_heatmap_plot,EMP_WGCNA_cor_analysis2-method}
\alias{EMP_heatmap.EMP_cor_analysis}
\alias{EMP_heatmap.WGCNA}
\alias{EMP_heatmap.EMP_assay_data}
\title{EMP heatmap plot}
\usage{
EMP_heatmap_plot(obj, ...)

\S4method{EMP_heatmap_plot}{EMP_assay_heatmap_union}(obj, ...)

\S4method{EMP_heatmap_plot}{EMP_cor_analysis}(obj, ...)

\S4method{EMP_heatmap_plot}{EMP_WGCNA_cor_analysis}(obj, ...)

\S4method{EMP_heatmap_plot}{EMP_WGCNA_cor_analysis2}(obj, ...)

EMP_heatmap.EMP_cor_analysis(
  obj,
  palette = c("steelblue", "white", "darkred"),
  clust_row = FALSE,
  clust_col = FALSE,
  dist_method = "euclidean",
  clust_method = "complete",
  tree_size = 0.1,
  show = "all",
  label_size = 4,
  mytheme = "theme()"
)

EMP_heatmap.WGCNA(
  obj,
  palette = c("steelblue", "white", "darkred"),
  clust_row = FALSE,
  clust_col = FALSE,
  dist_method = "euclidean",
  clust_method = "complete",
  tree_size = 0.1,
  show = "all",
  label_size = 4,
  mytheme = "theme()"
)

EMP_heatmap.EMP_assay_data(
  obj,
  palette = c("steelblue", "white", "darkred"),
  rotate = FALSE,
  scale = NULL,
  bySample = "default",
  logbase = 2,
  pseudocount = 1e-07,
  clust_row = FALSE,
  clust_col = FALSE,
  dist_method = "euclidean",
  clust_method = "complete",
  tree_size = 0.1,
  label_size = 4,
  mytheme = "theme()"
)
}
\arguments{
\item{obj}{EMPT or EMP object}

\item{...}{...}

\item{palette}{1-3 character string. Color palette. (default: steelblue, white, darkred)}

\item{clust_row}{A boolean. Whether the function clust the row or not. (default:FALSE)}

\item{clust_col}{A boolean. Whether the function clust the row or not. (default:FALSE)}

\item{dist_method}{A character string. More see stats::dist. (default: euclidean)}

\item{clust_method}{A character string. More see fastcluster::hclust (default: complete)}

\item{tree_size}{A number between 0 and 1. Set the clust tree size. (default:0.1)}

\item{show}{A character string. Show inluding all, sig and pvalue.}

\item{label_size}{A number. Set the label size. (default:4)}

\item{mytheme}{Modify components of a theme according to the \code{\link[ggplot2]{theme}} and \code{\link[ggplot2]{ggtheme}}.}

\item{rotate}{A boolean. Whether rotate the heatmap or not. (Only activated for EMP_assay_data)}

\item{scale}{A character string,eg relative, clr, standardize, etc. The parameter works in the same way as the method in \code{\link[EasyMultiProfiler]{EMP_decostand}} (Only activated for EMP_assay_data)}

\item{bySample}{A boolean. Whether the function decostand by the sample or feature. Detaled information in the EMP_decostand. (Only activated for EMP_assay_data)}

\item{logbase}{An interger. The logarithm base used in method = "log".(default=2). Detaled information in the EMP_decostand. (Only activated for EMP_assay_data)}

\item{pseudocount}{A number. The logarithm pseudocount used in method = "clr" or "alr".(default=0.0000001). (Only activated for EMP_assay_data)}
}
\value{
EMPT or EMP object
}
\description{
EMP heatmap plot
}
\section{Detaild about Palettes}{

The following palettes are available for use with these scales:
\describe{
  BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral,
  Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2, Set3
  Blues, BuGn, BuPu, GnBu, Greens, Greys, Oranges, OrRd, PuBu, 
  PuBuGn, PuRd, Purples, RdPu, Reds, YlGn, YlGnBu, YlOrBr, YlOrRd
}

Modify the palette through the `palette` argument.
}

\section{Detaild about scale}{

When the scale parameter is enabled, data normalization only takes effect within this function and does not affect the original data in the EMPT. 
If you need to transform the data, consider using the EMP_EMP_decostand function in the workflow.
}

\examples{
\dontrun{
data(MAE)
## for assay
MAE |>
  EMP_assay_extract('geno_ec') |>
  EMP_diff_analysis(method='DESeq2',.formula = ~Group) |>
  EMP_filter(feature_condition = pvalue<0.05 & abs(fold_change) >3.5) |>
  EMP_decostand(method = 'clr') |>
  EMP_heatmap_plot(rotate=FALSE,palette='Spectral')

MAE |>
  EMP_assay_extract('geno_ec') |>
  EMP_diff_analysis(method='DESeq2',.formula = ~Group) |>
  EMP_filter(feature_condition = pvalue<0.05 & abs(fold_change) >3.5) |>
  EMP_collapse(estimate_group = 'Group',collapse_by = 'col') |> # collapse the data by group
  EMP_heatmap_plot(rotate=TRUE,palette='Spectral')

## for cor analysis
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

(k1 + k2) |> EMP_cor_analysis(method = 'spearman') |>
  EMP_heatmap_plot() ## Visualization

## for WGCNA
MAE |>
  EMP_assay_extract('geno_ec')  |> 
  EMP_identify_assay(method = 'edgeR',estimate_group = 'Group') |>
  EMP_WGCNA_cluster_analysis(RsquaredCut = 0.85,mergeCutHeight=0.4)  |>
  EMP_WGCNA_cor_analysis(coldata_to_assay = c('BMI','PHQ9','GAD7','HAMD','SAS','SDS'),method='spearman') |>
  EMP_heatmap_plot(palette = 'Spectral') 
}
}
