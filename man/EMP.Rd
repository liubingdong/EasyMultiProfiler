\docType{package}
\name{EasyMultiPlot}
\alias{EMP-package}
\title{EasyMultiPlot: An Efficient and Convenient R package in Multi-omics Down-Stream Analysis and Visualization for Clinical Study}

\section{Description}{
\if{html}{\figure{logo.png}{options: align='right' alt='logo' width='120'}}

The EasyMulitiPlot package aims to offer a user-friendly and efficient multi-omics data analysis tool on the R platform. It facilitates various essential tasks related to microbiome, genome, and metabolite downstream analysis, providing a seamless workflow from beginning to end.
}

\section{Function}{
{Data Extract} : 
\itemize{
	\item \code{\link{EMP_assay_extract}}  
 	\item \code{\link{EMP_coldata_extract}}
 	\item \code{\link{EMP_rowdata_extract}}
}
{Data Preparetion} : 
\itemize{
	\item \code{\link{EMP_impute}} 
	\item \code{\link{EMP_rrarefy}}  
 	\item \code{\link{EMP_decostand}}
 	\item \code{\link{EMP_collapse}}
 	\item \code{\link{EMP_adjust_abudance}}
}
{Data Analysis} : 
\itemize{
	\item \code{\link{EMP_alpha_analysis}}  
	\item \code{\link{EMP_cluster_analysis}}
	\item \code{\link{EMP_cor_analysis}}
	\item \code{\link{EMP_diff_analysis}}
	\item \code{\link{EMP_dimension_analysis}}
	\item \code{\link{EMP_cor_analysis}}
	\item \code{\link{EMP_enrich_analysis}}
	\item \code{\link{EMP_GSEA}}
	\item \code{\link{EMP_WGCNA_cluster_analysis}}
	\item \code{\link{EMP_WGCNA_cor_analysis}}
	\item \code{\link{EMP_multienrich_analysis}}
 }
{Data Virtualization} : 
\itemize{
	\item \code{\link{EMP_boxplot}} 
	\item \code{\link{EMP_curveplot}}  
 	\item \code{\link{EMP_heatmap_plot}}
 	\item \code{\link{EMP_volcanol_plot}}
 }
{Data support} : 
\itemize{
	\item \code{\link{EMP_filter}} 
	\item \code{\link{EMP_summary}}
	\item \code{\link{EMP_feature_convert}}    
	\item \code{\link{EMP_taxonomy_import}} 
	\item \code{\link{EMP_function_import}} 
	\item \code{\link{EMP_normal_import}} 
 }
}

\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/XXXXX/}
}
}

