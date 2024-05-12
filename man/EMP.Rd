\docType{package}
\name{EasyMultiProfiler}
\alias{EasyMultiProfiler}
\title{EasyMultiProfiler: An Efficient and Convenient R package in Multi-omics Down-Stream Analysis and Visualization}

\description{
\if{html}{\figure{logo.png}{options: align='right' alt='logo' width='120'}}

The EasyMultiProfiler package aims to offer a user-friendly and efficient multi-omics data analysis tool on the R platform. It facilitates various essential tasks related to microbiome, genome, and metabolite downstream analysis, providing a seamless workflow from start to finish.
}

\section{Function}{
{Data Extract} : 
\itemize{
	\item \code{\link{EMP_assay_extract}}  
 	\item \code{\link{EMP_coldata_extract}}
 	\item \code{\link{EMP_rowdata_extract}}
}
{Data Preparation} : 
\itemize{
 	\item \code{\link{EMP_adjust_abudance}}
 	\item \code{\link{EMP_collapse}}
 	\item \code{\link{EMP_decostand}} 	
	\item \code{\link{EMP_feature_convert}}
	\item \code{\link{EMP_impute}} 
	\item \code{\link{EMP_identify_assay}}
	\item \code{\link{EMP_modify_assay}}		
	\item \code{\link{EMP_rrarefy}}  
}
{Data Analysis} : 
\itemize{
	\item \code{\link{EMP_alpha_analysis}}  
	\item \code{\link{EMP_cluster_analysis}}
	\item \code{\link{EMP_cor_analysis}}
	\item \code{\link{EMP_diff_analysis}}
	\item \code{\link{EMP_dimension_analysis}}
	\item \code{\link{EMP_enrich_analysis}}
	\item \code{\link{EMP_marker_analysis}}
	\item \code{\link{EMP_GSEA_analysis}}
	\item \code{\link{EMP_WGCNA_cluster_analysis}}
	\item \code{\link{EMP_WGCNA_cor_analysis}}
	\item \code{\link{EMP_multi_analysis}}
 }
{Data Virtualization} : 
\itemize{
	\item \code{\link{EMP_boxplot}} 
	\item \code{\link{EMP_scatterplot}}
	\item \code{\link{EMP_dotplot}}
	\item \code{\link{EMP_netplot}}
	\item \code{\link{EMP_curveplot}}  
 	\item \code{\link{EMP_heatmap_plot}}
 	\item \code{\link{EMP_sankey_plot}}
 	\item \code{\link{EMP_volcanol_plot}}

 }
{Data support} : 
\itemize{
	\item \code{\link{EMP_filter}} 
	\item \code{\link{EMP_history}} 
	\item \code{\link{EMP_result}}	
	\item \code{\link{EMP_summary}}
	\item \code{\link{EMP_taxonomy_import}} 
	\item \code{\link{EMP_function_import}} 
	\item \code{\link{EMP_normal_import}} 
	\item \code{\link{EMP_easy_import}} 
 }
}

\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/XXXXX/}
}
}


