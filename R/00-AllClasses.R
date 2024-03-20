##' Class "EMP"
##' This class represents the ...
##' @name EMP-class
##' @docType class
##' @slot ExperimentList ExperimentList.
##' @slot deposit deposit.
##' @slot plot_deposit plot_deposit.
##' @slot message_info message_info
##' @slot history history.
##' @slot method method.
##' @slot palette palette.
##' @slot plot_category plot_category.
##' @slot plot_specific plot_specific.
##' @slot plot_info plot_info.
##' @slot info info.
##' @exportClass EMP
##' @author Bingdong Liu
##' @keywords classes
setClass("EMP",
         contains = NULL,
         slots    = c(
           'ExperimentList',
           'deposit',
           'plot_deposit',
           'message_info',
           'history',
           'method',
           'palette',
           'plot_category',
           'plot_specific',
           'plot_info',
           'info'
         ),
         prototype = list(
           ExperimentList = list(),
           deposit = list(),
           plot_deposit =list(),
           history = list(total=list())
         )
)

##' Class "EMPT"
##' This class represents the ...
##'
##'
##' @name EMPT-class
##'
##' @docType class
##' @slot deposit deposit               
##' @slot deposit2 deposit2 
##' @slot plot_deposit plot_deposit
##' @slot deposit_append deposit_append
##' @slot deposit_info deposit_info
##' @slot experiment experiment
##' @slot assay_name assay_name
##' @slot estimate_group estimate_group
##' @slot estimate_group_info estimate_group_info
##' @slot message_info message_info
##' @slot formula formula
##' @slot method method
##' @slot algorithm algorithm
##' @slot history history
##' @slot palette palette
##' @slot plot_category plot_category
##' @slot plot_specific plot_specific
##' @slot plot_info plot_info
##' @slot info info
##' @exportClass EMPT
##' @author Bingdong Liu
##' @keywords classes
setClass("EMPT",
         contains = 'SummarizedExperiment',
         slots    = c(
           'deposit',
           'deposit2',
           'plot_deposit',
           'deposit_append',
           'deposit_info',
           'experiment',
           'assay_name',
           'estimate_group',
           'estimate_group_info',
           'message_info',
           'formula',
           'method',
           'algorithm',
           'history',
           'palette',
           'plot_category',
           'plot_specific',
           'plot_info',
           'info'
         ),
         prototype = list(
           deposit = list(),
           deposit2 = list(),
           plot_deposit = list(),
           deposit_append = list(),
           deposit_info = data.frame(Result= c('diversity_result',
                                              'diff_analysis_result',
                                              'distance_result',
                                              'sample_cluster_result',
                                              'feature_cluster_result',
                                              'feature_WGCNA_cluster_result',
                                              'enrich_data',
                                              'dimension_coordinate','dimension_axis','dimension_VIP'
                                              ),
                                    affect_when_sample_changed= c(0,1,0,1,1,1,1,1,1,1),
                                    affect_when_feature_changed=c(1,0,1,1,1,1,1,1,1,1),
                                    attribute=c('primary','feature','primary','primary','feature','feature','all','primary','all','feature'),
                                    attribute2=c('normal','normal','diagonal','normal','normal','normal','none','normal','none','normal'),
                                    source=c('EMP_alpha_analysis','EMP_diff_analysis','EMP_beta_analysis',
                                             'EMP_cluster_analysis','EMP_cluster_analysis','EMP_WGCNA_cluster_analysis',
                                             'EMP_GSEA_analysis or EMP_enrich_analysis',
                                             'EMP_dimension_analysis','EMP_dimension_analysis','EMP_dimension_analysis')
           ),
           experiment = NULL,
           assay_name = NULL,
           message_info = list(),
           formula = NULL,
           method = NULL,
           algorithm = NULL,
           history = list(),
           palette = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF",
                       "#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666"))
)


setClass("EMP_assay_data",contains = c("EMPT","SummarizedExperiment"))
setClass("EMP_assay_boxplot",contains = c("EMP_assay_data","EMPT","SummarizedExperiment"))

setClass("EMP_alpha_analysis",contains = c("EMPT","SummarizedExperiment"))
setClass("EMP_alpha_analysis_boxplot",contains = c("EMP_alpha_analysis","EMPT","SummarizedExperiment"))

#setClass("EMP_beta_analysis",contains = c("EMPT","SummarizedExperiment"))
#setClass("EMP_beta_analysis_boxplot",contains = c("EMP_beta_analysis","EMPT","SummarizedExperiment"))

setClass("EMP_diff_analysis",contains = c("EMPT","SummarizedExperiment"))
setClass("EMP_diff_volcanol_plot",contains = c("EMP_diff_analysis","EMPT","SummarizedExperiment"))

setClass("EMP_dimension_analysis",contains = c("EMPT","SummarizedExperiment"))
setClass("EMP_dimension_analysis_scatterplot",contains = c("EMPT","SummarizedExperiment"))


setClass("EMP_cluster_analysis",contains = c("EMPT","SummarizedExperiment"))
setClass("EMP_WGCNA_cluster_analysis",contains = c("EMPT","SummarizedExperiment"))
setClass("EMP_WGCNA_cor_analysis",contains = c("EMPT","SummarizedExperiment"))
setClass("EMP_WGCNA_cor_heatmap",contains = c("EMP_WGCNA_cor_analysis","EMPT","SummarizedExperiment"))



setClass("EMP_enrich_analysis",contains = c("EMP_diff_analysis","EMPT","SummarizedExperiment"))
setClass("EMP_enrich_analysis_dotplot",contains = c("EMP_diff_analysis","EMPT","SummarizedExperiment","EMP_enrich_analysis"))
setClass("EMP_enrich_analysis_netplot",contains = c("EMP_diff_analysis","EMPT","SummarizedExperiment","EMP_enrich_analysis"))
setClass("EMP_enrich_analysis_curveplot",contains = c("EMP_diff_analysis","EMPT","SummarizedExperiment","EMP_enrich_analysis"))


setClass("EMP_cor_analysis",contains = c("EMP"))
setClass("EMP_WGCNA_cor_analysis2",contains = c("EMP"))
setClass("EMP_WGCNA_cor_heatmap2",contains = c("EMP_WGCNA_cor_analysis2","EMP"))










#' Title
#'
#' @param obj obj
#' @param ... ...
#' @rdname EMP_assay_extract
#'
#' @return xx object
#' @export
#'
#' @examples
#' # xx
setGeneric("EMP_assay_extract",function(obj,...) standardGeneric("EMP_assay_extract"))

#' Title
#'
#' @param MultiAssayExperiment MultiAssayExperiment
#' @rdname EMP_assay_extract
#'
#' @return xx object
#' @export
#'
#' @examples
#' # 
setMethod("EMP_assay_extract","MultiAssayExperiment",function(obj,...){
  .EMP_assay_extract_EMP(obj,...)
})

#' Title
#'
#' @param EMPT xx
#' @rdname EMP_assay_extract
#'
#' @return xx object
#' @export
#'
#' @examples
#' #
setMethod("EMP_assay_extract","EMPT",function(obj,...){
  .EMP_assay_extract_EMPT(obj,...)
})



#' Title
#'
#' @param obj obj
#' @param ... ...
#' @rdname EMP_boxplot
#'
#' @return xx object
#' @export
#'
#' @examples
#' # xx
setGeneric("EMP_boxplot",function(obj,...) standardGeneric("EMP_boxplot"))

#' Title
#'
#' @rdname EMP_boxplot
#'
#' @return xx object
#' @examples
#' #
setMethod("EMP_boxplot","EMP_alpha_analysis",function(obj,...){
  EMP_boxplot_alpha(obj,...)
})


#' Title
#'
#' @rdname EMP_boxplot
#'
#' @return xx object
#' @examples
#' #
setMethod("EMP_boxplot","EMP_assay_data",function(obj,...){
  EMP_assay_boxplot(obj,...)
})


#' Title
#'
#' @param obj obj
#' @param ... ...
#' @rdname EMP_scatterplot
#'
#' @return xx object
#' @export
#'
#' @examples
#' # 
setGeneric("EMP_scatterplot",function(obj,...) standardGeneric("EMP_scatterplot"))

#' Title
#'
#' @param EMP_dimension_analysis xx
#' @rdname EMP_scatterplot
#'
#' @return xx object
#' @export
#'
#' @examples
#' #
setMethod("EMP_scatterplot","EMP_dimension_analysis",function(obj,...){
  EMP_scatterplot_reduce_dimension(obj,...)
})


#' Title
#'
#' @param obj obj
#' @param ... ...
#' @rdname EMP_dotplot
#'
#' @return xx object
#' @export
#'
#' @examples
#' # 
setGeneric("EMP_dotplot",function(obj,...) standardGeneric("EMP_dotplot"))

#' Title
#'
#' @param EMP_enrich_analysis xx
#' @rdname EMP_dotplot
#'
#' @return xx object
#' @export
#'
#' @examples
#' #
setMethod("EMP_dotplot","EMP_enrich_analysis",function(obj,...){
  EMP_dotplot_enrich(obj,...)
})



#' Title
#'
#' @param obj obj 
#' @param ... ...
#' @rdname EMP_netplot
#'
#' @return xx object
#' @export
#'
#' @examples
#' # 
setGeneric("EMP_netplot",function(obj,...) standardGeneric("EMP_netplot"))

#' Title
#'
#' @param EMP_enrich_analysis xx
#' @rdname EMP_netplot
#'
#' @return xx object
#' @export
#'
#' @examples
#' #
setMethod("EMP_netplot","EMP_enrich_analysis",function(obj,...){
  EMP_netplot_enrich(obj,...)
})


#' Title
#'
#' @param obj obj 
#' @param ... ...
#' @rdname EMP_curveplot
#'
#' @return xx object
#' @export
#'
#' @examples
#' #
setGeneric("EMP_curveplot",function(obj,...) standardGeneric("EMP_curveplot"))

#' Title
#'
#' @rdname EMP_curveplot
#'
#' @return xx object
#' @export
#'
#' @examples
#' #
setMethod("EMP_curveplot","EMP_enrich_analysis",function(obj,...){
  EMP_curveplot_enrich(obj,...)
})

#' Title
#'
#' @param e1 An object of class EMPT.
#' @param e2 An object of class EMPT or EMP.
#' @param ... ...
#' @export
#' @method + EMPT
`+.EMPT` <- function(e1, e2, ...) {
  data_list <- list(e1, e2, ...) %>% unlist()
  result <- as.EMP(data_list)
  return(result)
}

#' Title
#'
#' @param e1 An object of class EMP.
#' @param e2 An object of class EMPT or EMP.
#' @param ... ...
#' @export
#' @method + EMP
`+.EMP` <- function(e1, e2, ...) {
  data_list <- list(e1, e2, ...) %>% unlist()
  result <- as.EMP(data_list)
  return(result)
}




