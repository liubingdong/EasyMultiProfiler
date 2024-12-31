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
           palette = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF",
                       "#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666"),
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
                                                'sample_cluster_result',
                                                'feature_cluster_result',
                                                'feature_WGCNA_cluster_result',
                                                'enrich_data','enrich_data',
                                                'dimension_coordinate','dimension_axis','dimension_VIP',
                                                'Boruta_model','Boruta_feature_importance',
                                                'rf_model','rf_feature_importance',
                                                'xgb_model','xgb_feature_importance',
                                                'lasso_model','lasso_feature_importance'),
                                      affect_when_sample_changed= c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                                      affect_when_feature_changed=c(1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                                      attribute=c('primary',
                                                  'feature',
                                                  'primary',
                                                  'feature',
                                                  'feature',
                                                  'all','all',
                                                  'primary','all','feature',
                                                  'all','feature',
                                                  'all','feature',
                                                  'all','feature',
                                                  'all','feature'),
                                      attribute2=c('normal','normal','normal','normal','normal','none','none',
                                                   'normal','none','normal',
                                                   'none','normal',
                                                   'none','normal',
                                                   'none','normal',
                                                   'none','normal'),
                                      source=c('EMP_alpha_analysis','EMP_diff_analysis',
                                               'EMP_cluster_analysis','EMP_cluster_analysis','EMP_WGCNA_cluster_analysis',
                                               'EMP_GSEA_analysis','EMP_enrich_analysis',
                                               'EMP_dimension_analysis','EMP_dimension_analysis','EMP_dimension_analysis',
                                               'EMP_marker_analysis','EMP_marker_analysis',
                                               'EMP_marker_analysis','EMP_marker_analysis',
                                               'EMP_marker_analysis','EMP_marker_analysis',
                                               'EMP_marker_analysis','EMP_marker_analysis')
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
setClass("EMP_structure_plot",contains = c("EMP_assay_data","EMPT","SummarizedExperiment"))
setClass("EMP_fitline_plot",contains = c("EMP_assay_data","EMPT","SummarizedExperiment"))

setClass("EMP_alpha_analysis",contains = c("EMPT","SummarizedExperiment"))
setClass("EMP_alpha_analysis_boxplot",contains = c("EMP_alpha_analysis","EMPT","SummarizedExperiment"))

setClass("EMP_diff_analysis",contains = c("EMPT","SummarizedExperiment"))
setClass("EMP_diff_volcanol_plot",contains = c("EMP_diff_analysis","EMPT","SummarizedExperiment"))

setClass("EMP_dimension_analysis",contains = c("EMPT","SummarizedExperiment"))
setClass("EMP_dimension_analysis_scatterplot",contains = c("EMPT","SummarizedExperiment"))

setClass("EMP_decostand",contains = c("EMPT","SummarizedExperiment"))
setClass("EMP_marker_analysis",contains = c("EMPT","SummarizedExperiment"))
setClass("EMP_cluster_analysis",contains = c("EMPT","SummarizedExperiment"))

setClass("EMP_WGCNA_cluster_analysis",contains = c("EMPT","SummarizedExperiment"))
setClass("EMP_WGCNA_cor_analysis",contains = c("EMPT","SummarizedExperiment"))
setClass("EMP_WGCNA_cor_heatmap",contains = c("EMP_WGCNA_cor_analysis","EMPT","SummarizedExperiment"))

setClassUnion("EMP_assay_heatmap_union", c("EMP_assay_data","EMP_decostand","EMP_diff_analysis"))
setClassUnion("EMP_assay_boxplot_union", c("EMP_assay_data","EMP_decostand","EMP_diff_analysis"))


setClass("EMP_enrich_analysis",contains = c("EMP_diff_analysis","EMPT","SummarizedExperiment"))
setClass("EMP_enrich_analysis_dotplot",contains = c("EMP_diff_analysis","EMPT","SummarizedExperiment","EMP_enrich_analysis"))
setClass("EMP_enrich_analysis_netplot",contains = c("EMP_diff_analysis","EMPT","SummarizedExperiment","EMP_enrich_analysis"))
setClass("EMP_enrich_analysis_curveplot",contains = c("EMP_diff_analysis","EMPT","SummarizedExperiment","EMP_enrich_analysis"))


setClass("EMP_cor_analysis",contains = c("EMP"))
setClass("EMP_WGCNA_cor_analysis2",contains = c("EMP"))
setClass("EMP_WGCNA_cor_heatmap2",contains = c("EMP_WGCNA_cor_analysis2","EMP"))

setClass("EMP_assay_heatmap",contains = c("EMPT","SummarizedExperiment"))
setClass("EMP_cor_heatmap",contains = c("EMP"))
setClass("EMP_cor_sankey",contains = c("EMP"))
setClass("EMP_fitline_plot2",contains = c("EMP"))

setClass("EMP_multi_same_df",contains = c("EMP"))

setClass("EMP_multi_diff_enrich",contains = c("EMP"))
setClass("EMP_multi_diff_enrich_dotplot",contains = c("EMP","EMP_multi_diff_enrich"))
setClass("EMP_multi_diff_enrich_netplot",contains = c("EMP","EMP_multi_diff_enrich"))

setClass("EMP_multi_same_enrich",contains = c("EMP"))
setClass("EMP_multi_same_enrich_dotplot",contains = c("EMP","EMP_multi_same_enrich"))
setClass("EMP_multi_same_enrich_netplot",contains = c("EMP","EMP_multi_same_enrich"))


#' Boxplot for EMPT result
#'
#' @param ... Further parameters passed to the function ggsignif::geom_signif
#' @rdname EMP_boxplot
#'
#' @return EMPT object
#' @export
#'
#' @examples
#' data(MAE)
#' ## from assay
#' MAE |> 
#'   EMP_assay_extract('host_gene',pattern = 'A1BG',pattern_ref = 'feature') |>
#'   EMP_boxplot(method='t.test',estimate_group='Group')
#' 
#' ## from alpha analysis
#' MAE |> 
#'   EMP_assay_extract('taxonomy') |> 
#'   EMP_alpha_analysis()|>
#'   EMP_boxplot(method='t.test',estimate_group='Group')
setGeneric("EMP_boxplot",function(obj, ...) standardGeneric("EMP_boxplot"))


#' @param ... ...
#' @rdname EMP_boxplot

setMethod("EMP_boxplot","EMP_alpha_analysis",function(obj, ...){
  EMP_boxplot.EMP_alpha_analysis(obj, ...)
})

#' @param ... ...
#' @rdname EMP_boxplot

setMethod("EMP_boxplot","EMP_assay_boxplot_union",function(obj, ...){
  EMP_boxplot.EMP_assay_boxplot_union(obj, ...)
})


#' Scatterplot for EMPT result
#'
#' @param obj object
#' @param ... Further parameters passed to the function ggsignif::geom_signif
#' @rdname EMP_scatterplot
#'
#' @export
#'
#' @examples
#' data(MAE)
#' MAE |> 
#'   EMP_assay_extract('taxonomy') |> 
#'   EMP_collapse(estimate_group = 'Species',collapse_by = 'row')|>
#'   EMP_dimension_analysis(method = 'pcoa',distance = 'bray')|>
#'   EMP_scatterplot(show='p12html',estimate_group = 'Group') # eg. p12,p12html,p23,p23html
setGeneric("EMP_scatterplot",function(obj,...) standardGeneric("EMP_scatterplot"))

#' @param obj object
#' @param ... ...
#' @rdname EMP_scatterplot
#'
#' @export
#'

setMethod("EMP_scatterplot","EMP_dimension_analysis",function(obj,...){
  EMP_scatterplot.EMP_dimension_analysis(obj,...)
})


#' Dotplot for enrichment result
#'
#' @param obj object
#' @param ... ...
#' @rdname EMP_dotplot
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(MAE)
#' MAE |>
#'   EMP_GSEA_analysis(experiment = 'geno_ko',method='signal2Noise',
#'                     estimate_group = 'Group',
#'                     pvalueCutoff = 0.05,keyType = 'ko') |>
#'   EMP_dotplot(color='p.adjust',showCategory=10) 
#' }
setGeneric("EMP_dotplot",function(obj,...) standardGeneric("EMP_dotplot"))


#' @param EMP_enrich_analysis EMP_enrich_analysis
#' @rdname EMP_dotplot
#'
#' @export
#' @return Enrichment dotplot

setMethod("EMP_dotplot","EMP_enrich_analysis",function(obj,...){
  EMP_dotplot_enrich(obj,...)
})

#' @param EMP_enrich_analysis EMP_multi_enrich
#' @rdname EMP_dotplot
#'
#' @export
#'
setMethod("EMP_dotplot","EMP_multi_diff_enrich",function(obj,...){
  EMP_dotplot_enrich(obj,...)
})

#' @param EMP_enrich_analysis EMP_multi_enrich
#' @rdname EMP_dotplot
#'
#' @export
#'
setMethod("EMP_dotplot","EMP_multi_same_enrich",function(obj,...){
  EMP_dotplot_enrich(obj,...)
})




#' Netplot for enrichment result
#' @param obj object
#' @param ... ...
#' @rdname EMP_netplot
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(MAE)
#' MAE |>
#'   EMP_GSEA_analysis(experiment = 'geno_ko',method='signal2Noise',
#'                     estimate_group = 'Group',
#'                     pvalueCutoff = 0.05,keyType = 'ko') |>
#'   EMP_netplot(showCategory=10) 
#' }
setGeneric("EMP_netplot",function(obj,...) standardGeneric("EMP_netplot"))



#' @rdname EMP_netplot
#' @return Enrichment netplot object
#' @export
#'

setMethod("EMP_netplot","EMP_enrich_analysis",function(obj,...){
  EMP_netplot_enrich(obj,...)
})


#' @rdname EMP_netplot
#' @export
#'

setMethod("EMP_netplot","EMP_multi_same_enrich",function(obj,...){
  EMP_netplot_enrich(obj,...)
})


#' @rdname EMP_netplot
#' @export
#'

setMethod("EMP_netplot","EMP_multi_diff_enrich",function(obj,...){
  EMP_netplot_enrich(obj,...)
})




#' Curveplot for enrichment result
#'
#' @param obj object
#' @param ... ...
#' @rdname EMP_curveplot
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(MAE)
#' MAE |>
#'  EMP_GSEA_analysis(experiment = 'geno_ko',method='signal2Noise',
#'                    estimate_group = 'Group',
#'                    pvalueCutoff = 0.05,keyType = 'ko') |>
#'  EMP_curveplot(geneSetID='map00680')
#' }
setGeneric("EMP_curveplot",function(obj,...) standardGeneric("EMP_curveplot"))


#' @rdname EMP_curveplot
#'
#' @export
#'

setMethod("EMP_curveplot","EMP_enrich_analysis",function(obj,...){
  EMP_curveplot_enrich(obj,...)
})

#' Two variable regression plot
#'
#' @param obj object
#' @param ... ...
#' @rdname EMP_fitline_plot
#'
#' @export
#'
#' @examples
#' data(MAE)
#' ## For EMPT
#' MAE |>
#'   EMP_assay_extract('geno_ec') |>
#'   EMP_fitline_plot(var_select=c('1.1.1.1','BMI'))
#'
#' MAE |>
#'   EMP_assay_extract('geno_ec') |>
#'   EMP_fitline_plot(var_select=c('Weight','BMI'),estimate_group='Sex',show='html')
#' 
#' ## For EMP
#' k1 <- MAE |>
#'   EMP_assay_extract('geno_ec') |>
#'   EMP_decostand(method = 'clr')
#' 
#' k2 <- MAE |>
#'   EMP_assay_extract('taxonomy') |>
#'   EMP_collapse(collapse_by='row',estimate_group = 'Class',method = 'sum') |>
#'   EMP_decostand(method = 'relative')
#' 
#' (k1 + k2) |>
#'   EMP_fitline_plot(var_select=c('1.1.1.1','Bacilli'),
#'                    estimate_group='Group',eq_size=2.5)

setGeneric("EMP_fitline_plot",function(obj, ...) standardGeneric("EMP_fitline_plot"))

#' @rdname EMP_fitline_plot
setMethod("EMP_fitline_plot","EMPT",function(obj, ...){
  EMP_fitline_plot.EMPT(obj, ...)
})

#' @rdname EMP_fitline_plot
setMethod("EMP_fitline_plot","EMP",function(obj, ...){
  EMP_fitline_plot.EMP(obj, ...)
})



#' @export
`+.EMPT` <- function(e1, e2, ...) {
  data_list <- list(e1, e2, ...) %>% unlist()
  result <- as.EMP(data_list)
  return(result)
}


#' @export
`+.EMP` <- function(e1, e2, ...) {
  data_list <- list(e1, e2, ...) %>% unlist()
  result <- as.EMP(data_list)
  return(result)
}




