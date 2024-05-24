setGeneric(".get.info.EMP",function(obj) standardGeneric(".get.info.EMP"))
setMethod(".get.info.EMP","EMP",function(obj){
  obj@info
})
setGeneric(".get.info.EMP<-",function(obj,value) standardGeneric(".get.info.EMP<-"))
setMethod(".get.info.EMP<-","EMP",function(obj,value){
  obj@info <- value
  obj
})


setGeneric(".get.ExperimentList.EMP",function(obj) standardGeneric(".get.ExperimentList.EMP"))
setMethod(".get.ExperimentList.EMP","EMP",function(obj){
  obj@ExperimentList
})

setGeneric(".get.ExperimentList.EMP<-",function(obj,value) standardGeneric(".get.ExperimentList.EMP<-"))
setMethod(".get.ExperimentList.EMP<-","EMP",function(obj,value){
  obj@ExperimentList <- value
  obj
})

setGeneric(".get.method.EMP",function(obj) standardGeneric(".get.method.EMP"))
setMethod(".get.method.EMP","EMP",function(obj){
  obj@method
})
setGeneric(".get.method.EMP<-",function(obj,value) standardGeneric(".get.method.EMP<-"))
setMethod(".get.method.EMP<-","EMP",function(obj,value){
  obj@method <- value
  obj
})

setGeneric(".get.history.EMP",function(obj,replace=FALSE,all=FALSE) standardGeneric(".get.history.EMP"))
setMethod(".get.history.EMP","EMP",function(obj,replace=FALSE,all=FALSE){
  if(all == FALSE){
    obj@history[['total']]
  }else{
    obj@history
  }
  
})
setGeneric(".get.history.EMP<-",function(obj,replace=FALSE,all=FALSE,experiment,value) standardGeneric(".get.history.EMP<-"))
setMethod(".get.history.EMP<-","EMP",function(obj,replace=FALSE,all=FALSE,experiment,value){
  if(replace==FALSE){
    if(all == FALSE){
      obj@history[['total']] %<>% append(value)
      obj      
    }else{
      obj@history[[experiment]] %<>% append(value)
      obj         
    }
  }else{
    if(all == FALSE){ 
      obj@history[['total']] <- value
      obj
    }else{
      obj@history[[experiment]] <- value
      obj
    }
  }
})


setGeneric(".get.palette.EMP",function(obj) standardGeneric(".get.palette.EMP"))
setMethod(".get.palette.EMP","EMP",function(obj){
  obj@palette
})
setGeneric(".get.palette.EMP<-",function(obj,value) standardGeneric(".get.palette.EMP<-"))
setMethod(".get.palette.EMP<-","EMP",function(obj,value){
  obj@palette <- value
  obj
})


setGeneric(".get.plot_category.EMP",function(obj) standardGeneric(".get.plot_category.EMP"))
setMethod(".get.plot_category.EMP","EMP",function(obj){
  obj@plot_category
})
setGeneric(".get.plot_category.EMP<-",function(obj,value) standardGeneric(".get.plot_category.EMP<-"))
setMethod(".get.plot_category.EMP<-","EMP",function(obj,value){
  obj@plot_category <- value
  obj
})


setGeneric(".get.plot_specific.EMP",function(obj) standardGeneric(".get.plot_specific.EMP"))
setMethod(".get.plot_specific.EMP","EMP",function(obj){
  obj@plot_specific
})
setGeneric(".get.plot_specific.EMP<-",function(obj,value) standardGeneric(".get.plot_specific.EMP<-"))
setMethod(".get.plot_specific.EMP<-","EMP",function(obj,value){
  obj@plot_specific <- value
  obj
})

setGeneric(".get.plot_info.EMP",function(obj) standardGeneric(".get.plot_info.EMP"))
setMethod(".get.plot_info.EMP","EMP",function(obj){
  obj@plot_info
})
setGeneric(".get.plot_info.EMP<-",function(obj,value) standardGeneric(".get.plot_info.EMP<-"))
setMethod(".get.plot_info.EMP<-","EMP",function(obj,value){
  obj@plot_info <- value
  obj
})


#' EMP heatmap plot
#'
#' @param obj EMPT or EMP object
#' @param ... ...
#' @rdname EMP_heatmap_plot
#'
#' @return EMPT or EMP object
#' @export
#'
#' @examples
#' \dontrun{
#' data(MAE)
#' ## for assay
#' MAE %>%
#'   EMP_assay_extract('geno_ec') %>%
#'   EMP_diff_analysis(method='DESeq2',.formula = ~Group) %>%
#'   EMP_filter(feature_condition = pvalue<0.05 & abs(fold_change) >3.5) %>%
#'   EMP_decostand(method = 'clr') %>%
#'   EMP_heatmap_plot(rotate=FALSE,palette='Spectral')
#'
#' MAE %>%
#'   EMP_assay_extract('geno_ec') %>%
#'   EMP_diff_analysis(method='DESeq2',.formula = ~Group) %>%
#'   EMP_filter(feature_condition = pvalue<0.05 & abs(fold_change) >3.5) %>%
#'   EMP_collapse(estimate_group = 'Group',collapse_by = 'col') %>% # collapse the data by group
#'   EMP_heatmap_plot(rotate=TRUE,palette='Spectral')
#'
#' ## for cor analysis
#' k1 <- MAE |>
#'   EMP_assay_extract('taxonomy') |>
#'   EMP_collapse(estimate_group = 'Genus',collapse_by = 'row') |>
#'   EMP_diff_analysis(method='DESeq2', .formula = ~Group) |>
#'   EMP_filter(feature_condition = pvalue<0.05)
#' 
#' k2 <- MAE |>
#'   EMP_collapse(experiment = 'untarget_metabol',na_string=c('NA','null','','-'),
#'                estimate_group = 'MS2kegg',method = 'sum',collapse_by = 'row') |>
#'   EMP_diff_analysis(method='DESeq2', .formula = ~Group) |>
#'   EMP_filter(feature_condition = pvalue<0.05 & abs(fold_change) > 1.5)
#' 
#' (k1 + k2) %>% EMP_cor_analysis(method = 'spearman') |>
#'   EMP_heatmap_plot() ## Visualization
#'
#' ## for WGCNA
#' MAE |>
#'   EMP_assay_extract('geno_ec')  |> 
#'   EMP_identify_assay(method = 'edgeR',estimate_group = 'Group') |>
#'   EMP_WGCNA_cluster_analysis(RsquaredCut = 0.85,mergeCutHeight=0.4)  |>
#'   EMP_WGCNA_cor_analysis(coldata_to_assay = c('BMI','PHQ9','GAD7','HAMD','SAS','SDS'),method='spearman') |>
#'   EMP_heatmap_plot(palette = 'Spectral') 
#' }
setGeneric("EMP_heatmap_plot",function(obj, ...) standardGeneric("EMP_heatmap_plot"))

#' @rdname EMP_heatmap_plot
setMethod("EMP_heatmap_plot","EMP_assay_heatmap_union",function(obj, ...){
  EMP_heatmap.EMP_assay_data(obj, ...)
})

#' @rdname EMP_heatmap_plot
setMethod("EMP_heatmap_plot","EMP_cor_analysis",function(obj, ...){
  EMP_heatmap.EMP_cor_analysis(obj, ...)
})

#' @rdname EMP_heatmap_plot
setMethod("EMP_heatmap_plot","EMP_WGCNA_cor_analysis",function(obj, ...){
  EMP_heatmap.WGCNA(obj, ...)
})


#' @rdname EMP_heatmap_plot
setMethod("EMP_heatmap_plot","EMP_WGCNA_cor_analysis2",function(obj, ...){
  EMP_heatmap.WGCNA(obj, ...)
})

#' EMP_WGCNA_cor_analysis
#'
#' @param obj EMPT or MultiAssayExperiment object.
#' @rdname EMP_WGCNA_cor_analysis
#'
#' @return EMP object
#' @export
#'
#' @examples
#' \dontrun{
#' data(MAE)
#' ## from one experiment
#' WGCNA_COR_result <- MAE |>
#'   EMP_assay_extract('geno_ec')  |> 
#'   EMP_identify_assay(method = 'edgeR',estimate_group = 'Group') |>
#'   EMP_WGCNA_cluster_analysis(RsquaredCut = 0.85,mergeCutHeight=0.4)  |>
#'   EMP_WGCNA_cor_analysis(coldata_to_assay = c('BMI','PHQ9','GAD7','HAMD','SAS','SDS'),
#'                          method='spearman',action='add') # If want the detailed result, set action = 'get'
#' 
#' ## Visualization
#' MAE |>
#'   EMP_assay_extract('geno_ec')  |> 
#'   EMP_identify_assay(method = 'edgeR',estimate_group = 'Group') |>
#'   EMP_WGCNA_cluster_analysis(RsquaredCut = 0.85,mergeCutHeight=0.4)  |>
#'   EMP_WGCNA_cor_analysis(coldata_to_assay = c('BMI','PHQ9','GAD7','HAMD','SAS','SDS'),method='spearman') |>
#'   EMP_heatmap_plot(palette = 'Spectral')
#' 
#' ## Filter the interesting module and make the enrichment analysis
#' MAE |>
#'   EMP_assay_extract('geno_ec')  |> 
#'   EMP_identify_assay(method = 'edgeR',estimate_group = 'Group') |>
#'   EMP_WGCNA_cluster_analysis(RsquaredCut = 0.85,mergeCutHeight=0.4)  |>
#'   EMP_WGCNA_cor_analysis(coldata_to_assay = c('BMI','PHQ9','GAD7','HAMD','SAS','SDS'),method='spearman') |>
#'   EMP_heatmap_plot(palette = 'Spectral') |>
#'   EMP_filter(feature_condition = WGCNA_color == 'brown' ) |> 
#'   EMP_diff_analysis(method = 'DESeq2',.formula = ~Group) |>
#'   EMP_enrich_analysis(keyType = 'ec',KEGG_Type = 'MKEGG') |>
#'   EMP_dotplot()
#' 
#' ## from two different experiments
#' k1 <- MAE |>
#'   EMP_assay_extract('geno_ec')  |> 
#'   EMP_identify_assay(method = 'edgeR',estimate_group = 'Group') |>
#'   EMP_WGCNA_cluster_analysis(RsquaredCut = 0.85,mergeCutHeight=0.4)
#' 
#' k2 <- MAE |>
#'   EMP_assay_extract('host_gene',pattern = c('A1BG','A1CF','A2MP1','AACS'),pattern_ref = 'feature')
#' 
#' (k1 + k2) |>
#'   EMP_WGCNA_cor_analysis(method='spearman') |>
#'   EMP_heatmap_plot(palette = 'Spectral') 
#' }
setGeneric("EMP_WGCNA_cor_analysis",function(obj,...) standardGeneric("EMP_WGCNA_cor_analysis"))



#' @param method A character string. Methods include pearson (default), spearman.
#' @param coldata_to_assay A series of character strings. Select the column from coldata to caculate.
#' @param action A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param ... ...
#' @rdname EMP_WGCNA_cor_analysis

setMethod("EMP_WGCNA_cor_analysis","EMPT",function(obj,method='spearman',coldata_to_assay=NULL,action='add',...){
  .EMP_WGCNA_cor_analysis_EMPT_m(obj,method,coldata_to_assay,action,...)
})


#' @param select A character string. The experiment name in the EMP object.
#' @param method A character string. Methods include pearson (default), spearman.
#' @param ... ...
#' @rdname EMP_WGCNA_cor_analysis


setMethod("EMP_WGCNA_cor_analysis","EMP",function(obj,select=NULL,method='spearman',action='add',...){
  .EMP_WGCNA_cor_analysis_EMP_m(obj,select,method,action,...)
})


setGeneric(".get.plot_deposit.EMP",function(obj,info) standardGeneric(".get.plot_deposit.EMP"))
setMethod(".get.plot_deposit.EMP","EMP",function(obj,info){
  obj@plot_deposit[[info]]
})

setGeneric(".get.plot_deposit.EMP<-",function(obj,info,value) standardGeneric(".get.plot_deposit.EMP<-"))
setMethod(".get.plot_deposit.EMP<-","EMP",function(obj,info,value){
  obj@plot_deposit[[info]] <- value
  obj
})


setGeneric(".get.experiment.EMP",function(obj,info) standardGeneric(".get.experiment.EMP"))
setMethod(".get.experiment.EMP","EMP",function(obj,info){
  experiment_name <-names(obj@ExperimentList)
  experiment_name
})


#' EMP_sankey_plot
#'
#' @param obj EMPT
#' @param ... ...
#' @return EMP object
#' @export
#'
#' @examples
#' data(MAE)
#' k1 <- MAE |>
#'   EMP_assay_extract('taxonomy') |>
#'   EMP_collapse(estimate_group = 'Genus',collapse_by = 'row') |>
#'   EMP_diff_analysis(method='DESeq2', .formula = ~Group) |>
#'   EMP_filter(feature_condition = pvalue<0.05)
#' 
#' k2 <- MAE |>
#'   EMP_collapse(experiment = 'untarget_metabol',na_string=c('NA','null','','-'),
#'                estimate_group = 'MS2kegg',method = 'sum',collapse_by = 'row') |>
#'   EMP_diff_analysis(method='DESeq2', .formula = ~Group) |>
#'   EMP_filter(feature_condition = pvalue<0.05 & abs(fold_change) > 1.5)
#'
#' k3 <- MAE |>
#'   EMP_assay_extract('geno_ec') |>
#'   EMP_diff_analysis(method='DESeq2', .formula = ~Group) |>
#'   EMP_filter(feature_condition = pvalue<0.05 & abs(fold_change) > 2)
#'
#' (k1 + k3 + k2) |> EMP_cor_analysis() |>
#'   EMP_sankey_plot()
setGeneric("EMP_sankey_plot",function(obj, ...) standardGeneric("EMP_sankey_plot"))

#' @rdname EMP_sankey_plot
setMethod("EMP_sankey_plot","EMP_cor_analysis",function(obj, ...){
  EMP_sankey_plot.EMP_cor_analysis(obj, ...)
})







