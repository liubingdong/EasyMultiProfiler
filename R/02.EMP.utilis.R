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
#' #x
setGeneric("EMP_heatmap_plot",function(obj, ...) standardGeneric("EMP_heatmap_plot"))


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
#' #x
setGeneric("EMP_WGCNA_cor_analysis",function(obj,...) standardGeneric("EMP_WGCNA_cor_analysis"))



#' @param method A character string. Methods include pearson (default), spearman and kendall.
#' @param coldata_to_assay A series of character strings. Select the column from coldata to caculate.
#' @param ... Further parameters passed to the function agricolae::correlation
#' @rdname EMP_WGCNA_cor_analysis

setMethod("EMP_WGCNA_cor_analysis","EMPT",function(obj,method='spearman',coldata_to_assay=NULL,...){
  .EMP_WGCNA_cor_analysis_EMPT_m(obj,method,coldata_to_assay,...)
})


#' @param select A character string. The experiment name in the EMP object.
#' @param method A character string. Methods include pearson (default), spearman and kendall.
#' @param ... Further parameters passed to the function agricolae::correlation
#' @rdname EMP_WGCNA_cor_analysis


setMethod("EMP_WGCNA_cor_analysis","EMP",function(obj,select=NULL,method='spearman',...){
  .EMP_WGCNA_cor_analysis_EMP_m(obj,select,method,...)
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










