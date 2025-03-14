#' @importFrom SummarizedExperiment colData
setGeneric(".get.mapping.EMPT",function(obj) standardGeneric(".get.mapping.EMPT"))
setMethod(".get.mapping.EMPT","EMPT",function(obj){
  colData(obj) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column('primary') %>%
  tibble::as_tibble()
})
setGeneric(".get.mapping.EMPT<-",function(obj,value) standardGeneric(".get.mapping.EMPT<-"))
setMethod(".get.mapping.EMPT<-","EMPT",function(obj,value){
  primary <- NULL
  check_primary_exist <- 'primary' %in% colnames(value)
  if (!check_primary_exist) {
    stop('New coldata must contain primary column!')
  }

  new_sample_name <- value %>% dplyr::pull(primary)
  raw_sample_name <- colData(obj)  %>% rownames()
  is_raw_contains_new <- all(new_sample_name %in% raw_sample_name) 
  if (!is_raw_contains_new) {
    stop('New primary must be included in the old coldata!')
  }

  coldata <- value %>% tibble::column_to_rownames('primary') %>% DataFrame()
  rowdata<- rowData(obj)

  assay_content <- assay(obj)
  assay_content <- assay_content %>% as.data.frame() %>% dplyr::select(dplyr::all_of(!!new_sample_name))

  coldata_reorder <- coldata[match(colnames(assay_content), rownames(coldata)),,drop = FALSE] # make sure the order match the colnames of assay!
  rowdata_reorder <- rowdata[match(rownames(assay_content), rowdata$feature),,drop = FALSE] # make sure the order match the colnames of assay!

  data.se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=as.matrix(assay_content)),
                                                        rowData=rowdata_reorder, colData = coldata_reorder)
  obj@colData <- data.se@colData
  obj@assays <-data.se@assays
  obj@NAMES <-data.se@NAMES
  obj@elementMetadata <-data.se@elementMetadata
  obj@metadata <-data.se@metadata
  obj
})


setGeneric(".get.row_info.EMPT",function(obj) standardGeneric(".get.row_info.EMPT"))
setMethod(".get.row_info.EMPT","EMPT",function(obj){
  feature <- NULL
  rowData(obj) %>% as.data.frame() %>% 
    dplyr::select(feature,everything()) %>% tibble::as_tibble()
})
setGeneric(".get.row_info.EMPT<-",function(obj,value) standardGeneric(".get.row_info.EMPT<-"))
setMethod(".get.row_info.EMPT<-","EMPT",function(obj,value){
  feature <- NULL
  check_feature_exist <- 'feature' %in% colnames(value)
  if (!check_feature_exist) {
    stop('New row_info data must contain feature column!')
  }

  new_feature_name <- value  %>% dplyr::pull(feature)
  old_feature_name <- rowData(obj) %>% rownames()
  is_raw_contains_new <- all(new_feature_name %in% old_feature_name) 
  if (!is_raw_contains_new) {
    stop('New feature must be included in the old row_info!')
  }


  assay_content <- assay(obj) %>% as.data.frame()
  assay_content <- assay_content %>% dplyr::filter(rownames(assay_content) %in% !!new_feature_name)

  rowdata <- value

  coldata <- colData(obj)

  coldata_reorder <- coldata[match(colnames(assay_content), rownames(coldata)),,drop = FALSE ] # make sure the order match the colnames of assay!
  rowdata_reorder <- rowdata[match(rownames(assay_content), rowdata$feature),,drop = FALSE] # make sure the order match the colnames of assay!

  data.se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=as.matrix(assay_content)),
                                                        rowData=rowdata_reorder, colData = coldata_reorder)
  obj@colData <- data.se@colData
  obj@assays <-data.se@assays
  obj@NAMES <-data.se@NAMES
  obj@elementMetadata <-data.se@elementMetadata
  obj@metadata <-data.se@metadata
  obj
})



setGeneric(".get.deposit_info.EMPT",function(obj) standardGeneric(".get.deposit_info.EMPT"))
setMethod(".get.deposit_info.EMPT","EMPT",function(obj){
  obj@deposit_info
})



#' Extract the existed result or inject external result into EMPT
#'
#' @param obj EMPT or EMP object.
#' @param info A character string. Result or analysis name in the EMPT or EMP object.
#' @export
#' @rdname EMP_result
#' @examples
#' 
#' data(MAE)
#' \dontrun{
#' ## obtain the result from EMPT
#' MAE |>
#'   EMP_assay_extract('geno_ec') |>
#'   EMP_alpha_analysis() |> 
#'   EMP_diff_analysis(method = 'DESeq2',.formula = ~Group) |>
#'   EMP_enrich_analysis(pvalue<0.05,keyType='ec',pvalueCutoff=0.05) -> result
#' 
#' diff_re <- result |> EMP_result(info = 'EMP_diff_analysis')
#' alpha_re <- result |> EMP_result(info = 'EMP_alpha_analysis')
#' enrich_re <- result |> EMP_result(info = 'EMP_enrich_analysis')
#' }
#' @return list
setGeneric("EMP_result",function(obj,info=NULL) standardGeneric("EMP_result"))

#' @rdname EMP_result
setMethod("EMP_result","EMPT",function(obj,info=NULL){
  if (is.null(info)) {
    info <- .get.info.EMPT(obj)
  }
  if (info %in% c('EMP_assay_data','EMP_decostand')) {
    result <- .get.assay.EMPT(obj)
    return(result)
  }else{
    result_list <- list()
    deposilt_info <- .get.deposit_info.EMPT(obj)
    if (info %in% deposilt_info$Result) {
      result_list <- obj@deposit[[info]]
    }else if(info %in% deposilt_info$source){
      real_info <- deposilt_info$Result[deposilt_info$source %in% info]
      result_list <- list()
      for (i in real_info) {
         result_list[[i]] <-obj@deposit[[i]]
      }
      if (length(result_list) == 1) {
        result_list <- result_list[[1]]
      }
    }else{
      warning("please check the info!")
    }
    return(result_list)    
  }
})

#' @rdname EMP_result
setMethod("EMP_result","EMP",function(obj,info=NULL){
  if (is.null(info)) {
    info <- .get.info.EMP(obj)
  }
  result_list <-  .get.result.EMP(obj,info=info)
  return(result_list)
})

#' @param obj EMPT or EMP object.
#' @param value A data frame or tibble from the external result.
#' @param value_name A character string. Set the name of external result to inject into EMPT project.
#' @param affect_when_sample_changed 0 or 1. 0 means that the result is not influenced by sample changes, while 1 means the contrary.
#' @param affect_when_feature_changed 0 or 1. 0 means that the result is not influenced by feature changes, while 1 means the contrary.
#' @param attribute A character string inculding primary, feature, all. This parameter indicates whether the result is about primary, feature or all.
#' @param attribute2 A character string inculding normal, diagonal or none. This parameter indicates the format of result.
#' @param source A character string. Set the name of the analysis which generate the result. (default: user_import)
#' @rdname EMP_result
#'
#' @export
#'
#' @examples
#' ## inject external result into EMPT
#' ### get a EMPT object
#' MAE |> 
#'   EMP_assay_extract('taxonomy') |>
#'   EMP_collapse(estimate_group = 'Genus',method = 'sum',
#'                collapse_by = 'row',action = 'add') -> obj  
#' 
#' ### get the raw data from the EMPT
#' MAE |> 
#'   EMP_assay_extract('taxonomy') |>
#'   EMP_collapse(estimate_group = 'Genus',method = 'sum',
#'                collapse_by = 'row',action = 'get') -> assay_data  
#' 
#' ### caculate the result from other packages
#' assay_data <- assay_data |> tibble::column_to_rownames('primary')
#' shannon_index <- vegan::diversity(assay_data,index = 'shannon') 
#' new_result <- tibble::tibble(primary=names(shannon_index),new_shannon=shannon_index)
#' 
#' ### inject the new result into EMPT object
#' EMP_result(obj,
#'            value_name = 'new_alpha',
#'            affect_when_sample_changed=0,
#'            affect_when_feature_changed=1,
#'            attribute='primary',
#'            attribute2='normal',source='user_import') <- new_result
#' 
#' obj |> 
#'   EMP_filter(sample_condition  = new_shannon >2)
setGeneric("EMP_result<-",function(obj,value_name,affect_when_sample_changed,affect_when_feature_changed,attribute,attribute2,source,value) standardGeneric("EMP_result<-"))


#' @rdname EMP_result
#'
#' @return EMPT object
setMethod("EMP_result<-","EMPT",function(obj,value_name,affect_when_sample_changed,affect_when_feature_changed,attribute,attribute2,source,value){
 deposilt_info <- .get.deposit_info.EMPT(obj)
 
 if (attribute2 == 'normal') {
   check_attribute <- attribute %in% (value %>% as.data.frame() %>% colnames())
   if (!check_attribute) {
     stop("When attribute2 == 'normal', new value must contain primary or feature or attribute must be none!")
   }
 }else if(attribute2 == 'diagonal'){
   check_diagonal <- colnames(as.data.frame(value)) == rownames(as.data.frame(value))
   if (!check_diagonal) {
     stop("New value is not a diagonal matrix or data.frame,please check!")
   }
 }else if(attribute2 == 'none'){
   if (attribute != 'all') {
     stop("When attribute2 == 'none', attribute must be all!")
   }
 }else {
  stop("Parameter attribute2 must be normal,diagonal or none!")
 }
 
 
 if (!affect_when_sample_changed %in% 0:1) {
   stop("Paramter affect_when_sample_changed must be 0 or 1!")
 }
 
 if (!affect_when_feature_changed %in% 0:1) {
   stop("Paramter affect_when_feature_changed must be 0 or 1!")
 }
 
 if (is.null(source)) {
   source <- 'user_import'
 }
 
 new_value_info <- data.frame(Result=value_name,
                              affect_when_sample_changed=affect_when_sample_changed,
                              affect_when_feature_changed=affect_when_feature_changed,
                              attribute=attribute,
                              attribute2=attribute2,
                              source=source)
 if (value_name %in% deposilt_info$Result) {
  stop('An analysis result named ',value_name,' already exists in the data!')
 }else{
  deposilt_info %<>% dplyr::bind_rows(new_value_info)
 }

 obj@deposit_info <- deposilt_info
 obj@deposit[[value_name]] <- value
 return(obj)
})






setGeneric(".get.deposit_append.EMPT",function(obj,info) standardGeneric(".get.deposit_append.EMPT"))
setMethod(".get.deposit_append.EMPT","EMPT",function(obj,info){
  obj@deposit_append[[info]]
})

setGeneric(".get.deposit_append.EMPT<-",function(obj,info,value) standardGeneric(".get.deposit_append.EMPT<-"))
setMethod(".get.deposit_append.EMPT<-","EMPT",function(obj,info,value){
  obj@deposit_append[[info]] <- value
  obj
})



setGeneric(".get.plot_deposit.EMPT",function(obj,info) standardGeneric(".get.plot_deposit.EMPT"))
setMethod(".get.plot_deposit.EMPT","EMPT",function(obj,info){
  obj@plot_deposit[[info]]
})

setGeneric(".get.plot_deposit.EMPT<-",function(obj,info,value) standardGeneric(".get.plot_deposit.EMPT<-"))
setMethod(".get.plot_deposit.EMPT<-","EMPT",function(obj,info,value){
  obj@plot_deposit[[info]] <- value
  obj
})


setGeneric(".get.rowRanges.EMPT",function(obj) standardGeneric(".get.rowRanges.EMPT"))
setMethod(".get.rowRanges.EMPT","EMPT",function(obj){
  obj@rowdata[["rowRanges"]]
})
setGeneric(".get.rowRanges.EMPT<-",function(obj,value) standardGeneric(".get.rowRanges.EMPT<-"))
setMethod(".get.rowRanges.EMPT<-","EMPT",function(obj,value){
  obj@rowdata[["rowRanges"]] <- value
  obj
})



setGeneric(".get.experiment.EMPT",function(obj) standardGeneric(".get.experiment.EMPT"))
setMethod(".get.experiment.EMPT","EMPT",function(obj){
  obj@experiment
})
setGeneric(".get.experiment.EMPT<-",function(obj,value) standardGeneric(".get.experiment.EMPT<-"))
setMethod(".get.experiment.EMPT<-","EMPT",function(obj,value){
  obj@experiment <- value
  obj
})


setGeneric(".get.assay_name.EMPT",function(obj) standardGeneric(".get.assay_name.EMPT"))
setMethod(".get.assay_name.EMPT","EMPT",function(obj){
  obj@assay_name
})
setGeneric(".get.assay_name.EMPT<-",function(obj,value) standardGeneric(".get.assay_name.EMPT<-"))
setMethod(".get.assay_name.EMPT<-","EMPT",function(obj,value){
  obj@assay_name <- value
  obj
})

setGeneric(".get.SE.EMPT",function(obj) standardGeneric(".get.SE.EMPT"))
setMethod(".get.SE.EMPT","EMPT",function(obj){
  count.da <- assay(obj) 
  sample.da<-colData(obj)
  row.da <- rowData(obj)

  coldata_reorder <- sample.da[match(colnames(count.da), rownames(sample.da)),,drop = FALSE] # make sure the order match the colnames of assay!
  rowdata_reorder <- row.da[match(rownames(count.da), row.da$feature),,drop = FALSE] # make sure the order match the colnames of assay!

  SE <- SummarizedExperiment::SummarizedExperiment(assays=list(counts = as.matrix(count.da)), colData = coldata_reorder,rowData =rowdata_reorder) 
  SE
})
setGeneric(".get.SE.EMPT<-",function(obj,value) standardGeneric(".get.SE.EMPT<-"))
setMethod(".get.SE.EMPT<-","EMPT",function(obj,value){
  obj@colData <- value@colData
  obj@assays <-value@assays
  obj@NAMES <-value@NAMES
  obj@elementMetadata <-value@elementMetadata
  obj@metadata <-value@metadata
  obj
})


setGeneric(".get.assay.EMPT",function(obj) standardGeneric(".get.assay.EMPT"))
setMethod(".get.assay.EMPT","EMPT",function(obj){
  primary <- NULL
  obj %>% 
  assay() %>% 
  as.data.frame() %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column('primary') %>% 
  dplyr::arrange(primary) %>%
  tibble::as_tibble()
})
setGeneric(".get.assay.EMPT<-",function(obj,value) standardGeneric(".get.assay.EMPT<-"))
setMethod(".get.assay.EMPT<-","EMPT",function(obj,value){
  primary <- feature <- . <- NULL
  sample_name <- value %>% dplyr::pull(primary)
  feature_name <- colnames(value)[-1]

  rowdata <- rowData(obj) %>% as.data.frame() %>% 
    dplyr::filter(feature %in% feature_name) 

  coldata <- colData(obj) %>% as.data.frame() %>% 
    dplyr::filter(rownames(.) %in% sample_name) %>% 
    dplyr::arrange(rownames(.))  ## necessary

  value <- value %>% 
    dplyr::arrange(primary) %>% ## necessary
    tibble::column_to_rownames('primary') %>% 
    t() 

  coldata_reorder <- coldata[match(colnames(value), rownames(coldata)),,drop = FALSE] # make sure the order match the colnames of assay!
  rowdata_reorder <- rowdata[match(rownames(value), rowdata$feature),,drop = FALSE] # make sure the order match the colnames of assay!
  
  data.se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=as.matrix(value)),
                                                    rowData=rowdata_reorder, colData = coldata_reorder)

  obj@colData <- data.se@colData
  obj@assays <-data.se@assays
  obj@NAMES <-data.se@NAMES
  obj@elementMetadata <-data.se@elementMetadata
  obj@metadata <-data.se@metadata
  obj
})

setGeneric(".get.estimate_group.EMPT",function(obj) standardGeneric(".get.estimate_group.EMPT"))
setMethod(".get.estimate_group.EMPT","EMPT",function(obj){
  obj@estimate_group
})
setGeneric(".get.estimate_group.EMPT<-",function(obj,value) standardGeneric(".get.estimate_group.EMPT<-"))
setMethod(".get.estimate_group.EMPT<-","EMPT",function(obj,value){
  obj@estimate_group <- value
  obj
})


setGeneric(".get.estimate_group_info.EMPT",function(obj) standardGeneric(".get.estimate_group_info.EMPT"))
setMethod(".get.estimate_group_info.EMPT","EMPT",function(obj){
  obj@estimate_group_info
})
setGeneric(".get.estimate_group_info.EMPT<-",function(obj,value) standardGeneric(".get.estimate_group_info.EMPT<-"))
setMethod(".get.estimate_group_info.EMPT<-","EMPT",function(obj,value){
  obj@estimate_group_info <- value
  obj
})



setGeneric(".get.formula.EMPT",function(obj) standardGeneric(".get.formula.EMPT"))
setMethod(".get.formula.EMPT","EMPT",function(obj){
  obj@formula
})
setGeneric(".get.formula.EMPT<-",function(obj,value) standardGeneric(".get.formula.EMPT<-"))
setMethod(".get.formula.EMPT<-","EMPT",function(obj,value){
  obj@formula <- value
  obj
})


setGeneric(".get.method.EMPT",function(obj) standardGeneric(".get.method.EMPT"))
setMethod(".get.method.EMPT","EMPT",function(obj){
  obj@method
})
setGeneric(".get.method.EMPT<-",function(obj,value) standardGeneric(".get.method.EMPT<-"))
setMethod(".get.method.EMPT<-","EMPT",function(obj,value){
  obj@method <- value
  obj
})


setGeneric(".get.message_info.EMPT",function(obj) standardGeneric(".get.message_info.EMPT"))
setMethod(".get.message_info.EMPT","EMPT",function(obj){
  for (str in obj@message_info) {
       message_wrap(str)
  }       
})
setGeneric(".get.message_info.EMPT<-",function(obj,value) standardGeneric(".get.message_info.EMPT<-"))
setMethod(".get.message_info.EMPT<-","EMPT",function(obj,value){
  obj@message_info <- value
  obj
})

setGeneric(".get.message_info.EMP",function(obj) standardGeneric(".get.message_info.EMP"))
setMethod(".get.message_info.EMP","EMP",function(obj){
  for (str in obj@message_info) {
       message_wrap(str)
  }       
})
setGeneric(".get.message_info.EMP<-",function(obj,value) standardGeneric(".get.message_info.EMP<-"))
setMethod(".get.message_info.EMP<-","EMP",function(obj,value){
  obj@message_info <- value
  obj
})


setGeneric(".get.algorithm.EMPT",function(obj) standardGeneric(".get.algorithm.EMPT"))
setMethod(".get.algorithm.EMPT","EMPT",function(obj){
  obj@algorithm
})
setGeneric(".get.algorithm.EMPT<-",function(obj,value) standardGeneric(".get.algorithm.EMPT<-"))
setMethod(".get.algorithm.EMPT<-","EMPT",function(obj,value){
  obj@algorithm <- value
  obj
})


setGeneric(".get.history.EMPT",function(obj,replace=FALSE) standardGeneric(".get.history.EMPT"))
setMethod(".get.history.EMPT","EMPT",function(obj,replace=FALSE){
  obj@history
})
setGeneric(".get.history.EMPT<-",function(obj,replace=FALSE,value) standardGeneric(".get.history.EMPT<-"))
setMethod(".get.history.EMPT<-","EMPT",function(obj,replace=FALSE,value){
  if(replace==FALSE){
  obj@history %<>% append(value)
  obj
  }else{
  obj@history <- value
  obj
  }
})


setGeneric(".get.palette.EMPT",function(obj) standardGeneric(".get.palette.EMPT"))
setMethod(".get.palette.EMPT","EMPT",function(obj){
  obj@palette
})
setGeneric(".get.palette.EMPT<-",function(obj,value) standardGeneric(".get.palette.EMPT<-"))
setMethod(".get.palette.EMPT<-","EMPT",function(obj,value){
  obj@palette <- value
  obj
})


setGeneric(".get.plot_category.EMPT",function(obj) standardGeneric(".get.plot_category.EMPT"))
setMethod(".get.plot_category.EMPT","EMPT",function(obj){
  obj@plot_category
})
setGeneric(".get.plot_category.EMPT<-",function(obj,value) standardGeneric(".get.plot_category.EMPT<-"))
setMethod(".get.plot_category.EMPT<-","EMPT",function(obj,value){
  obj@plot_category <- value
  obj
})


setGeneric(".get.plot_specific.EMPT",function(obj) standardGeneric(".get.plot_specific.EMPT"))
setMethod(".get.plot_specific.EMPT","EMPT",function(obj){
  obj@plot_specific
})
setGeneric(".get.plot_specific.EMPT<-",function(obj,value) standardGeneric(".get.plot_specific.EMPT<-"))
setMethod(".get.plot_specific.EMPT<-","EMPT",function(obj,value){
  obj@plot_specific <- value
  obj
})


setGeneric(".get.plot_info.EMPT",function(obj) standardGeneric(".get.plot_info.EMPT"))
setMethod(".get.plot_info.EMPT","EMPT",function(obj){
  obj@plot_info
})
setGeneric(".get.plot_info.EMPT<-",function(obj,value) standardGeneric(".get.plot_info.EMPT<-"))
setMethod(".get.plot_info.EMPT<-","EMPT",function(obj,value){
  obj@plot_info <- value
  obj
})


setGeneric(".get.info.EMPT",function(obj) standardGeneric(".get.info.EMPT"))
setMethod(".get.info.EMPT","EMPT",function(obj){
  obj@info
})
setGeneric(".get.info.EMPT<-",function(obj,value) standardGeneric(".get.info.EMPT<-"))
setMethod(".get.info.EMPT<-","EMPT",function(obj,value){
  obj@info <- value
  obj
})







