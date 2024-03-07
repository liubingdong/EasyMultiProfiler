#' @importFrom SummarizedExperiment colData
setGeneric(".get.mapping.EMPT",function(obj,...) standardGeneric(".get.mapping.EMPT"))
setMethod(".get.mapping.EMPT","EMPT",function(obj,...){
  colData(obj) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column('primary') %>%
  tibble::as_tibble()
})
setGeneric(".get.mapping.EMPT<-",function(obj,value,...) standardGeneric(".get.mapping.EMPT<-"))
setMethod(".get.mapping.EMPT<-","EMPT",function(obj,value){
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
  data.se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=as.matrix(assay_content)),
                                                        rowData=rowdata, colData = coldata)
  obj@colData <- data.se@colData
  obj@assays <-data.se@assays
  obj@NAMES <-data.se@NAMES
  obj@elementMetadata <-data.se@elementMetadata
  obj@metadata <-data.se@metadata
  obj
})


setGeneric(".get.row_info.EMPT",function(obj,...) standardGeneric(".get.row_info.EMPT"))
setMethod(".get.row_info.EMPT","EMPT",function(obj,...){
  rowData(obj) %>% as.data.frame() %>% 
    dplyr::select(feature,everything()) %>% tibble::as_tibble()
})
setGeneric(".get.row_info.EMPT<-",function(obj,value,...) standardGeneric(".get.row_info.EMPT<-"))
setMethod(".get.row_info.EMPT<-","EMPT",function(obj,value){
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
  assay_content <- assay_content %>% dplyr::filter(rownames(.) %in% !!new_feature_name)

  rowdata <- value

  coldata <- colData(obj)

  data.se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=as.matrix(assay_content)),
                                                        rowData=rowdata, colData = coldata)
  obj@colData <- data.se@colData
  obj@assays <-data.se@assays
  obj@NAMES <-data.se@NAMES
  obj@elementMetadata <-data.se@elementMetadata
  obj@metadata <-data.se@metadata
  obj
})



setGeneric(".get.deposit_info.EMPT",function(obj,...) standardGeneric(".get.deposit_info.EMPT"))
setMethod(".get.deposit_info.EMPT","EMPT",function(obj,...){
  obj@deposit_info
})

setGeneric(".get.deposit.EMPT",function(obj,info,...) standardGeneric(".get.deposit.EMPT"))
setMethod(".get.deposit.EMPT","EMPT",function(obj,info,...){
  deposilt_info <- .get.deposit_info.EMPT(obj)
  if (info %in% deposilt_info$Result) {
    obj@deposit[[info]]
  }else if(info %in% deposilt_info$source){
    real_info <- deposilt_info$Result[deposilt_info$source %in% info]
    obj@deposit[[real_info]]
  }else{
    warning("please check the info!")
  }
})

setGeneric(".get.deposit.EMPT<-",function(obj,result,affect_when_sample_changed,affect_when_feature_changed,attribute,attribute2,source,...) standardGeneric(".get.deposit.EMPT<-"))
setMethod(".get.deposit.EMPT<-","EMPT",function(obj,result,result_name,affect_when_sample_changed,affect_when_feature_changed,attribute,attribute2,source,...){
 deposilt_info <- .get.deposit_info.EMPT(obj)
 
 if (attribute2 == 'normal') {
   check_attribute <- attribute %in% (result %>% as.data.frame() %>% colnames())
   if (!check_attribute) {
     stop("When attribute2 == 'normal', new result must contain primary or feature or attribute must be none!")
   }
 }else if(attribute2 == 'diagonal'){
   check_diagonal <- colnames(as.data.frame(result)) == rownames(as.data.frame(result))
   if (!check_diagonal) {
     stop("New result is not a diagonal matrix or data.frame,please check!")
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
 
 new_result_info <- data.frame(Result=result_name,
                              affect_when_sample_changed=affect_when_sample_changed,
                              affect_when_feature_changed=affect_when_feature_changed,
                              attribute=attribute,
                              attribute2=attribute2,
                              source=source)
 deposilt_info %<>% dplyr::bind_rows(new_result_info)
 obj@deposit_info <- deposilt_info
 obj@deposit[[result_name]] <- result
 return(obj)
})


setGeneric(".get.deposit_append.EMPT",function(obj,info,...) standardGeneric(".get.deposit_append.EMPT"))
setMethod(".get.deposit_append.EMPT","EMPT",function(obj,info,...){
  obj@deposit_append[[info]]
})

setGeneric(".get.deposit_append.EMPT<-",function(obj,info,value,...) standardGeneric(".get.deposit_append.EMPT<-"))
setMethod(".get.deposit_append.EMPT<-","EMPT",function(obj,info,value,...){
  obj@deposit_append[[info]] <- value
  obj
})



setGeneric(".get.plot_deposit.EMPT",function(obj,info,...) standardGeneric(".get.plot_deposit.EMPT"))
setMethod(".get.plot_deposit.EMPT","EMPT",function(obj,info,...){
  obj@plot_deposit[[info]]
})

setGeneric(".get.plot_deposit.EMPT<-",function(obj,info,value,...) standardGeneric(".get.plot_deposit.EMPT<-"))
setMethod(".get.plot_deposit.EMPT<-","EMPT",function(obj,info,value,...){
  obj@plot_deposit[[info]] <- value
  obj
})


setGeneric(".get.rowRanges.EMPT",function(obj,...) standardGeneric(".get.rowRanges.EMPT"))
setMethod(".get.rowRanges.EMPT","EMPT",function(obj,...){
  obj@rowdata[["rowRanges"]]
})
setGeneric(".get.rowRanges.EMPT<-",function(obj,value,...) standardGeneric(".get.rowRanges.EMPT<-"))
setMethod(".get.rowRanges.EMPT<-","EMPT",function(obj,value){
  obj@rowdata[["rowRanges"]] <- value
  obj
})



setGeneric(".get.experiment.EMPT",function(obj,...) standardGeneric(".get.experiment.EMPT"))
setMethod(".get.experiment.EMPT","EMPT",function(obj,...){
  obj@experiment
})
setGeneric(".get.experiment.EMPT<-",function(obj,value,...) standardGeneric(".get.experiment.EMPT<-"))
setMethod(".get.experiment.EMPT<-","EMPT",function(obj,value){
  obj@experiment <- value
  obj
})


setGeneric(".get.assay_name.EMPT",function(obj,...) standardGeneric(".get.assay_name.EMPT"))
setMethod(".get.assay_name.EMPT","EMPT",function(obj,...){
  obj@assay_name
})
setGeneric(".get.assay_name.EMPT<-",function(obj,value,...) standardGeneric(".get.assay_name.EMPT<-"))
setMethod(".get.assay_name.EMPT<-","EMPT",function(obj,value){
  obj@assay_name <- value
  obj
})

setGeneric(".get.SE.EMPT",function(obj,...) standardGeneric(".get.SE.EMPT"))
setMethod(".get.SE.EMPT","EMPT",function(obj,...){
  count.da <- assay(obj) 
  sample.da<-colData(obj)
  row.da <- rowData(obj)
  SE <- SummarizedExperiment::SummarizedExperiment(assays=list(counts = as.matrix(count.da)), colData = sample.da,rowData =row.da) 
  SE
})
setGeneric(".get.SE.EMPT<-",function(obj,value,...) standardGeneric(".get.SE.EMPT<-"))
setMethod(".get.SE.EMPT<-","EMPT",function(obj,value){
  obj@colData <- value@colData
  obj@assays <-value@assays
  obj@NAMES <-value@NAMES
  obj@elementMetadata <-value@elementMetadata
  obj@metadata <-value@metadata
  obj
})


setGeneric(".get.assay.EMPT",function(obj,...) standardGeneric(".get.assay.EMPT"))
setMethod(".get.assay.EMPT","EMPT",function(obj,...){
  obj %>% 
  assay() %>% 
  as.data.frame() %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column('primary') %>% 
  tibble::as_tibble()
})
setGeneric(".get.assay.EMPT<-",function(obj,value,...) standardGeneric(".get.assay.EMPT<-"))
setMethod(".get.assay.EMPT<-","EMPT",function(obj,value){
  
  sample_name <- value %>% dplyr::pull(primary)
  feature_name <- colnames(value)[-1]

  rowdata <- rowData(obj) %>% as.data.frame() %>% 
    dplyr::filter(feature %in% feature_name) 

  coldata <- colData(obj)
  coldata <- coldata[rownames(coldata) %in% sample_name, ] %>% as.data.frame()

  value %<>% tibble::column_to_rownames('primary') %>% t() %>% DataFrame()

  data.se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=as.matrix(value)),
                                                    rowData=rowdata, colData = coldata)

  obj@colData <- data.se@colData
  obj@assays <-data.se@assays
  obj@NAMES <-data.se@NAMES
  obj@elementMetadata <-data.se@elementMetadata
  obj@metadata <-data.se@metadata
  obj
})

setGeneric(".get.estimate_group.EMPT",function(obj,...) standardGeneric(".get.estimate_group.EMPT"))
setMethod(".get.estimate_group.EMPT","EMPT",function(obj,...){
  obj@estimate_group
})
setGeneric(".get.estimate_group.EMPT<-",function(obj,value,...) standardGeneric(".get.estimate_group.EMPT<-"))
setMethod(".get.estimate_group.EMPT<-","EMPT",function(obj,value){
  obj@estimate_group <- value
  obj
})


setGeneric(".get.estimate_group_info.EMPT",function(obj,...) standardGeneric(".get.estimate_group_info.EMPT"))
setMethod(".get.estimate_group_info.EMPT","EMPT",function(obj,...){
  obj@estimate_group_info
})
setGeneric(".get.estimate_group_info.EMPT<-",function(obj,value,...) standardGeneric(".get.estimate_group_info.EMPT<-"))
setMethod(".get.estimate_group_info.EMPT<-","EMPT",function(obj,value){
  obj@estimate_group_info <- value
  obj
})



setGeneric(".get.formula.EMPT",function(obj,...) standardGeneric(".get.formula.EMPT"))
setMethod(".get.formula.EMPT","EMPT",function(obj,...){
  obj@formula
})
setGeneric(".get.formula.EMPT<-",function(obj,value,...) standardGeneric(".get.formula.EMPT<-"))
setMethod(".get.formula.EMPT<-","EMPT",function(obj,value){
  obj@formula <- value
  obj
})


setGeneric(".get.method.EMPT",function(obj,...) standardGeneric(".get.method.EMPT"))
setMethod(".get.method.EMPT","EMPT",function(obj,...){
  obj@method
})
setGeneric(".get.method.EMPT<-",function(obj,value,...) standardGeneric(".get.method.EMPT<-"))
setMethod(".get.method.EMPT<-","EMPT",function(obj,value){
  obj@method <- value
  obj
})


setGeneric(".get.message_info.EMPT",function(obj,...) standardGeneric(".get.message_info.EMPT"))
setMethod(".get.message_info.EMPT","EMPT",function(obj,...){
  for (str in obj@message_info) {
       message_wrap(str)
  }       
})
setGeneric(".get.message_info.EMPT<-",function(obj,value,...) standardGeneric(".get.message_info.EMPT<-"))
setMethod(".get.message_info.EMPT<-","EMPT",function(obj,value){
  obj@message_info <- value
  obj
})

setGeneric(".get.message_info.EMP",function(obj,...) standardGeneric(".get.message_info.EMP"))
setMethod(".get.message_info.EMP","EMP",function(obj,...){
  for (str in obj@message_info) {
       message_wrap(str)
  }       
})
setGeneric(".get.message_info.EMP<-",function(obj,value,...) standardGeneric(".get.message_info.EMP<-"))
setMethod(".get.message_info.EMP<-","EMP",function(obj,value){
  obj@message_info <- value
  obj
})


setGeneric(".get.algorithm.EMPT",function(obj,...) standardGeneric(".get.algorithm.EMPT"))
setMethod(".get.algorithm.EMPT","EMPT",function(obj,...){
  obj@algorithm
})
setGeneric(".get.algorithm.EMPT<-",function(obj,value,...) standardGeneric(".get.algorithm.EMPT<-"))
setMethod(".get.algorithm.EMPT<-","EMPT",function(obj,value){
  obj@algorithm <- value
  obj
})


setGeneric(".get.history.EMPT",function(obj,replace=FALSE,...) standardGeneric(".get.history.EMPT"))
setMethod(".get.history.EMPT","EMPT",function(obj,replace=FALSE,...){
  obj@history
})
setGeneric(".get.history.EMPT<-",function(obj,value,replace=FALSE,...) standardGeneric(".get.history.EMPT<-"))
setMethod(".get.history.EMPT<-","EMPT",function(obj,value,replace=FALSE){
  if(replace==FALSE){
  obj@history %<>% append(value)
  obj
  }else{
  obj@history <- value
  obj
  }
})


setGeneric(".get.palette.EMPT",function(obj,...) standardGeneric(".get.palette.EMPT"))
setMethod(".get.palette.EMPT","EMPT",function(obj,...){
  obj@palette
})
setGeneric(".get.palette.EMPT<-",function(obj,value,...) standardGeneric(".get.palette.EMPT<-"))
setMethod(".get.palette.EMPT<-","EMPT",function(obj,value){
  obj@palette <- value
  obj
})


setGeneric(".get.plot_category.EMPT",function(obj,...) standardGeneric(".get.plot_category.EMPT"))
setMethod(".get.plot_category.EMPT","EMPT",function(obj,...){
  obj@plot_category
})
setGeneric(".get.plot_category.EMPT<-",function(obj,value,...) standardGeneric(".get.plot_category.EMPT<-"))
setMethod(".get.plot_category.EMPT<-","EMPT",function(obj,value){
  obj@plot_category <- value
  obj
})


setGeneric(".get.plot_specific.EMPT",function(obj,...) standardGeneric(".get.plot_specific.EMPT"))
setMethod(".get.plot_specific.EMPT","EMPT",function(obj,...){
  obj@plot_specific
})
setGeneric(".get.plot_specific.EMPT<-",function(obj,value,...) standardGeneric(".get.plot_specific.EMPT<-"))
setMethod(".get.plot_specific.EMPT<-","EMPT",function(obj,value){
  obj@plot_specific <- value
  obj
})


setGeneric(".get.plot_info.EMPT",function(obj,...) standardGeneric(".get.plot_info.EMPT"))
setMethod(".get.plot_info.EMPT","EMPT",function(obj,...){
  obj@plot_info
})
setGeneric(".get.plot_info.EMPT<-",function(obj,value,...) standardGeneric(".get.plot_info.EMPT<-"))
setMethod(".get.plot_info.EMPT<-","EMPT",function(obj,value){
  obj@plot_info <- value
  obj
})


setGeneric(".get.info.EMPT",function(obj,...) standardGeneric(".get.info.EMPT"))
setMethod(".get.info.EMPT","EMPT",function(obj,...){
  obj@info
})
setGeneric(".get.info.EMPT<-",function(obj,value,...) standardGeneric(".get.info.EMPT<-"))
setMethod(".get.info.EMPT<-","EMPT",function(obj,value){
  obj@info <- value
  obj
})







