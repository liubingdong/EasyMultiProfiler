.perform_operation <- function(x, method,na.rm=TRUE,...) {
    result <- switch(method,
                     "mean" = mean(x,na.rm=na.rm),
                     "sum" = sum(x,na.rm=na.rm),
                     "median" = median(x,na.rm=na.rm),
                     "max" = max(x,na.rm=na.rm),
                     "min" = min(x,na.rm=na.rm))
    return(result)
}

.pattern_Dectect_multi <- function(pattern_ref,pattern,exact=F){
  if (length(pattern) ==1) {
    .pattern_Dectect(pattern_ref = pattern_ref,pattern = pattern,exact=exact)
  }else if (length(pattern) > 1){
    id_detect <-list()
    for (j in pattern) {
      id_temp <- .pattern_Dectect(pattern_ref = pattern_ref,pattern = j,exact=exact)
      id_detect[[j]] <- id_temp
    }
    .combine_logical_vectors(id_detect)
  }
}


.combine_logical_vectors <- function(vector_list) {
  combined_vec <- Reduce(`|`, vector_list)

  return(combined_vec)
}


.pattern_Dectect <- function(pattern_ref,pattern,exact=F){
  if (exact == F) {
    pattern_ref <- tolower(pattern_ref)
    pattern <- tolower(pattern)
    id <- grepl(x=pattern_ref, pattern = pattern,fixed = TRUE)
  }else {
    pattern_ref <- tolower(pattern_ref)
    id <- c(pattern_ref == tolower(pattern))
  }
  return(id)
}




.as.EMPT <- function(x,experiment,estimate_group = NULL) {

  sampleMap <- x@sampleMap %>% tibble::as_tibble() %>%
    dplyr::filter(assay %in% experiment) %>%
    dplyr::select(colname,primary)

  assay_name <- names(assays(x[[experiment]]))[1]

  assay_data <- assays(x[[experiment]])[[1]] %>% t() %>% as.data.frame() %>%
    tibble::rownames_to_column(var = 'colname') %>%
    dplyr::left_join(.,sampleMap,by = 'colname') %>%
    dplyr::arrange(primary) %>%  ## confrim the sample order
    tibble::column_to_rownames('primary') %>%
    dplyr::select(-colname) %>% t()


  mapping <- x %>%
      EMP_coldata_extract(experiment = experiment) %>%
      dplyr::arrange(primary) %>%  ## confrim the sample order
      tibble::column_to_rownames('primary') %>% DataFrame()

  rowdata <- rowData(x[[experiment]])

  data.se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=as.matrix(assay_data)),
                      rowData=rowdata, colData = mapping)


  EMPT <- new(Class = 'EMPT',
              data.se,
              assay_name = assay_name,
              experiment = experiment,
              info = 'EMP_assay_data',
              estimate_group = estimate_group)
  return(EMPT)
}

as.EMP <- function(object,select=NULL) {
  call <- match.call()
  deposit <- new(Class = 'EMP')

  if (inherits(object,'MultiAssayExperiment')) {
    if (is.null(select)) {
      stop("When input data is MultiAssayExperiment, parameter select need to be specified! ")
    }else{
      data_list <- list()
      for (i in select) {
        data_list[[i]] <- object %>% .as.EMPT(i)
      }
    }
  }else if(inherits(object,'list')) {
    data_list <- object
  }else {
    stop("Please check the input data!")
  }

  experiment_name <- c()
  for (i in data_list) {
    experiment_name <- append(experiment_name,.get.experiment.EMPT(i))
  }

  experiment_name <- .check_duplicate(experiment_name)

  for (i in 1:length(data_list)) {
    .get.ExperimentList.EMP(deposit)[[experiment_name[i]]] <- data_list[[i]]
    .get.history.EMP(deposit,all = TRUE,experiment=experiment_name[[i]]) <- .get.history.EMPT(data_list[[i]])
   }
  .get.history.EMP(deposit,all = FALSE) <- call
  .get.info.EMP(deposit) <- 'EMP_list_data'
  return(deposit)
}

.group_level_modified <- function(EMPT,estimate_group,group_level = 'default') {

  name_group <- unique(.get.mapping.EMPT(EMPT)[[estimate_group]])
  group_level_check=all(name_group%in%group_level)&all(group_level%in%name_group)
  if (group_level_check==T) {
    .get.mapping.EMPT(EMPT)[[estimate_group]] <- factor(.get.mapping.EMPT(EMPT)[[estimate_group]],levels = group_level)
  }else{
    if (group_level[1] != 'default') {
      warning('group level can not match, and pipe will follow the default level !')
    }
  }
  EMPT
}


## This '.return_wrap' and 'message_wrap' code are from package 'MicrobiotaProcess' for message print
#' @importFrom glue glue_collapse
 .return_wrap <- function(...){
   msg <- paste(..., collapse = "", sep = "")
   wrapped <- strwrap(msg, width = getOption("width") - 2) %>%
     glue::glue_collapse(., "\n", last = "\n")
   wrapped
 }
 message_wrap <- function(...){
   msg <- .return_wrap(...)
   message(msg)
 }

.check_estimate_group.EMPT <- function(EMPT,estimate_group) {
  if (is.null(estimate_group)) {
    if (!is.null(.get.estimate_group.EMPT(EMPT))) {
      estimate_group <- .get.estimate_group.EMPT(EMPT)
    }else {
      stop('Please input the estimate_group paramater!')
    }
  }else {
    estimate_group = estimate_group
  }
  return(estimate_group)
}

.check_group_level.EMPT <- function(EMPT,group_level,estimate_group) {
  if (is.null(group_level)) {
    if (!is.null(.get.estimate_group_info.EMPT(EMPT))) {
      group_level <- .get.estimate_group_info.EMPT(EMPT) %>%
        strsplit(' vs ') %>%
        unlist() %>% trimws()
    }else {
      estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)
      group_level <- .get.mapping.EMPT(EMPT) %>% dplyr::pull(!!estimate_group) %>% unique() %>% rev()
    }
  }else {
    group_level = group_level
  }
  return(group_level)
}

.get.kegg_data <- function(type,use_cached=T) {
  switch(type,
         "KO" = {
           if (!use_cached){
             memoise::forget(gson_KO_pathway) %>% invisible()
           }
           gason_data <- gson_KO_pathway()
         },
         "KO_module" = {
           if (!use_cached){
             memoise::forget(gson_KO_module) %>% invisible()
           }
           gason_data <- gson_KO_module()
         },
         "EC" = {
           if (!use_cached){
             memoise::forget(gson_EC_pathway) %>% invisible()
           }
           gason_data <- gson_EC_pathway()
         },
         "EC_module" = {
           if (!use_cached){
             memoise::forget(gson_EC_module) %>% invisible()
           }
           gason_data <- gson_EC_module()
         },
         "compound" = {
           if (!use_cached){
             memoise::forget(gson_cpd_pathway) %>% invisible()
           }
           gason_data <- gson_cpd_pathway()
         },
         {
           stop('Parameter type must be one of KO,KO_module,EC,EC_module,compound!')
         }
  )
  message('KEGG database version: ',gason_data@version)
  return(gason_data)
}

.check_duplicate <- function(string_vector){
  dup_indices <- duplicated(string_vector)
  if (any(dup_indices)) {
    # 在重复元素后添加序号
    string_vector[dup_indices] <- paste0(string_vector[dup_indices],
                                           ave(seq_along(string_vector),
                                               string_vector,
                                               FUN = seq_along)[dup_indices])

  }
  return(string_vector)
}

