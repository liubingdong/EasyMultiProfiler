.perform_operation <- function(x, method,na.rm=TRUE,...) {
    result <- switch(method,
                     "mean" = mean(x,na.rm=na.rm, ...),
                     "sum" = sum(x,na.rm=na.rm, ...),
                     "median" = median(x,na.rm=na.rm, ...),
                     "max" = max(x,na.rm=na.rm, ...),
                     "min" = min(x,na.rm=na.rm, ...))
    return(result)
}

#' Title
#'
#' @param pattern_ref Object in EMPT format.
#' @param pattern A character string. Experiment name in the MultiAssayExperiment object.
#' @param exact A boolean. Whether the extract use exact search method.
#' @return xx object
#' @export
#'
#' @examples
#' # add example
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



#' @importFrom SummarizedExperiment assays
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment rowData
.as.EMPT <- function(x,experiment,estimate_group = NULL) {
  assay <- colname <- primary <- NULL
  sampleMap <- x@sampleMap %>% tibble::as_tibble() %>%
    dplyr::filter(assay %in% experiment) %>%
    dplyr::select(colname,primary)

  assay_name <- names(assays(x[[experiment]]))[1]

  assay_data <- assays(x[[experiment]])[[1]] %>% t() %>% as.data.frame() %>%
    tibble::rownames_to_column(var = 'colname') %>%
    dplyr::left_join(sampleMap,by = 'colname') %>%
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
    if (inherits(i,'EMP')) {
      each_experiment_name <- .get.experiment.EMP(i)
    }else if (inherits(i,'EMPT')) {
      each_experiment_name <- .get.experiment.EMPT(i)
    }
    experiment_name  <- append(experiment_name,each_experiment_name)
  }
  experiment_name <- .check_duplicate(experiment_name)
  
  real_data_list <- list()
  for (i in data_list) {
    if (inherits(i,'EMP')) {
      each_data <- i@ExperimentList
    }else if (inherits(i,'EMPT')) {
      each_data <- i
    }
    counts <- length(real_data_list) +1
    real_data_list[[counts]] <- each_data
  }
  real_data_list %<>% unlist()
  
  for (i in 1:length(real_data_list)) {
    .get.ExperimentList.EMP(deposit)[[experiment_name[i]]] <- real_data_list[[i]]
    .get.history.EMP(deposit,all = TRUE,experiment=experiment_name[[i]]) <- .get.history.EMPT(real_data_list[[i]])
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
     glue::glue_collapse("\n", last = "\n")
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


.check_duplicate <- function(string_vector){
  dup_indices <- duplicated(string_vector)
  if (any(dup_indices)) {
    # 在重复元素后添加序号
    for(i in which(dup_indices)){
      string_original <- string_vector[i]
      string_counter <- 1
      while(string_original %in% string_vector){
        string_counter <- string_counter + 1
        string_original <- paste0(string_vector[i], string_counter)
      }
      string_vector[i] <- string_original
    }
  }
  return(string_vector)
}


#' Get history from EMPT or EMP object
#'
#' @param obj Object in EMPT or EMP format.
#'
#' @return EMPT object
#' @export
#'
#' @examples
#' data(MAE)
#' ## from EMPT
#' MAE |>
#'   EMP_assay_extract(experiment = 'taxonomy') |>
#'   EMP_alpha_analysis() |>
#'   EMP_history()
#' 
#' ## from EMP
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
#' (k1 + k2) |> EMP_cor_analysis(method = 'spearman') |>
#'   EMP_heatmap_plot()  |>
#'   EMP_history()  ## get all the history about what happened to the object

EMP_history <- function(obj) {
  if (inherits(obj,"EMP")) {
    his_info2 <- .get.history.EMP(obj,all = T) 
    his_names <- names(his_info2)
    for (i in 1:length(his_info2)) {
      his_info_each <- his_info2[[i]] 
      for (j in 1:length(his_info_each)) {
        his_info_each[j] <- gsub('\"', '', his_info_each[j])
      }
      his_info2[[his_names[i]]] <- his_info_each %>% unlist()
    }
    return(his_info2)
  }else if(inherits(obj,'EMPT')) {
    his_info <- .get.history.EMPT(obj) %>% as.character()
    for (i in 1:length(his_info)) {
      his_info[i] <-  gsub('\"', '', his_info[i])
    }
    return(his_info)
  }
}


.merge_tax_value <- function(row) {
  last_non_na <- NA
  unclassified_count <- 0
  for (i in 1:length(row)) {
    if (!is.na(row[i])) {
      last_non_na <- row[i]
      unclassified_count <- 0
    } else {
      if (!is.na(last_non_na)) {
        unclassified_count <- unclassified_count + 1
        unclassified <- paste(rep("unclassified", unclassified_count), collapse = "_")
        row[i] <- paste(last_non_na, unclassified, sep = "_")
      }
    }
  }
  return(row)
}

.impute_tax <- function(df) {
  # 逐行应用函数
  df_merged <- t(apply(df, 1, .merge_tax_value))
  
  # 更新数据框
  df <- as.data.frame(df_merged)
  
  return(df)
}


## Enhance the print for EMP_assay_data
modify_tbl_format_setup <- function(x, totalX, ...){ 
  tmpxx <- x$tbl_sum %>%
    strsplit(" ") %>%
    unlist()
  tmpxx <- paste(getFromNamespace("big_mark", "pillar")(totalX), tmpxx[2], tmpxx[3])
  names(tmpxx) <- c("A EMPT-tibble (EMPT object) size")
  x$tbl_sum <- tmpxx
  x$rows_total <- totalX
  x$rows_missing <- totalX - nrow(x$df)
  return(x)
}

modify_tbl_format_footer <- function(x,EMPT,...){
    result_num <-  length(EMPT@deposit)
    advice <- paste0("The obeject contains ", result_num," analysis result.")
    x <- c(x, pillar::style_subtle(advice))

  return(x)
}


#' @importFrom pillar tbl_format_setup
#' @importFrom pillar tbl_format_header
#' @importFrom pillar tbl_format_footer
#' @method print EMP_assay_data
#' @export
print.EMP_assay_data <- function(EMPT, ..., n = NULL, width = NULL, 
                                 max_extra_cols = NULL, max_footer_lines=NULL) {
  
  . <- NULL
  
  result_num <-  length(EMPT@deposit)
  
  x <- .get.assay.EMPT(EMPT)
  
  total_nrows <-  nrow(x)
  formatted_EMPT_setup <- pillar::tbl_format_setup(x = x, width = width,n = n, 
                                           max_extra_cols = max_extra_cols, 
                                           max_footer_lines = max_footer_lines)
  
  formatted_EMPT_setup <- modify_tbl_format_setup(formatted_EMPT_setup, totalX = total_nrows)
  
  format_comment <- getFromNamespace("format_comment", "pillar")
  
  subtitle <- sprintf(" Sample=%s | Feature=%s",
                      nrow(x),
                      ncol(x)
  ) %>% 
    format_comment(width=nchar(.) + 5) %>% 
    pillar::style_subtle()
  
  header <- pillar::tbl_format_header(x, formatted_EMPT_setup) %>%
    append(subtitle, after=1)
  body <- pillar::tbl_format_body(x, formatted_EMPT_setup)
  footer <- pillar::tbl_format_footer(x, formatted_EMPT_setup)
  footer <- modify_tbl_format_footer(footer,EMPT)
  writeLines(c(header, body, footer))
  invisible(x)
}
