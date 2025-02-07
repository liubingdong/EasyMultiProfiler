perform_operation <- function(x, method,na.rm=TRUE,...) {
    result <- switch(method,
                     "mean" = mean(x,na.rm=na.rm, ...),
                     "sum" = sum(x,na.rm=na.rm, ...),
                     "median" = median(x,na.rm=na.rm, ...),
                     "max" = max(x,na.rm=na.rm, ...),
                     "min" = min(x,na.rm=na.rm, ...))
    return(result)
}

#' Multi-detect the presence/absence of a match
#'
#' @param string A character vector.
#' @param pattern Pattern to look for.
#' @param exact A boolean. Whether search the pattern matched completely.
#' @return Boolean vector
#' @export
#' @author Bingdong Liu
#' @examples
#' text <- c('Bacilli_unclassfiled','Bacteroidia_uncuture','Other')
#' str_detect_multi(text,c('Bacilli','bacteroidia'),exact=FALSE) # Ignore the capital letter
#' str_detect_multi(text,c('Bacilli','Bacteroidia'),exact=TRUE) # Set the matched completely
str_detect_multi <- function(string,pattern,exact=FALSE){
  if (length(pattern) ==1) {
    .pattern_Dectect(pattern_ref = string,pattern = pattern,exact=exact)
  }else if (length(pattern) > 1){
    id_detect <-list()
    for (j in pattern) {
      id_temp <- .pattern_Dectect(pattern_ref = string,pattern = j,exact=exact)
      id_detect[[j]] <- id_temp
    }
    .combine_logical_vectors(id_detect)
  }
}


.combine_logical_vectors <- function(vector_list) {
  combined_vec <- Reduce(`|`, vector_list)

  return(combined_vec)
}


.pattern_Dectect <- function(pattern_ref,pattern,exact=FALSE){
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
  
  if(is.null(experiment)){
    stop('Parameter experiment should be specified!')
  }

  if (!experiment %in% names(x)) {
    stop('Parameter experiment is not in the object, please check!')
  }

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

  mapping_reorder <- mapping[match(colnames(assay_data), rownames(mapping)),,drop = FALSE] # make sure the order match the colnames of assay!
  rowdata_reorder <- rowdata[match(rownames(assay_data), rowdata$feature),,drop = FALSE] # make sure the order match the colnames of assay!

  data.se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=as.matrix(assay_data)),
                      rowData=rowdata_reorder, colData = mapping_reorder)

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
  
  
  if (is(object,'MultiAssayExperiment')) {
    if (is.null(select)) {
      stop("When input data is MultiAssayExperiment, parameter select need to be specified! ")
    }else{
      data_list <- list()
      for (i in select) {
        data_list[[i]] <- object %>% .as.EMPT(i)
      }
    }
  }else if(is(object,'list')) {
    data_list <- object
  }else {
    stop("Please check the input data for as.EMP!")
  }
  
  experiment_name <- c()
  for (i in data_list) {
    if (is(i,'EMP')) {
      each_experiment_name <- .get.experiment.EMP(i)
    }else if (is(i,'EMPT')) {
      each_experiment_name <- .get.experiment.EMPT(i)
    }
    experiment_name  <- append(experiment_name,each_experiment_name)
  }
  experiment_name <- .check_duplicate(experiment_name)
  
  real_data_list <- list()
  for (i in data_list) {
    if (is(i,'EMP')) {
      each_data <- i@ExperimentList
    }else if (is(i,'EMPT')) {
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
#' @return list
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
  if (is(obj,"EMP")) {
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
  }else if(is(obj,'EMPT')) {
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

# double name for microbial tax (deprecated)
.double_tax_name <- function(df,sep=';') {

  Domain <- Kindom <- Phylum <- Class <- Order <- Family <- Genus <- Species <- Strain <- NULL

  if ('Strain' %in% colnames(df) & 'Species' %in% colnames(df)) {
    df <- df |>
      dplyr::mutate(Strain = paste0(Species,sep,Strain))
  }
  
  if ('Species' %in% colnames(df) & 'Genus' %in% colnames(df)) {
    df <- df |>
      dplyr::mutate(Species = paste0(Genus,sep,Species))
  }
  
  if ('Genus' %in% colnames(df) & 'Family' %in% colnames(df)){
    df <- df |>
      dplyr::mutate(Genus = paste0(Family,sep,Genus))
  }
  
  if ('Family' %in% colnames(df) & 'Order' %in% colnames(df)) {
    df <- df |>
      dplyr::mutate(Family = paste0(Order,sep,Family))
  }
  
  if ('Order' %in% colnames(df) & 'Class' %in% colnames(df)) {
    df <- df |>
      dplyr::mutate(Order = paste0(Class,sep,Order))
  }
    
  if ('Class' %in% colnames(df) & 'Phylum' %in% colnames(df)) {
    df <- df |>
      dplyr::mutate(Class = paste0(Phylum,sep,Class))
  }
  
  if ('Phylum' %in% colnames(df) & 'Kindom' %in% colnames(df)) {
    df <- df |>
      dplyr::mutate(Phylum = paste0(Kindom,sep,Phylum))
  }
  
  if ('Kindom' %in% colnames(df) & 'Domain' %in% colnames(df)) {
    df <- df |>
      dplyr::mutate(Kindom = paste0(Domain,sep,Kindom))
  }
  
  return(df) 
}

# full name for microbial tax 
.tax_to_full <- function(df,sep=';') {
  
  Domain <- Kindom <- Phylum <- Class <- Order <- Family <- Genus <- Species <- Strain <- NULL
  
  total_tax_info <- c('Domain','Kindom','Phylum','Class','Order','Family','Genus','Species','Strain')

  raw_col <- colnames(df)
  
  if ('Strain' %in% colnames(df)) {
    need_info <- total_tax_info[1: match('Strain', total_tax_info)]
    df <- df |>
      tidyr::unite(col = "Strain", dplyr::any_of(need_info), sep = sep, remove = FALSE)
  }
  if ('Species' %in% colnames(df)) {
    need_info <- total_tax_info[1: match('Species', total_tax_info)]
    df <- df |>
      tidyr::unite(col = "Species", dplyr::any_of(need_info), sep = sep, remove = FALSE)
  }  
  if ('Genus' %in% colnames(df)) {
    need_info <- total_tax_info[1: match('Genus', total_tax_info)]
    df <- df |>
      tidyr::unite(col = "Genus", dplyr::any_of(need_info), sep = sep, remove = FALSE)
  }  
  if ('Family' %in% colnames(df)) {
    need_info <- total_tax_info[1: match('Family', total_tax_info)]
    df <- df |>
      tidyr::unite(col = "Family", dplyr::any_of(need_info), sep = sep, remove = FALSE)
  }
  if ('Order' %in% colnames(df)) {
    need_info <- total_tax_info[1: match('Order', total_tax_info)]
    df <- df |>
      tidyr::unite(col = "Order", dplyr::any_of(need_info), sep = sep, remove = FALSE)
  }
  if ('Class' %in% colnames(df)) {
    need_info <- total_tax_info[1: match('Class', total_tax_info)]
    df <- df |>
      tidyr::unite(col = "Class", dplyr::any_of(need_info), sep = sep, remove = FALSE)
  }
  if ('Phylum' %in% colnames(df)) {
    need_info <- total_tax_info[1: match('Phylum', total_tax_info)]
    df <- df |>
      tidyr::unite(col = "Phylum", dplyr::any_of(need_info), sep = sep, remove = FALSE)
  }
  if ('Kindom' %in% colnames(df)) {
    need_info <- total_tax_info[1: match('Kindom', total_tax_info)]
    df <- df |>
      tidyr::unite(col = "Kindom", dplyr::any_of(need_info), sep = sep, remove = FALSE)
  }

  df <- df %>% dplyr::select(dplyr::all_of(raw_col))
  return(df) 
}



# single name for microbial tax 
.tax_to_single <- function(df,sep=';') {
  
  raw_col <- colnames(df)
  need_col <- raw_col[-1:-2] # del feature and first level
  
  pattern <- paste0("(?<=", sep, ")[^", sep, "]*$")
  for (i in need_col) {
    df <- df %>% 
      dplyr::mutate(!!i := stringr::str_extract(df[[i]], pattern))
  }

  df <- df %>% dplyr::select(dplyr::all_of(raw_col))
  return(df) 
}


# Check the EMPT from microbial data
#' @importFrom SummarizedExperiment rowData
.check_is_tax <- function(EMPT,skip_old_feature=FALSE){
  row_names <- rowData(EMPT) |> colnames()
  if (skip_old_feature) {
      total_tax_name <- c('feature','Kindom','Phylum','Class','Order','Family','Genus','Species','Strain')
  }else{
      total_tax_name <- c('feature','old_feature','Kindom','Phylum','Class','Order','Family','Genus','Species','Strain')
  }

  if (all(row_names %in% total_tax_name)) {
    flag <- TRUE
  }else{
    flag <- FALSE
  }
  return(flag)
}

# Check the tax anotation is full or single
#' @importFrom SummarizedExperiment rowData
#' @importFrom stringr str_count
.check_is_tax_full <- function (EMPT,sep=';') {
  total_tax_name <- c('feature','Kindom','Phylum','Class','Order','Family','Genus','Species','Strain')
  row_data <- rowData(EMPT) 
  row_name <- row_data|> colnames()
  # check the two example in the last column to confirm the tax anotation method
  flag <- all(str_count(row_data[1,length(row_name)],sep) == which(total_tax_name %in%  row_name[length(row_name)]) -2,
              str_count(row_data[2,length(row_name)],sep) == which(total_tax_name %in%  row_name[length(row_name)]) -2)
  return(flag)
}


# Check the confilct in the collpase for microbial data 
#' @importFrom SummarizedExperiment rowData
.check_tax_collapse_confilct <- function(EMPT,estimate_group){
  
  EMPT_full <- EMPT |> EMP_feature_convert(from = 'tax_single',add = 'tax_full')
  rowdata_single <- rowData(EMPT) 
  rowdata_full <- rowData(EMPT_full)
  
  flag <- rowdata_single[[estimate_group]] |> unique() |> length() != rowdata_full[[estimate_group]] |> unique() |> length()
  return(flag)
}

# make choice for users to deal with the confilct in the collpase
#' @importFrom utils menu
.choose_tax_anotation <- function() {
  selection <- c("Keep single-level tax anotation.", "Add full-level tax anotation.")
  choice <- menu(selection, title = "Detected the confilcts between single-level and full-level annotation:")
  
  if (choice == 0) {
    return("single")
  }else if (choice == 1) {
    return("single")
  }else {
    return("full")
  }
}



#' Transfer microbial data from EasyMultiProfiler to EasyMciroPlot
#'
#' @param obj EMPT object.
#' @param estimate_group A character string. Select the column in the coldata to add in the mapping file.
#' @rdname EMP_to_EMP1
#' @return list
#' @export
#' @examples
#' data(MAE)
#' # Get the data from EasyMultiProfiler
#' MAE |>
#'   EMP_assay_extract('taxonomy') |>
#'   EMP_feature_convert(from = 'tax_single',add = 'tax_full') |>
#'   EMP_to_EMP1(estimate_group = 'Group') -> deposit
#' \dontrun{ 
#' # Work in the EasyMicroPlot
#' library(EasyMicroPlot)
#' alpha_re <- alpha_plot(data = deposit$data,
#'                        design = deposit$mapping,
#'                        min_relative = 0.001,min_ratio = 0.7,method = 'ttest')
#' }
#' @importFrom stringr str_detect
EMP_to_EMP1 <- function(obj,estimate_group){
  
  primary <- feature <- NULL

  deposit <- list()
  
  if (!is(obj,"EMPT")) {
    stop('Please input EMPT object!')
  }else if(is(obj,'EMPT')) {
    EMPT <- obj
  }
  
  rowdata <- EMPT |>
    EMP_rowdata_extract() 
  check_feature <- str_detect(rowdata$feature,';') |> all()
  if(!check_feature) {
    stop("This function need full taxonomy with ';', please check input data and EMP_feature_convert!")
  }
  
  if ('old_feature' %in% colnames(rowdata)) {
    stop("This function only work before EMP_collapse!")
  }
  
  mapping <- EMPT |>
    EMP_coldata_extract() |>
    dplyr::select(primary,{{estimate_group}}) |>
    dplyr::rename(SampleID = primary,Group={{estimate_group}}) |> as.data.frame()
  
  meta_data <-  EMPT |>
    EMP_coldata_extract() |>
    dplyr::select(primary,where(is.numeric)) |>
    dplyr::rename(SampleID = primary) |> as.data.frame()
  
  
  tax_level <- obj |> EMP_rowdata_extract() |> dplyr::select(-feature) |> colnames()
  tax_ava <- c('Phylum','Class','Order','Family','Genus','Species')
  real_tax <- intersect(tax_level,tax_ava)
  
  for (tax_select in real_tax) {
    assay_df <- obj |> 
      EMP_collapse(estimate_group = tax_select,method = 'sum',collapse_by = 'row',action = 'get') |>
      tibble::column_to_rownames('primary') |>
      t() |> as.data.frame() |> suppressMessages() |>
      tibble::rownames_to_column('SampleID') 
    deposit[['data']][[tax_select]] <- assay_df
  }
  
  deposit[['mapping']] <- mapping
  deposit[['meta_data']] <- meta_data
  
  return(deposit)
}

## Enhance the print for EMP_assay_data
## These code below is modified from MicrobiotaProcess
modify_tbl_format_setup <- function(x, totalX,totalY,by,...){ 
  tmpxx <- x$tbl_sum %>%
    strsplit(" ") %>%
    unlist()
  tmpxx <- paste(getFromNamespace("big_mark", "pillar")(totalX), tmpxx[2], getFromNamespace("big_mark", "pillar")(totalY))
  names(tmpxx) <- c("A EMPT-tibble (EMPT object) size")
  x$tbl_sum <- tmpxx
  x$rows_total <- totalX
  x$cols_total <- totalY
  if (by=='primary') {
    x$rows_missing <- totalX - nrow(x$df)
  }else if (by=='feature') {
    x$rows_missing <- totalY - nrow(x$df)
  }else{
    "please check the by in modify_tbl_format_setup!"
  }
  return(x)
}

modify_tbl_format_footer <- function(x,EMPT,...){
    result_num <-  length(EMPT@deposit)
    result_info <- paste0("The obeject contains ", result_num," analysis result.")
    x <- c(x, pillar::style_subtle(result_info))

  return(x)
}

#' @importFrom SummarizedExperiment assay
#' @importFrom pillar tbl_format_setup
#' @importFrom pillar tbl_format_header
#' @importFrom pillar tbl_format_footer
enhance_print <- function(EMPT, ..., n = NULL, width = NULL, 
                                 max_extra_cols = NULL, max_footer_lines=NULL) {
  
  . <- assay_dim <- NULL
  
  result_num <-  length(EMPT@deposit)
  
  # Here is for speed up
  info <- .get.info.EMPT(EMPT)
  if (info %in% c('EMP_assay_data','EMP_rrarefy','EMP_decostand')) {
      assay_dim <- EMPT %>% dim() %>% rev()
      if (any(assay_dim > 500)) {

        select_col <- ifelse(assay_dim[1]<20,assay_dim[1],20)
        select_row <- ifelse(assay_dim[2]<20,assay_dim[2],20)
        
        x <- assay(EMPT)[1:select_row,1:select_col] %>% t() %>% as.data.frame() %>%
          tibble::rownames_to_column('primary') %>% tibble::as_tibble()
      }else{
        x <- .get.assay.EMPT(EMPT)
      }
  }else{
      x <- .get.result.EMPT(EMPT) %>% suppressMessages()
      assay_dim <- EMPT %>% dim() %>% rev()  
  }
  check_flag <- colnames(x)[1]
  total_nrows <-  assay_dim[1]
  total_cols <-  assay_dim[2]
  formatted_EMPT_setup <- pillar::tbl_format_setup(x = x, width = width,n = n, 
                                           max_extra_cols = max_extra_cols, 
                                           max_footer_lines = max_footer_lines)
  
  formatted_EMPT_setup <- modify_tbl_format_setup(formatted_EMPT_setup,by=check_flag,
                                           totalX = total_nrows,totalY=total_cols)
  
  format_comment <- getFromNamespace("format_comment", "pillar")
  
  subtitle <- sprintf(" Sample=%s | Feature=%s",
                      assay_dim[1],
                      (assay_dim[2])
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


#' Save result into environment 
#' @param data data. 
#' @param var_name  a name to be assigned to data.
#' @param envir the \link{environment} to use.
#' @param ... Addtional parameters, see also \code{\link[base]{assign}}.
#' @rdname EMP_save_var
#' @return EMPT object
#' @export
#' @examples
#' data(MAE)
#' MAE |> 
#'   EMP_assay_extract('host_gene',pattern = 'A1BG',pattern_ref = 'feature') |>
#'   EMP_save_var('temp_data') |>
#'   EMP_collapse(estimate_group = 'Group',collapse_by = 'col') |>
#'   EMP_heatmap_plot()
#' 
#' temp_data
EMP_save_var <- function(data, var_name,envir = .GlobalEnv,...) {
  assign(var_name, data, envir = envir,...)  
  return(data) 
}

#' Info output with style and color
#' @param str data. 
#' @param order order.
#' @param color color.
#' @param pure logical. Default message, warning and stop.
#' @param show character including message, warning and stop.
#' @param call. logical, indicating if the call should become part of the warning or error message.
#' @param ... See also \code{\link[base]{message}}, \code{\link[base]{warning}} and \code{\link[base]{stop}}.
#' @export
#' @author Bingdong Liu
#' @section Detaild :
#' This function performed followed by the ANSI.
#'
#' For color example:
#'
#' 30: Black
#'
#' 31: Red
#'
#' 32: Green
#'
#' 33: Yellow
#'
#' 34: Blue
#'
#' 35: Magenta
#'
#' 36: Cyan
#'
#' 37: White
#'
#'
#' For order example:
#'
#' 0: Regular
#'
#' 1: Bold
#'
#' 2: blurred
#'
#' 3: Italic
#'
#' 4: Underline
#'
#' 5: Blinking text
#'
#' 9: Strikethrough
#'
#' 40: Black background
#'
#' 41: Red background
#'
#' 42: Green background
#'
#' 43: Yellow background
#'
#' 44: Blue background
#'
#' 45: Magenta background
#'
#' 46: Cyan background
#'
#' 47: White background
#'
#' @examples
#' \dontrun{
#' # blue
#' EMP_message('Hello World',color = 94,order = 0)
#' # green 
#' EMP_message('Hello World',color = 32,order = 0)
#' # green and bold
#' EMP_message('Hello World',color = 32,order = 1)
#' # background  
#' EMP_message('Hello World',color = 31,order = 7)
#' # underline 
#' EMP_message('Hello World',color = 32,order = 4)
#' # italics 
#' EMP_message('Hello World',color = 32,order = 3)
#' # blurred
#' EMP_message('Hello World',color = 32,order = 2)
#' 
#' # warning with background  
#' EMP_message('Hello World',color = 31,order = 7,show='warning')
#' 
#' # stop with background  
#' EMP_message('Hello World',color = 31,order = 7,show='stop')
#' }
EMP_message <- function(str,order=0,color=32,pure=FALSE,show='message',call.=FALSE,...) {
  stopifnot(is.logical(pure))
  show <- match.arg(show,c('message','stop','warning'))
  if (pure) {
    info <- str
  }else{
    info <- paste0("\033[",order,";",color,"m",str,"\033[0m")
  }
  switch (show,
    'message' = {message(info,...)},
    'warning' = {warning(info,call. = call.,...)},
    'stop' = {stop(info,call. = call.,...)}
  )
}



#' Detect the top or bottom number
#' @param x vec. 
#' @param n This parameter requires a number greater than 0. If it is greater than 1, selection is based on the numeric value. If it is between 0 and 1, selection is based on a percentage.
#' @param type top or bottom.
#' @param index output the index or the number.
#' @export
#' @author Bingdong Liu
#' @examples
#' top_detect(1:10,n = 3,index = FALSE)
#' 
#' top_detect(1:10,n = 3,index = TRUE)
#' 
#' data(MAE)
#'
#' MAE |>
#'   EMP_assay_extract(experiment='taxonomy')|> 
#'   EMP_alpha_analysis() |>
#'   EMP_filter(top_detect(shannon,3))
top_detect<- function(x, n, type = 'top',index=TRUE) {
  type <- match.arg(type, c("top", "bottom"))  
  
  if (length(x) < n) {
    stop("The input n must not exceed the scale of the data!")
  }
  
  # 确保 n 合法（大于 0 且不超过向量长度）
  if (n <= 0) {
    stop("n must be greater than 0.")
  } else if (n < 1) {
    n <- ceiling(length(x) * n)  
  } else {
    n <- min(n, length(unique(x)))  
  }
  
  if (type == "top") {
    threshold <- sort(unique(x), decreasing = TRUE)[1:n]  
  } else {
    threshold <- sort(unique(x))[1:n]  
  }

  if (index) {
    result <- x %in% threshold
    return(result)
  }else{
    return(threshold)
  }

}
