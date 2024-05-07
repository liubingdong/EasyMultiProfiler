#' @importFrom tidybulk tidybulk
# ' @importFrom dplyr select
#' @importFrom dplyr rename
#' @importFrom tidyr drop_na
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @noRd 
EMP_collapse_byrow <- function(x,experiment,estimate_group=NULL,method='sum',na_string=c('NA','null',''),
    collapse_sep=' ',action='add',...) {
  `.sample` <- counts <- feature <- primary <- old_feature <- NULL 
  call <- match.call()
  if (inherits(x,"MultiAssayExperiment")) {
    EMPT <- .as.EMPT(x,
                     experiment = experiment)
    estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)
    .get.estimate_group.EMPT(EMPT) <- estimate_group
    .get.experiment.EMPT(EMPT) <- experiment
    class(EMPT) <- 'EMP_assay_data'
  }else if(inherits(x,'EMPT')){
    EMPT <- x
    estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)
    experiment <- .get.experiment.EMPT(EMPT)
    class(EMPT) <- 'EMP_assay_data'
  }else {
    stop('Please check the input data')
  }
  
  # check the df attr
  df_attr_row <- EMPT@deposit2[['df_attr_row']]
  if (!is.null(df_attr_row)) {
    if (df_attr_row$feature == estimate_group) {
      stop('estimate_group parameter should be different in a serious of EMP_collapse_byrow! ')
    }
  }

  new_assay_data <- EMPT %>% 
    tidybulk::tidybulk() %>%
    dplyr::select(.sample,counts,!!estimate_group) %>%
    dplyr::rename(primary=.sample,feature=!!estimate_group) %>%
    tidyr::drop_na() %>%
    dplyr::filter(!feature %in% !!na_string) %>% 
    dplyr::group_by(primary, feature) %>%
    dplyr::summarise(counts = .perform_operation(counts,method),.groups='drop') %>% 
    tidyr::pivot_wider(names_from = 'feature',values_from = 'counts') %>%
    tibble::column_to_rownames('primary') %>% t()
  

  if (is.null(df_attr_row)) {
    new_row_data <- .get.row_info.EMPT(EMPT)  %>% 
      tidyr::drop_na(!!estimate_group) %>%
      dplyr::filter(!(!!dplyr::sym(estimate_group) %in% !!na_string)) %>%  # filter the missing value
      .collpseBygroup.tibble(estimate_group = estimate_group,method=method,collapse_by='row',collapse_sep=collapse_sep,...)
    EMPT@deposit2[['df_attr_row']] <- lapply(new_row_data, attr, "raw_info")%>% as.data.frame()
  }else {
    row_data <- .get.row_info.EMPT(EMPT)
    
    ## Recovery the raw name
    change_name <- EMPT@deposit2[['df_attr_row']] %>% dplyr::select(feature,old_feature)
    feature_origin <- change_name %>% dplyr::pull(feature)
    old_feature_origin <- change_name %>% dplyr::pull(old_feature)
    row_data <- row_data %>% dplyr::rename({{feature_origin}} := feature,{{old_feature_origin}} := old_feature)

    new_row_data <- row_data  %>% 
      tidyr::drop_na(!!estimate_group) %>%
      dplyr::filter(!(!!dplyr::sym(estimate_group) %in% !!na_string)) %>%  # filter the missing value
      .collpseBygroup.tibble(estimate_group = estimate_group,method=method,collapse_by='row',collapse_sep=collapse_sep,...)
    EMPT@deposit2[['df_attr_row']] <- lapply(new_row_data, attr, "raw_info")%>% as.data.frame()
  }
  
  col_data <- colData(EMPT)
  
  message_info <- list()
  message_info %<>% append(paste0('Feature changed!'))
  message_info %<>% append(paste0('Current feature: ',estimate_group))
  
  data.se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=as.matrix(new_assay_data)),
                                                       rowData=new_row_data, colData = col_data)
  
  EMPT@colData <- data.se@colData
  EMPT@assays <-data.se@assays
  EMPT@NAMES <-data.se@NAMES
  EMPT@elementMetadata <-data.se@elementMetadata
  EMPT@metadata <-data.se@metadata

  
  .get.estimate_group.EMPT(EMPT) <- estimate_group
  .get.message_info.EMPT(EMPT) <- message_info
  .get.history.EMPT(EMPT) <- call
  .get.method.EMPT(EMPT) <- 'collapse'
  .get.algorithm.EMPT(EMPT) <- 'collapse_byrow'
  .get.info.EMPT(EMPT) <- 'EMP_assay_data'
  
  if (action == 'add') {
    return(EMPT)
  }else if(action == 'get') {
    return(.get.result.EMPT(EMPT))
  }else{
    warning('action should be one of add or get!')
  }
}

.EMP_collapse_byrow_m <- memoise::memoise(EMP_collapse_byrow)


#' @importFrom tidybulk tidybulk
#' @importFrom dplyr rename
#' @importFrom tidyr drop_na
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @return xx object
#' @noRd 
EMP_collapse_bycol <- function(x,experiment,estimate_group=NULL,method='sum',na_string=c('NA','null',''),collapse_sep=' ',action='add',...) {
  `.feature` <- counts <- primary <- feature <- old_feature <- NULL
  call <- match.call()
  if (inherits(x,"MultiAssayExperiment")) {
    EMPT <- .as.EMPT(x,
                     experiment = experiment)
    estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)
    .get.estimate_group.EMPT(EMPT) <- estimate_group
    .get.experiment.EMPT(EMPT) <- experiment
    class(EMPT) <- 'EMP_assay_data'
  }else if(inherits(x,'EMPT')){
    EMPT <- x
    estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)
    experiment <- .get.experiment.EMPT(EMPT)
    class(EMPT) <- 'EMP_assay_data'
  }else {
    stop('Please check the input data')
  }
  
  # check the df attr
  df_attr_col <- EMPT@deposit2[['df_attr_col']]
  if (!is.null(df_attr_col)) {
    if (df_attr_col$primary == estimate_group) {
      stop('estimate_group parameter should be different in a serious of EMP_collapse_byrow! ')
    }
  }
  
  new_assay_data <- EMPT %>% 
    tidybulk::tidybulk() %>%
    dplyr::select(.feature,counts,!!estimate_group) %>%
    dplyr::rename(primary=!!estimate_group,feature=.feature) %>%
    tidyr::drop_na() %>%
    dplyr::filter(!primary %in% !!na_string) %>% 
    dplyr::group_by(primary, feature) %>%
    dplyr::summarise(counts = .perform_operation(counts,method),.groups='drop') %>%
    tidyr::pivot_wider(names_from = 'feature',values_from = 'counts') %>%
    tibble::column_to_rownames('primary') %>% t()
  
  
  if (is.null(df_attr_col)) {
    new_col_data <- .get.mapping.EMPT(EMPT)  %>% 
      tidyr::drop_na(!!estimate_group) %>%
      dplyr::filter(!(!!dplyr::sym(estimate_group) %in% !!na_string)) %>%  # filter the missing value
      .collpseBygroup.tibble(estimate_group = estimate_group,collapse_sep=collapse_sep,collapse_by = 'col',method=method,...)
    EMPT@deposit2[['df_attr_col']] <- lapply(new_col_data, attr, "raw_info")%>% as.data.frame()
    
  }else {
    col_data <- .get.mapping.EMPT(EMPT)

    ## Recovery the raw name
    change_name <- EMPT@deposit2[['df_attr_row']] %>% dplyr::select(feature,old_feature)
    feature_origin <- change_name %>% dplyr::pull(feature)
    old_feature_origin <- change_name %>% dplyr::pull(old_feature)
    row_data <- row_data %>% dplyr::rename({{feature_origin}} := feature,{{old_feature_origin}} := old_feature)

    new_col_data <- col_data  %>% 
      tidyr::drop_na(!!estimate_group) %>%
      dplyr::filter(!(!!dplyr::sym(estimate_group) %in% !!na_string)) %>%  # filter the missing value
      .collpseBygroup.tibble(estimate_group = estimate_group,collapse_sep=collapse_sep,collapse_by = 'col',method=method,...)
    EMPT@deposit2[['df_attr_col']] <- lapply(new_col_data, attr, "raw_info")%>% as.data.frame()
  }
  
  row_data <- EMPT %>% EMP_rowdata_extract()
  new_col_data %<>% tibble::column_to_rownames('primary')
    
  message_info <- list()
  message_info %<>% append(paste0('Primary changed!'))
  message_info %<>% append(paste0('Current primary: ',estimate_group))
  
  data.se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=as.matrix(new_assay_data)),
                                                        rowData=row_data, colData = new_col_data)
  
  
  EMPT@colData <- data.se@colData
  EMPT@assays <-data.se@assays
  EMPT@NAMES <-data.se@NAMES
  EMPT@elementMetadata <-data.se@elementMetadata
  EMPT@metadata <-data.se@metadata
  
  
  .get.estimate_group.EMPT(EMPT) <- estimate_group
  .get.message_info.EMPT(EMPT) <- message_info
  .get.history.EMPT(EMPT) <- call
  .get.method.EMPT(EMPT) <- 'collapse'
  .get.algorithm.EMPT(EMPT) <- 'collapse_bycol'
  .get.info.EMPT(EMPT) <- 'EMP_assay_data'
  
  if (action == 'add') {
    return(EMPT)
  }else if(action == 'get') {
    return(.get.result.EMPT(EMPT))
  }else{
    warning('action should be one of add or get!')
  }
}


.EMP_collapse_bycol_m <- memoise::memoise(EMP_collapse_bycol)


#' @importFrom purrr reduce
#' @noRd
.collpseBygroup.tibble <- function(df,estimate_group,method,collapse_sep=' ',collapse_by,...) {
  feature <- primary <- NULL
  idx <- df  %>% dplyr::select(-!!estimate_group) %>% colnames()
  old_col_name <- df %>% colnames()
  data_deposit <- list()


  idx_numeric <- df  %>% dplyr::select(-!!estimate_group) %>% dplyr::select(where(is.numeric)) %>% colnames()
  idx_character <- df  %>% dplyr::select(-!!estimate_group) %>% dplyr::select(where(is.character)) %>% colnames()
  
  for (i in idx_character) {
    df %>%
      dplyr::group_by(!!dplyr::sym(estimate_group)) %>%
      dplyr::summarise(!!i := paste0(unique(!!dplyr::sym(i)), collapse = collapse_sep ) ) -> data_deposit[[i]]
  }

  for (i in idx_numeric) {
    df %>%
      dplyr::group_by(!!dplyr::sym(estimate_group)) %>%
      dplyr::summarise(!!i := .perform_operation(!!dplyr::sym(i),method,...),.groups='drop') -> data_deposit[[i]]
  }

  combined_df <- purrr::reduce(data_deposit, dplyr::full_join, by = estimate_group)
  
  for (i in old_col_name) {
    attr(combined_df[[i]], "raw_info")  <- i
  }
  
  if (collapse_by == 'row') {
      combined_df %<>% dplyr::rename(old_feature = feature,
                                 feature = !!dplyr::sym(estimate_group))
  }else if (collapse_by == 'col') {
      combined_df %<>% dplyr::rename(old_primary = primary,
                                 primary = !!dplyr::sym(estimate_group))
  }

  
  return(combined_df)
}



#' Aggregates abundace or experssion from the same attributes in the coldata or rowdata
#'
#' @param x Object in EMPT or MultiAssayExperiment format.
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object.
#' @param estimate_group A character string. Select the column in the rowdata or coldata to collapse. 
#' @param method  A character string. Methods include mean, sum, median, min, max.
#' @param na_string A series of character strings. Indicate which characters can be considered missing values.
#' @param collapse_sep A character string. The linking symbol used when strings are combined.
#' @param action A character string. A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param collapse_by A character string. Methods include col or row.
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
#' @param ... Further parameters passed to the function mean, sum, median, min, max in the base package.
#' @importFrom tidybulk tidybulk
#' @importFrom dplyr rename
#' @importFrom tidyr drop_na
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @return EMPT object
#' @export
#'
#' @examples
#' data(MAE)
#' ##  merge assay data accoding to duplicate coldata.
#' MAE |>
#'   EMP_collapse(experiment = 'untarget_metabol',collapse_by='col',
#'                estimate_group = 'Group',method = 'mean',collapse_sep = '+') 
#' ##  merge assay data accoding to duplicate rowdata.
#' MAE |> EMP_rowdata_extract('untarget_metabol')
#' MAE |>
#'   EMP_collapse(experiment = 'untarget_metabol',collapse_by='row',na_string = c("NA", "null", "","-"),
#'                estimate_group = 'MS2kegg',method = 'mean',collapse_sep = '+') 
#' ## combie two collapse method
#' MAE |>
#'   EMP_collapse(experiment = 'untarget_metabol',collapse_by='row',na_string = c("NA", "null", "","-"),
#'                estimate_group = 'MS2kegg',method = 'mean',collapse_sep = '+') |>
#'   EMP_collapse(collapse_by='col',estimate_group = 'Group',method = 'mean',collapse_sep = '+')
EMP_collapse <- function (x,experiment=NULL,estimate_group=NULL,method='sum',na_string=c('NA','null',''),collapse_by,collapse_sep=' ',action='add',use_cached=TRUE,...) {

  if (use_cached == FALSE) {
    memoise::forget(.EMP_collapse_byrow_m) %>% invisible()
    memoise::forget(.EMP_collapse_bycol_m) %>% invisible()
  }  

  if (collapse_by == 'row') {
    deposit <- .EMP_collapse_byrow_m(x,experiment,
                                  estimate_group,method,
                                  na_string,collapse_sep,action,...)
  }else if (collapse_by == 'col') {
    deposit <- .EMP_collapse_bycol_m(x,experiment,
                                  estimate_group,method,
                                  na_string,collapse_sep,action,...)
  }else{
    stop("Please set parameter collapse_by (row or col)! ")
  }
  return(deposit)
}







