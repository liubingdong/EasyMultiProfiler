.EMP_assay_extract_EMP <- function (x,experiment,
                                    pattern_ref = 'Name',pattern=NULL,
                                    exact=F,action='add') {
  call <- match.call()
  if (inherits(x,"MultiAssayExperiment")) {
    EMPT <- .as.EMPT(x, experiment = experiment)
  }else if(inherits(x,'EMPT')) {
    EMPT <-x
  }

  EMPT %<>%.EMP_assay_extract_EMPT(pattern_ref = pattern_ref,pattern=pattern,
                                   exact=exact,action = 'add') ## action must be add here!
  if (action == 'add') {
    return(EMPT)
  } else if(action == 'get') {
    return(.get.result.EMPT(EMPT))
  }else{
    warning('action should be one of add or get!')
  }
}

.EMP_assay_extract_EMPT <- function (EMPT,
                               pattern_ref = 'Name',pattern = NULL,
                               exact=F,action = 'add') {
  feature <- primary <- NULL
  call <- match.call()

  assay_content <-.get.assay.EMPT(EMPT)
  #assay_name <- .get.assay_name.EMPT(EMPT)

  if (!is.null(pattern)) {
    .get.row_info.EMPT(EMPT) %>%
      dplyr::filter(str_detect_multi(!!dplyr::sym(pattern_ref),pattern,exact=exact)) %>%
      dplyr::pull(feature)  -> id

    id_real <- intersect(colnames(assay_content),id) ## in case that no complete matched id in the EMP_rrarefy
    assay_content %<>% dplyr::select(primary,!!id_real)
  }

  .get.assay.EMPT(EMPT) <- assay_content
  .get.info.EMPT(EMPT) <- 'EMP_assay_data'
  .get.history.EMPT(EMPT) <- call
  class(EMPT) <- 'EMP_assay_data'

  if (action == 'add') {
    return(EMPT)
  }else if(action == 'get') {
    return(.get.result.EMPT(EMPT))
  }else{
    warning('action should be one of add or get!')
  }
}



#' Extract rowdata
#'
#' @param obj EMPT or MultiAssayExperiment object.
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object. 
#' @param pattern_ref A character string. Select which column in the rowdata to extract rowdata from.
#' @param pattern A character string. Select which pattern in the feature to extract rowdata.
#' @param exact A boolean. Whether the extract use exact search method.(default:FALSE)
#' @importFrom dplyr select_if
#'
#' @return table
#' @export
#'
#' @examples
#' data(MAE)
#' ## Extract the rowdata of one existed experiment from MultiAssaayExperiment
#' MAE |>
#'   EMP_rowdata_extract('taxonomy')
#'
#' MAE |>
#'   EMP_rowdata_extract('geno_ko')  
#'   
#' ## Extract the rowdata of MultiAssaayExperiment
#' MAE |>
#'   EMP_rowdata_extract(experiment = NULL) -> total_row_data
#' dim(total_row_data)
EMP_rowdata_extract <- function(obj,experiment=NULL,pattern_ref = 'Name',pattern = NULL,exact=FALSE){

  if (inherits(obj,"MultiAssayExperiment")) {
    if (!is.null(experiment)) {
      deposit <- rowData(obj[[experiment]]) %>% as.data.frame() %>% tibble::as_tibble()
      if (!is.null(pattern)) {
        deposit %<>% dplyr::filter(str_detect_multi(!!dplyr::sym(pattern_ref),pattern,exact=exact))
      }
    }else {
      rowdata_list <- list()
      experiment_name <- names(obj)
      for (i in experiment_name) {
        rowdata_list[[i]] <- obj %>% EMP_rowdata_extract(experiment = i)
      }
      deposit <- do.call(dplyr::bind_rows, c(rowdata_list,.id='experiment_name'))
      if (!is.null(pattern)) {
        deposit %<>%
          dplyr::filter(str_detect_multi(!!dplyr::sym(pattern_ref),pattern,exact=exact)) %>%
          dplyr::select_if(~!all(is.na(.))) ## Delete columns that are all NA
      }
    }

  }else if(inherits(obj,'EMPT')) {
    EMPT <- obj
    deposit <- rowData(EMPT) %>% as.data.frame() %>% tibble::as_tibble()

    if (!is.null(pattern)) {
      deposit %<>% dplyr::filter(str_detect_multi(!!dplyr::sym(pattern_ref),pattern,exact=exact))
    }
  }else {
    stop("Please check the input data!")
  }
  return(deposit)
}

#' Extract coldata
#'
#' @param obj EMPT or MultiAssayExperiment object.
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object. 
#' @param coldata_to_assay A series of character strings. The coldata_to_assay used in action = "add". Select which columns in the coldata to transfer into assay.
#' @param assay_name A character string. The assay_name used in action = "add".(default:undefined)
#' @param action A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr any_of
#' @importFrom methods new
#'
#' @return EMPT object or table
#' @export
#'
#' @examples
#' data(MAE)
#' ## Extract the coldata/meta-data/sample-info/patient-info of one existed experiment from MultiAssaayExperiment
#' MAE |>
#'   EMP_coldata_extract('taxonomy')
#' 
#' ## Extract all coldata
#' MAE |>
#'   EMP_coldata_extract(action = 'get')  # when action = get, output is a tibble.
#' 
#' MAE |>
#'   EMP_coldata_extract(action = 'add')  # when action = add, output is a EMPT object for downstreanm analysis.

EMP_coldata_extract <- function(obj,experiment=NULL,coldata_to_assay=NULL,assay_name='undefined',action='get'){
    assay <- colname <- primary <- NULL
    call <- match.call()
    
    if (inherits(obj,'EMPT')) {
      experiment <- .get.experiment.EMPT(obj)
      coldata <- colData(obj) %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('primary') %>% 
        dplyr::arrange(primary) %>%
        tibble::as_tibble() 
    }else if (inherits(obj,'MultiAssayExperiment')){
      colData(obj) %>% as.data.frame()%>% tibble::rownames_to_column('primary') %>% tibble::as_tibble() -> coldata
      if (!is.null(experiment)) {
        real_colname <- obj[[experiment]] %>% colnames()
        sampleMap <- obj@sampleMap %>% tibble::as_tibble() %>%
          dplyr::filter(assay %in% experiment) %>%
          dplyr::select(colname,primary) %>%
          dplyr::filter(colname %in% real_colname) %>%
          dplyr::pull(primary) -> real_sample
        coldata %<>% dplyr::filter(primary %in% real_sample) %>% dplyr::select_if(~!all(is.na(.))) ## Delete any columns when all values are NA
      }
    }else {
      stop("Please check the input data!")
    }

    if (action == 'get') {
      return(coldata)
    }else if(action == 'add'){
      if (is.null(coldata_to_assay)) {
        assay_data <- coldata %>% dplyr::select(primary,where(is.numeric)) %>%
          tibble::column_to_rownames('primary') %>% t() %>% DataFrame()
      }else{
        assay_data <- coldata %>% dplyr::select(primary,any_of(!!coldata_to_assay)) %>%
          tibble::column_to_rownames('primary') %>% t() %>% DataFrame()
      }

      coldata %<>%  tibble::column_to_rownames('primary') %>% DataFrame()
      rowdata <- data.frame(feature = rownames(assay_data),Name = rownames(assay_data))

      data.se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=as.matrix(assay_data)),
                                                            rowData=rowdata, colData = coldata)

      EMPT <- new(Class = 'EMPT')

      if(is.null(experiment)){
        experiment <- 'ALL_coldata'
      }else{
        experiment <- paste0(experiment,'_coldata')
      }

      EMPT@colData <- data.se@colData
      EMPT@assays <-data.se@assays
      EMPT@NAMES <-data.se@NAMES
      EMPT@elementMetadata <-data.se@elementMetadata
      EMPT@metadata <-data.se@metadata
      EMPT@info <- 'EMP_assay_data'
      class(EMPT) <- 'EMP_assay_data'
      .get.history.EMPT(EMPT) <- call
      .get.assay_name.EMPT(EMPT) <- assay_name
      .get.experiment.EMPT(EMPT) <- experiment
      return(EMPT)
    }
}

