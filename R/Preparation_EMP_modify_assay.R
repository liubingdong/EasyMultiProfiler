#' Modify assay in condition for some special cases
#'
#' @param obj MAE or EMPT objetc
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object.
#' @param condition A character string. Set the condition data will modify.(default: '==0')
#' @param select_sample A character string. Choose which sample to modify.
#' @param select_feature A character string. Choose which feature to modify.
#' @param pseudocount A number or function. Set the pseudocount data will change into. (default: '0.0001')
#' @param action A character string.A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @importFrom SummarizedExperiment assay
#' @rdname EMP_modify_assay
#' @return EMPT object
#' @export
#'
#' @examples
#' data(MAE)
#' ## For some special cases, to modify the assay data and EMP_identify_assay works better in most cases.
#' ## Change the expression value which is 0 into 0.0001
#' MAE |>
#'   EMP_assay_extract('geno_ec') |>
#'   EMP_modify_assay('==0',pseudocount=0.0001)
#' 
#' ## Change the counts which is below 
#' MAE |>
#'   EMP_assay_extract('taxonomy') |>
#'   EMP_modify_assay('<10',pseudocount=0) 
#' 
#' ## suport external function
#' MAE |>
#'   EMP_assay_extract('taxonomy') |>
#'   EMP_modify_assay('==0',pseudocount=function(x) x +1 ) 
#' 
#' MAE |>
#'   EMP_assay_extract('taxonomy') |>
#'   EMP_modify_assay('>0',pseudocount=function(x) log(x) ) 

EMP_modify_assay <- function(obj,condition='==0',experiment,select_sample='all',select_feature='all',pseudocount=0.0001,action='add') {
  call <- match.call()
  
  . <- NULL

  if (is(obj,"MultiAssayExperiment")) {
    EMPT <- .as.EMPT(obj,
                     experiment = experiment)
  }else if(is(obj,'EMPT')) {
    EMPT <- obj
  }
  
  replace_condition <- function(x, condition) {
    eval(parse(text = paste0("x ", condition)))
  }
  
  if(is(pseudocount,'numeric') | is(pseudocount,'integer')) {
    check_funtion <- FALSE
  }else if (is(pseudocount,'function')){
    check_funtion <- TRUE
  }else{
    stop("Parameter pseudocount should be a number or function!")
  }

  assay_content <- assay(EMPT) 
  
  if (all(select_sample =='all' & select_feature =='all')) {
    assay_content <- assay_content %>%
      as.data.frame() %>%
      dplyr::mutate_all( ~ ifelse(replace_condition(., condition), if (is.function(pseudocount)) pseudocount(.) else pseudocount, .)) %>%
      t() %>% as.data.frame() %>% # after t(), it truns into a martix, so we need  as.data.frame
      tibble::rownames_to_column('primary')
  }else if (all(select_sample =='all' & select_feature !='all')) {
    featureID <- rownames(EMPT)
    check_select_feature <- all(select_feature %in% featureID)
    
    if (!check_select_feature) {
      stop("Please make sure select_feature in the data!")
    }

    assay_content <- assay_content %>% t() %>%
      as.data.frame() %>% 
      dplyr::mutate_if(names(.) %in% select_feature, ~ ifelse(replace_condition(., condition), if (is.function(pseudocount)) pseudocount(.) else pseudocount, .)) %>%
      as.data.frame() %>%
      tibble::rownames_to_column('primary')
    
  }else if (all(select_sample !='all' & select_feature =='all')) {
    sampleID <- colnames(EMPT)
    check_select_sample <- all(select_sample %in% sampleID)
    
    if (!check_select_sample) {
      stop("Please make sure select_sample in the data!")
    }
    
    assay_content <- assay_content %>%
      as.data.frame() %>%
      dplyr::mutate_if(names(.) %in% select_sample, ~ ifelse(replace_condition(., condition), if (is.function(pseudocount)) pseudocount(.) else pseudocount, .)) %>%
      t() %>% as.data.frame() %>%
      tibble::rownames_to_column('primary')
  }else if (all(select_sample !='all' & select_feature !='all')) {
    
    featureID <- rownames(EMPT)
    check_select_feature <- all(select_feature %in% featureID)
    
    if (!check_select_feature) {
      stop("Please make sure select_feature in the data!")
    }
    
    sampleID <- colnames(EMPT)
    check_select_sample <- all(select_sample %in% sampleID)
    
    if (!check_select_sample) {
      stop("Please make sure select_sample in the data!")
    }
    
    assay_content <- assay_content %>%
      as.data.frame() %>%
      dplyr::mutate_if(names(.) %in% select_sample, ~ ifelse(replace_condition(., condition), if (is.function(pseudocount)) pseudocount(.) else pseudocount, .))
    
    assay_content <- assay_content %>% t() %>%
      as.data.frame() %>% 
      dplyr::mutate_if(names(.) %in% select_feature, ~ ifelse(replace_condition(., condition), if (is.function(pseudocount)) pseudocount(.) else pseudocount, .)) %>%
      as.data.frame() %>%
      tibble::rownames_to_column('primary')
    
  }

   .get.assay.EMPT(EMPT) <- assay_content
  
   .get.method.EMPT(EMPT) <- 'modify_assay'
   if (check_funtion) {
    .get.algorithm.EMPT(EMPT) <- paste(deparse(pseudocount), collapse = " ")
   }else{
    .get.algorithm.EMPT(EMPT) <- paste0(condition,' into ',pseudocount)
   }
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