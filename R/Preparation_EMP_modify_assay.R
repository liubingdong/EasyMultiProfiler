#' Modify assay in condition for some special cases
#'
#' @param obj MAE or EMPT objetc
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object.
#' @param condition A character string. Set the condition data will modify.(default: '==0')
#' @param select_sample A character string. Choose which sample to modify.
#' @param select_feature A character string. Choose which feature to modify.
#' @param pseudocount A number. Set the pseudocount data will change into. (default: '0.0001')
#' @param action A character string.A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @importFrom SummarizedExperiment assay
#' @rdname EMP_modify_assay
#' @return EMPT object
#' @export
#'
#' @examples
#' # add example

EMP_modify_assay <- function(obj,experiment,condition='==0',select_sample='all',select_feature='all',pseudocount=0.0001,action='add') {
  call <- match.call()
  
  . <- NULL

  if (inherits(obj,"MultiAssayExperiment")) {
    EMPT <- .as.EMPT(obj,
                     experiment = experiment)
  }else if(inherits(obj,'EMPT')) {
    EMPT <- obj
  }
  
  replace_condition <- function(x, condition) {
    eval(parse(text = paste0("x ", condition)))
  }
  
  assay_content <- assay(EMPT) 
  
  if (all(select_sample =='all' & select_feature =='all')) {
    assay_content <- assay_content %>%
      as.data.frame() %>%
      dplyr::mutate_all( ~ ifelse(replace_condition(., condition), pseudocount, .)) %>%
      t() %>% as.data.frame() %>% # after t(), it truns into a martix, so we need  as.data.frame
      tibble::rownames_to_column('primary')
  }else if (all(select_sample =='all' & select_feature !='all')) {
    featureID <- rownames(EMPT)
    check_select_feature <- all(select_feature %in% featureID)
    
    if (!check_select_feature) {
      stop("Please make sure select_feature in the data! ")
    }

    assay_content <- assay_content %>% t() %>%
      as.data.frame() %>% 
      dplyr::mutate_if(names(.) %in% select_feature, ~ ifelse(replace_condition(., condition), pseudocount, .)) %>%
      as.data.frame() %>%
      tibble::rownames_to_column('primary')
    
  }else if (all(select_sample !='all' & select_feature =='all')) {
    sampleID <- colnames(EMPT)
    check_select_sample <- all(select_sample %in% sampleID)
    
    if (!check_select_sample) {
      stop("Please make sure select_sample in the data! ")
    }
    
    assay_content <- assay_content %>%
      as.data.frame() %>%
      dplyr::mutate_if(names(.) %in% select_sample, ~ ifelse(replace_condition(., condition), pseudocount, .)) %>%
      t() %>% as.data.frame() %>%
      tibble::rownames_to_column('primary')
  }else if (all(select_sample !='all' & select_feature !='all')) {
    
    featureID <- rownames(EMPT)
    check_select_feature <- all(select_feature %in% featureID)
    
    if (!check_select_feature) {
      stop("Please make sure select_feature in the data! ")
    }
    
    sampleID <- colnames(EMPT)
    check_select_sample <- all(select_sample %in% sampleID)
    
    if (!check_select_sample) {
      stop("Please make sure select_sample in the data! ")
    }
    
    assay_content <- assay_content %>%
      as.data.frame() %>%
      dplyr::mutate_if(names(.) %in% select_sample, ~ ifelse(replace_condition(., condition), pseudocount, .))
    
    assay_content <- assay_content %>% t() %>%
      as.data.frame() %>% 
      dplyr::mutate_if(names(.) %in% select_feature, ~ ifelse(replace_condition(., condition), pseudocount, .)) %>%
      as.data.frame() %>%
      tibble::rownames_to_column('primary')
    
  }

   .get.assay.EMPT(EMPT) <- assay_content
  
   .get.method.EMPT(EMPT) <- 'modify_assay'
   .get.algorithm.EMPT(EMPT) <- paste0(condition,' into ',pseudocount)
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