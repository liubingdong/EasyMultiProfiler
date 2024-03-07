#' @importFrom vegan decostand
#' @noRd
.EMP_decostand <- function(EMPT,bySample=T,logbase,pseudocount=pseudocount,...){

  MARGIN <- ifelse(bySample, 2, 1)

  .get.assay.EMPT(EMPT) -> assay_data

  .get.method.EMPT(EMPT) -> method

  if (method=="log"){
    method_name <- paste0(method, logbase)
  }else{
    method_name <- method
  }


  if (method == 'relative') {
    method  <- 'total'
  }

  ## .calc_rclr unused argument pseudocount
  if (method == 'rclr') {
    assay_decostand_data <- assay_data %>%  tibble::column_to_rownames('primary') %>% t() %>%
    vegan::decostand(., method = method, MARGIN = MARGIN,logbase=logbase,...)
  }else if(method == 'integer'){
    assay_decostand_data <- assay_data %>%  tibble::column_to_rownames('primary') %>% t() %>%
     round(digits = 0)
  }else {
    assay_decostand_data <- assay_data %>%  tibble::column_to_rownames('primary') %>% t() %>%
    vegan::decostand(., method = method, MARGIN = MARGIN,logbase=logbase,pseudocount=pseudocount,...)
  }


  .get.assay.EMPT(EMPT) <- assay_decostand_data %>% t() %>% as.data.frame() %>% tibble::rownames_to_column('primary') %>% tibble::as_tibble() %>% dplyr::select(primary,everything())
  .get.method.EMPT(EMPT) <- method_name
  .get.assay_name.EMPT(EMPT) <- method_name
  .get.algorithm.EMPT(EMPT) <- 'EMP_decostand'
  .get.info.EMPT(EMPT) <- 'EMP_decostand'
  EMPT

}

.EMP_decostand_m <- memoise::memoise(.EMP_decostand)

#' Title
#'
#' @param x wait_for_add
#' @param experiment wait_for_add
#' @param method wait_for_add
#' @param bySample wait_for_add
#' @param logbase wait_for_add
#' @param use_cached wait_for_add
#' @param pseudocount wait_for_add
#' @param action wait_for_add
#' @param ... wait_for_add
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_decostand <- function(x,experiment,method,bySample=T,logbase =2,use_cached = T,pseudocount=0.0000001,action='add',...){
  call <- match.call()
  if (inherits(x,"MultiAssayExperiment")) {
    x <- .as.EMPT(x,
                     experiment = experiment)
    .get.method.EMPT(x) <- method
  }else if(inherits(x,'EMPT')){
    x <- x
    .get.method.EMPT(x) <- method
    class(x) <- 'EMP_assay_data'
  }else {
    stop('Please check the input data')
  }
  if (use_cached == F) {
    memoise::forget(.EMP_decostand_m) %>% invisible()
  }
  EMPT <- .EMP_decostand_m(EMPT=x,bySample=bySample,logbase=logbase,pseudocount=pseudocount,...)
  .get.history.EMPT(EMPT) <- call
  class(EMPT) <- 'EMP_cluster_analysis'
  if (action == 'add') {
    return(EMPT)
  }else if(action == 'get') {
    return(.get.result.EMPT(EMPT))
  }else{
    warning('action should be one of add or get!')
  }
}
