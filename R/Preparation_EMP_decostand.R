#' @importFrom vegan decostand
#' @noRd
.EMP_decostand <- function(EMPT,bySample=TRUE,logbase,pseudocount=pseudocount,...){
  primary <- NULL
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
    vegan::decostand(method = method, MARGIN = MARGIN,logbase=logbase,...)
  }else if(method == 'integer'){
    assay_decostand_data <- assay_data %>%  tibble::column_to_rownames('primary') %>% t() %>%
     round(digits = 0)
  }else {
    assay_decostand_data <- assay_data %>%  tibble::column_to_rownames('primary') %>% t() %>%
    vegan::decostand(method = method, MARGIN = MARGIN,logbase=logbase,pseudocount=pseudocount,...)
  }


  .get.assay.EMPT(EMPT) <- assay_decostand_data %>% t() %>% as.data.frame() %>% tibble::rownames_to_column('primary') %>% tibble::as_tibble() %>% dplyr::select(primary,everything())
  .get.method.EMPT(EMPT) <- method_name
  .get.assay_name.EMPT(EMPT) <- method_name
  .get.algorithm.EMPT(EMPT) <- 'EMP_decostand'
  .get.info.EMPT(EMPT) <- 'EMP_decostand'
  EMPT

}

.EMP_decostand_m <- memoise::memoise(.EMP_decostand)

#' Standardization Methods
#'
#' @param x Object in EMPT or MultiAssayExperiment format.
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object.
#' @param method A character string. Standardization method inluding relative(Equivalent to total in the vegan::decostand), integer, log, clr, alr, aclr,and more. Details see vegan::decostand.
#' @param bySample A boolean. Whether the function decostand by the sample or feature.
#' @param logbase An interger. The logarithm base used in method = "log".(default=2)
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
#' @param pseudocount A number. The logarithm pseudocount used in method = "clr".(default=0.0000001)
#' @param action A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param ... Further parameters passed to the function vegan::decostand.
#'
#' @return EMPT object
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