#' @importFrom tidybulk adjust_abundance
#' @importFrom tibble rownames_to_column
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom stats aov
#' @importFrom stats as.dendrogram
#' @importFrom stats as.dist
#' @importFrom stats ave
#' @importFrom stats cor
#' @importFrom stats median
#' @importFrom stats na.omit
#' @importFrom stats predict
#' @importFrom stats pt
#' @importFrom stats sd
#' @importFrom stats setNames
#' @importFrom utils combn
#' @importFrom utils read.table
#' @importFrom SummarizedExperiment assay
#' @noRd
.EMP_adjust_abundance <- function(EMPT,method,.factor_unwanted,.factor_of_interest,...) {

  primary <- NULL
  
  coldata <- .get.mapping.EMPT(EMPT) %>% dplyr::select(primary,!!.factor_unwanted,!!.factor_of_interest)
  
  if(any(is.na(coldata))) {
    stop('Column ',.factor_unwanted,' or ',.factor_of_interest,' has beed deteced missing value, please check and filter them!')
  }  

  EMPT %>%
      tidybulk::adjust_abundance( .factor_unwanted = !!.factor_unwanted, .factor_of_interest =  !!.factor_of_interest,
                                  method=method,...) %>% suppressWarnings() %>% suppressMessages() -> EMPT_adjust


  new_assay <- assay(EMPT_adjust,i = 'counts_adjusted') %>% t() %>% as.data.frame %>%
    tibble::rownames_to_column('primary') %>%
    tibble::as_tibble()

  .get.assay.EMPT(EMPT) <- new_assay

  return(EMPT)
}

#' @importFrom memoise memoise
.EMP_adjust_abundance_m <- memoise::memoise(.EMP_adjust_abundance,cache = cachem::cache_mem(max_size = 4096 * 1024^2))


#' Adjust experssion or abundance for unexpected bias or batch effect
#'
#' @param obj Object in EMPT or SummarisedExperiment format.
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object.
#' @param method A character string. Methods include combat_seq (default), combat and limma_remove_batch_effect.
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
#' @param .factor_unwanted A tidy select, e.g. column names without double quotation. c(batch, country) These are the factor that we want to adjust for, including unwanted batcheffect, and unwanted biological effects.
#' @param .factor_of_interest A tidy select, e.g. column names without double quotation. c(treatment) These are the factor that we want to preserve.
#' @param action A character string.A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param ... Further parameters passed to the function \code{\link[tidybulk]{adjust_abundance}}.
#' @importFrom memoise forget
#'
#' @return EMPT object
#' @export
#'
#' @examples
#' data(MAE)
## combat_seq method 
#' MAE |>
#'   EMP_assay_extract(experiment='geno_ko') |>
#'   EMP_adjust_abundance(.factor_unwanted = 'Region',.factor_of_interest = 'Group',
#'                       method = 'combat_seq',action = 'add') 
#' ## combat method 
#' MAE |>
#'   EMP_assay_extract(experiment='geno_ko') |>
#'   EMP_adjust_abundance(.factor_unwanted = 'Region',.factor_of_interest = 'Group',
#'                       method = 'combat',action = 'add') 
#' ## limma_remove_batch_effect
#' MAE |>
#'   EMP_assay_extract(experiment='geno_ko') |>
#'   EMP_adjust_abundance(.factor_unwanted = 'Region',.factor_of_interest = 'Group',
#'                       method = 'limma_remove_batch_effect') 
EMP_adjust_abundance <- function(obj,experiment,
                                method='combat_seq',use_cached=TRUE,
                                .factor_unwanted,.factor_of_interest,action='add',...) {
  call <- match.call()
  if (is(obj,"MultiAssayExperiment")) {
    x <- .as.EMPT(obj,
                  experiment = experiment)
    .get.method.EMPT(x) <- method
  }else if(is(obj,'EMPT')){
    x <- obj
    .get.method.EMPT(x) <- method
    class(x) <- 'EMP_assay_data'
  }else {
    stop('Please check the input data for EMP_adjust_abundance!')
  }
  if (use_cached == FALSE) {
    memoise::forget(.EMP_adjust_abundance_m) %>% invisible()
  }

  # avoid the factor error after filter
  try(x[[.factor_unwanted]] <- droplevels(x[[.factor_unwanted]]),silent = TRUE)
  try(x[[.factor_of_interest]] <- droplevels(x[[.factor_of_interest]]),silent = TRUE)
  
  EMPT <- .EMP_adjust_abundance_m(EMPT=x,method=method,.factor_unwanted,.factor_of_interest,...)
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
