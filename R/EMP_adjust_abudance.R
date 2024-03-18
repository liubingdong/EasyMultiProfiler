#' Title
#'
#' @param EMPT object in EMPT format.
#' @param method A character string. Methods include combat_seq (default), combat and limma_remove_batch_effect.
#' @param .factor_unwanted A tidy select, e.g. column names without double quotation. c(batch, country) These are the factor that we want to adjust for, including unwanted batcheffect, and unwanted biological effects.
#' @param .factor_of_interest A tidy select, e.g. column names without double quotation. c(treatment) These are the factor that we want to preserve.
#' @param ... Further parameters passed to the function tidybulk::adjust_abundance
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
.EMP_adjust_abudance <- function(EMPT,method,.factor_unwanted,.factor_of_interest,...) {

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
.EMP_adjust_abudance_m <- memoise::memoise(.EMP_adjust_abudance)


#' Title
#'
#' @param x Object in EMPT format.
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object.
#' @param method A character string. Methods include combat_seq (default), combat and limma_remove_batch_effect.
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
#' @param .factor_unwanted A tidy select, e.g. column names without double quotation. c(batch, country) These are the factor that we want to adjust for, including unwanted batcheffect, and unwanted biological effects.
#' @param .factor_of_interest A tidy select, e.g. column names without double quotation. c(treatment) These are the factor that we want to preserve.
#' @param action A character string.A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param ... Further parameters passed to the function tidybulk::adjust_abundance
#' @importFrom memoise forget
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_adjust_abudance <- function(x,experiment,
                                method='combat_seq',use_cached=T,
                                .factor_unwanted,.factor_of_interest,action='add',...) {
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
    memoise::forget(.EMP_adjust_abudance_m) %>% invisible()
  }
  EMPT <- .EMP_adjust_abudance_m(EMPT=x,method=method,.factor_unwanted,.factor_of_interest,...)
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
