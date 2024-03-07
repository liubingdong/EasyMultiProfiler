#' Title
#'
#' @param EMPT wait_for_add
#' @param method wait_for_add
#' @param ... wait_for_add
#' @importFrom tibble column_to_rownames
#' @importFrom vegan vegdist
#' @noRd
.EMP_beta_analysis <- function(EMPT,method,...) {
  assay_data  <- EMPT %>%
    .get.assay.EMPT() %>% tibble::column_to_rownames('primary')

  distance_result <- vegan::vegdist(assay_data,method=method,...) %>% as.matrix()

  EMPT@deposit[['distance_result']] <- distance_result
  EMPT@method <- method
  EMPT@algorithm <- 'beta_analysis'
  .get.info.EMPT(EMPT) <- 'EMP_beta_analysis'
  EMPT

}

EMP_beta_analysis_m <- memoise::memoise(.EMP_beta_analysis)


#' Title
#'
#' @param x wait_for_add
#' @param experiment wait_for_add
#' @param method wait_for_add
#' @param use_cached wait_for_add
#' @param action wait_for_add
#' @param ... wait_for_add
#' @importFrom memoise forget
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_beta_analysis <- function(x,experiment,method,use_cached = T,action='add',...) {
  call <- match.call()
  if (inherits(x,"MultiAssayExperiment")) {

    EMPT <- .as.EMPT(x,
                     experiment = experiment)

  }else if(inherits(x,'EMPT')) {
    EMPT <- x
  }else {
    stop('Please check the input data!')
  }

  if (use_cached == F) {
    memoise::forget(.EMP_alpha_analysis) %>% invisible()
  }

  EMPT <- EMP_beta_analysis_m(EMPT = EMPT,method = method,...)
  .get.history.EMPT(EMPT) <- call
  class(EMPT) <- 'EMP_beta_analysis'
  if (action == 'add') {
    return(EMPT)
  }else if(action == 'get') {
    return(.get.result.EMPT(EMPT))
  }else{
    warning('action should be one of add or get!')
  }
}
