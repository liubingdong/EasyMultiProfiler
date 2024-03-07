#' Title
#'
#' @param x wait_for_add
#' @param experiment wait_for_add
#' @param coldata wait_for_add
#' @param assay wait_for_add
#' @param rowdata wait_for_add
#' @param pmm.k wait_for_add
#' @param num.trees wait_for_add
#' @param seed wait_for_add
#' @param verbose wait_for_add
#' @param action wait_for_add
#' @param ... wait_for_add
#' @importFrom missRanger missRanger
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_impute <- function(x,experiment,coldata = T,assay = F, rowdata = F,
                       pmm.k = 10, num.trees = 1000, seed = 123,verbose = 0,action='add',...){
  call <- match.call()
  if (inherits(x,"MultiAssayExperiment")) {
    EMPT <- .as.EMPT(x,
                     experiment = experiment)
  }else if(inherits(x,'EMPT')) {
    EMPT <-x
  }


  if (coldata == T) {
    if (any(is.na(.get.mapping.EMPT(EMPT)))) {
      .get.mapping.EMPT(EMPT) <- .get.mapping.EMPT(EMPT) %>%
        missRanger::missRanger(pmm.k = pmm.k, num.trees = num.trees, seed = seed,verbose = verbose,...)
    }else{
      message('Coldata has not NA value!')
    }
  }

  if (rowdata == T) {
    if (any(is.na(.get.row_info.EMPT(EMPT)))) {
      .get.row_info.EMPT(EMPT) <- .get.row_info.EMPT(EMPT) %>%
        missRanger::missRanger(pmm.k = pmm.k, num.trees = num.trees, seed = seed,verbose = verbose,...)
    }else{
      message('Rowdata has not NA value!')
    }
  }

  if (assay == T) {
    if (any(is.na(.get.assay.EMPT(EMPT)))) {
      .get.assay.EMPT(EMPT) <- .get.assay.EMPT(EMPT) %>%
        missRanger::missRanger(pmm.k = pmm.k, num.trees = num.trees, seed = seed,verbose = verbose,...)
    }else{
      message('Assay data has not NA value!')
    }
  }

  if (action == 'add') {
    return(EMPT)
  }else if(action == 'get'){
    tbl_deposit <- list()
    if (coldata == T) {
      tbl_deposit[['coldata']] <- .get.mapping.EMPT(EMPT)
    }
    if (rowdata == T) {
      tbl_deposit[['row_info']] <- .get.row_info.EMPT(EMPT)
    }
    if (assay == T) {
      tbl_deposit[['assay']] <- .get.assay.EMPT(EMPT)
    }
    return(tbl_deposit)
  }else {
    stop("action should be one of add or get")
  }
}
