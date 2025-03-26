#' @importFrom missRanger missRanger
.EMP_impute_EMPT <- function(obj,experiment,coldata = TRUE,assay = FALSE, rowdata = FALSE,.formula=. ~ .,
                       pmm.k = 10, num.trees = 1000, seed = 123,verbose = 0,action='add',...){
  #call <- match.call()
  if (is(obj,"MultiAssayExperiment")) {
    EMPT <- .as.EMPT(obj,
                     experiment = experiment)
  }else if(is(obj,'EMPT')) {
    EMPT <- obj
  }


  if (coldata == T) {
    if (any(is.na(.get.mapping.EMPT(EMPT)))) {
      .get.mapping.EMPT(EMPT) <- .get.mapping.EMPT(EMPT) %>%
        missRanger::missRanger(pmm.k = pmm.k, num.trees = num.trees, seed = seed,verbose = verbose,formula=.formula,...)
    }else{
      EMP_message('Coldata has no NA value!',color=31,order=1,show='warning')
    }
  }

  if (rowdata == T) {
    if (any(is.na(.get.row_info.EMPT(EMPT)))) {
      .get.row_info.EMPT(EMPT) <- .get.row_info.EMPT(EMPT) %>%
        missRanger::missRanger(pmm.k = pmm.k, num.trees = num.trees, seed = seed,verbose = verbose,formula=.formula,...)
    }else{
      EMP_message('Rowdata has no NA value!',color=31,order=1,show='warning')
    }
  }

  if (assay == T) {
    if (any(is.na(.get.assay.EMPT(EMPT)))) {
      .get.assay.EMPT(EMPT) <- .get.assay.EMPT(EMPT) %>%
        missRanger::missRanger(pmm.k = pmm.k, num.trees = num.trees, seed = seed,verbose = verbose,formula=.formula,...)
    }else{
      EMP_message('Assay data has no NA value!',color=31,order=1,show='warning')
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


#' @importFrom missRanger missRanger
#' @importFrom S4Vectors DataFrame
.EMP_impute_mae <- function(obj,.formula=. ~ .,
                       pmm.k = 10, num.trees = 1000, seed = 123,verbose = 0,action='add',...){
  coldata <- obj |> EMP_coldata_extract()
  if (any(is.na(coldata))) {
    new_coldata <- coldata %>%
      missRanger::missRanger(pmm.k = pmm.k, num.trees = num.trees, seed = seed,verbose = verbose,formula=.formula,...)
  }else{
    EMP_message('Coldata has no NA value!',color=31,order=1,show='warning')
  }

  if (action == 'add') {
    new_coldata <- new_coldata |>
      tibble::column_to_rownames('primary') |>
      DataFrame()
    obj@colData <- new_coldata
    return(obj)
  }else if(action == 'get'){
    return(new_coldata)
  }else {
    stop("action should be one of add or get")
  }
}




#' @importFrom memoise memoise
.EMP_impute_EMPT_m <- memoise::memoise(.EMP_impute_EMPT,cache = cachem::cache_mem(max_size = 4096 * 1024^2))
.EMP_impute_MAE_m <- memoise::memoise(.EMP_impute_mae,cache = cachem::cache_mem(max_size = 4096 * 1024^2))

#' Fast Imputation of Missing Values by Chained Random Forests
#'
#' @param obj Object in EMPT or MultiAssayExperiment format.
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object.
#' @param coldata A boolean. Whether the function use the coldata to impute or not.(default:TRUE)
#' @param assay A boolean. Whether the function use the assay to impute or not.(default:FALSE)
#' @param rowdata A boolean. Whether the function use the rowdata to impute or not.(default:FALSE)
#' @param .formula A formula.A two-sided formula specifying variables to be imputed (left hand side) and variables used to impute (right hand side). Defaults to . ~ ., i.e., use all variables to impute all variables. For instance, if all variables (with missings) should be imputed by all variables except variable "ID", use . ~ . - ID. Note that a "." is evaluated separately for each side of the formula. Further note that variables with missings must appear in the left hand side if they should be used on the right hand side.
#' @param pmm.k Number of candidate non-missing values to sample from in the predictive mean matching steps. 0 to avoid this step.
#' @param num.trees Tree number.
#' @param seed Integer seed to initialize the random generator.
#' @param verbose Controls how much info is printed to screen. 0 to print nothing.
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
#' @param action A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param ... Additional parameters, see also \code{\link[missRanger]{missRanger}}
#' @importFrom missRanger missRanger
#'
#' @return EMPT or MultiAssayExperiment object
#' @export
#'
#' @examples
#' data(MAE)
#' 
#' ## For MultiAssayExperiment
#' MAE |>
#'   EMP_impute(.formula = HAMA+HAMD ~ .)
#' 
#' ## For EMPT
#' ## For coldata
#' MAE |>
#'   EMP_assay_extract('geno_ec') |>
#'   EMP_impute(assay=FALSE,coldata=TRUE,rowdata=FALSE) |>
#'   EMP_coldata_extract()
#' 
#' ### Support formula, such as only impute SAS and SDS
#' MAE |>
#'   EMP_assay_extract('geno_ec') |>
#'   EMP_impute(.formula = HAMA+HAMD ~ .) |>
#'   EMP_coldata_extract() 
#' 
#' ## For assay
#' MAE |>
#'   EMP_assay_extract('geno_ec') |>
#'   EMP_impute(assay=TRUE,coldata=FALSE,rowdata=FALSE)
#' 
#' ## For rowdata (Not Not recommended)
#' MAE |>
#'   EMP_assay_extract('geno_ec') |>
#'   EMP_impute(assay=FALSE,coldata=FALSE,rowdata=TRUE)
EMP_impute <- function (obj,experiment=NULL,coldata = TRUE,assay = FALSE, rowdata = FALSE,.formula=. ~ .,
                       pmm.k = 10, num.trees = 1000, seed = 123,verbose = 0,use_cached=TRUE,action='add',...) {
  call <- match.call()
  deposit <- list()

  if (use_cached == FALSE) {
    memoise::forget(.EMP_impute_EMPT_m) %>% invisible()
    memoise::forget(.EMP_impute_MAE_m) %>% invisible()
  }
  
  if (is(obj,"MultiAssayExperiment") &  is.null(experiment)) {
    deposit <- .EMP_impute_MAE_m(obj=obj,.formula=.formula,
                       pmm.k = pmm.k, num.trees = num.trees, seed = seed,verbose = verbose,action=action,...)
  }else if(is(obj,"MultiAssayExperiment") &  !is.null(experiment)) {
    deposit <- .EMP_impute_EMPT_m(obj=obj,experiment=experiment,coldata = coldata,assay = assay, rowdata = rowdata,.formula=.formula,
                       pmm.k = pmm.k, num.trees = num.trees, seed = seed,verbose = verbose,action=action,...)
    if (action=='add') {
     .get.history.EMPT(deposit) <- call
    }
  }else if(is(obj,'EMPT')) {
    deposit <- .EMP_impute_EMPT_m(obj=obj,experiment=experiment,coldata = coldata,assay = assay, rowdata = rowdata,.formula=.formula,
                       pmm.k = pmm.k, num.trees = num.trees, seed = seed,verbose = verbose,action=action,...)
    if (action=='add') {
     .get.history.EMPT(deposit) <- call
    }    
  }
  return(deposit)
}


