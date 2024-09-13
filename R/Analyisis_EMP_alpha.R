#' @param EMPT wait_for_add
#' @param ... wait_for_add
#' @importFrom tibble column_to_rownames
#' @importFrom vegan diversity
#' @importFrom dplyr select
#' @importFrom dplyr everything
#' @noRd
.EMP_alpha_analysis <- function(EMPT,...) {
  primary <- NULL
  assay_data  <- assay(EMPT) %>% t()
  alpha_shannon <- vegan::diversity(assay_data,,index = 'shannon',...)
  alpha_simpson <- vegan::diversity(assay_data,index = 'simpson',...)
  alpha_invsimpson <- vegan::diversity(assay_data,index = 'invsimpson',...)
  diversity_result <- data.frame(shannon=alpha_shannon,
                       simpson=alpha_simpson,
                       invsimpson=alpha_invsimpson) %>%
    tibble::rownames_to_column('primary') %>%
    dplyr::select(primary,everything()) %>% tibble::as_tibble()
  if (.get.assay_name.EMPT(EMPT) == 'counts') {
    assay_data_abs <- assay_data %>% round(digits = 0) %>% vegan::estimateR()
    observerd_index <- assay_data_abs[1,]
    chao1 <- assay_data_abs[2,]
    ACE <- assay_data_abs[4,]

    diversity_result <- data.frame(shannon=alpha_shannon,
                         simpson=alpha_simpson,
                         invsimpson=alpha_invsimpson,
                         observerd_index,
                         chao1,
                         ACE) %>%
      tibble::rownames_to_column('primary') %>%
      dplyr::select(primary,everything()) %>% tibble::as_tibble()
  }
  EMPT@deposit[['diversity_result']] <- diversity_result
  .get.algorithm.EMPT(EMPT) <- 'alpha_analysis'
  .get.info.EMPT(EMPT) <- 'EMP_alpha_analysis'
  EMPT
}

#' @importFrom memoise memoise
EMP_alpha_analysis_m <- memoise::memoise(.EMP_alpha_analysis,cache = cachem::cache_mem(max_size = 2048 * 1024^2))

#' Diversity Indices
#'
#' @param obj Object in EMPT format.
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object.
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
#' @param action A character string.A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param ... Further parameters passed to the function vegan::diversity
#' @importFrom memoise forget
#'
#' @return EMPT object
#' @export
#'
#' @examples
#' data(MAE)
#' MAE |>
#'   EMP_assay_extract(experiment='taxonomy')|> 
#'   EMP_alpha_analysis()
#' 
#' MAE |>
#'   EMP_assay_extract(experiment='geno_ec') |> 
#'   EMP_alpha_analysis()
EMP_alpha_analysis <- function(obj,experiment,use_cached = TRUE,action='add',...) {
  call <- match.call()
  if (is(obj,"MultiAssayExperiment")) {

    EMPT <- .as.EMPT(obj,
                     experiment = experiment)

  }else if(is(obj,'EMPT')) {
    EMPT <- obj
  }else {
    stop('Please check the input data for EMP_alpha_analysis!')
  }


  if (use_cached == FALSE) {
    memoise::forget(.EMP_alpha_analysis) %>% invisible()
  }

  ## only integer, counts and relative data could used in alpha_analysis
  if(!.get.assay_name.EMPT(EMPT) %in% c('counts','relative','integer')) {
    stop('only integer, counts and relative data could used in alpha_analysis')
  }

  EMPT <- EMP_alpha_analysis_m(EMPT = EMPT,...)
  .get.history.EMPT(EMPT) <- call
  class(EMPT) <- 'EMP_alpha_analysis'
  if (action == 'add') {
    return(EMPT)
  }else if(action == 'get') {
    return(.get.result.EMPT(EMPT))
  }else{
    warning('action should be one of add or get!')
  }
}
