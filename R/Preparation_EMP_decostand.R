#' @importFrom vegan decostand
#' @noRd
.EMP_decostand <- function(EMPT,method,bySample='default',logbase,pseudocount=pseudocount,...){
  primary <- NULL

  if (method == 'relative') {
    method  <- 'total'
  }

  .get.assay.EMPT(EMPT) -> assay_data

  if (bySample != 'default' & !is.logical(bySample)) {
   stop('Paramer bySample should be default, TRUE or FALSE!')
  }

  if (bySample != 'default' & is.logical(bySample)) {
    MARGIN <- ifelse(bySample, 1, 2)
  }else {
    switch(method,
           "total" = {
            MARGIN <- 1
           },
           "max" = {
             MARGIN <- 2
           },
           "frequency" = {
             MARGIN <- 2
           },
           "normalize" = {
             MARGIN <- 1
           },
           "range" = {
             MARGIN <- 2
           },
           "rank" = {
             MARGIN <- 1
           },
           "standardize" = {
             MARGIN <- 2
           },
           "pa" = {
             MARGIN <- 1
           },
           "chi.square" = {
             MARGIN <- 1
           },
           "hellinger" = {
             MARGIN <- 1
           },
           "log" = {
             MARGIN <- 1
           },
           "alr" = {
             MARGIN <- 1
           },
           "clr" = {
             MARGIN <- 1
           },
           "rclr" = {
             MARGIN <- 1
           }
    )
  }

  # Only the log, pa and integer algorithm do not need to consider the direction of standardization.
  if (!method %in% c('log','pa','integer')) {
    decostand_info <- c('Sample','Feature')
    message_info <- list()
    message_info %<>% append(paste0('Standardization method proceed by ',decostand_info[MARGIN],'!')) 
    .get.message_info.EMPT(EMPT) <- message_info
  }


  if (method=="log"){
    method_name <- paste0(method, logbase)
  }else{
    method_name <- method
  }

  ## .calc_rclr unused argument pseudocount
  if (method == 'rclr') {
    assay_decostand_data <- assay_data %>%  tibble::column_to_rownames('primary') %>% 
    vegan::decostand(method = method, MARGIN = MARGIN,logbase=logbase,...)
  }else if(method == 'integer'){
    assay_decostand_data <- assay_data %>%  tibble::column_to_rownames('primary') %>% 
     round(digits = 0)
  }else {
    assay_decostand_data <- assay_data %>%  tibble::column_to_rownames('primary') %>% 
    vegan::decostand(method = method, MARGIN = MARGIN,logbase=logbase,pseudocount=pseudocount,...)
  }


  .get.assay.EMPT(EMPT) <- assay_decostand_data %>% tibble::rownames_to_column('primary') %>% tibble::as_tibble() %>% dplyr::select(primary,everything())
  .get.method.EMPT(EMPT) <- method_name
  .get.assay_name.EMPT(EMPT) <- method_name
  .get.algorithm.EMPT(EMPT) <- 'EMP_decostand'
  .get.info.EMPT(EMPT) <- 'EMP_decostand'
  EMPT

}

.EMP_decostand_m <- memoise::memoise(.EMP_decostand)

#' Standardization Methods
#'
#' @param obj Object in EMPT or MultiAssayExperiment format.
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object.
#' @param method A character string. Details see vegan::decostand.
#' @param bySample A boolean. Whether the function decostand by the sample or feature.
#' @param logbase An interger. The logarithm base used in method = "log".(default=2)
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
#' @param pseudocount A number. The logarithm pseudocount used in method = "clr" or "alr".(default=0.0000001)
#' @param action A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param ... Further parameters passed to the function vegan::decostand.
#' @section Detaild about method:
#' The following method from vegan::decostand are available：
#' \describe{
#'   relative(total), max, frequency, normalize, range,
#'   rank, standardize, pa, chi.square,
#'   hellinger, log, alr, clr, rclr
#' }
#'
#' When standardizing data, it's important to consider whether the data processing direction should be by sample or by feature. 
#' All methods have a default parameter bySample in EMP_decostand function. Users could change this parameter according to study design. 
#' For more detailed information, please refer to the decostand function in the vegan package.
#'
#' Here is a brife introduction of document from vegan package:
#'
#' 1. relative(total): divide by margin total, also called relative abundance. (default bySample = TRUE).
#'
#' 2. max: divide by margin maximum (default bySample = FALSE).
#'
#' 3. frequency: divide by margin total and multiply by the number of non-zero items, so that the average of non-zero entries is one (Oksanen 1983; default bySample = FALSE).
#'
#' 4. normalize: make margin sum of squares equal to one (default bySample = TRUE).
#'
#' 5. range: standardize values into range 0 ... 1 (default bySample = FALSE). If all values are constant, they will be transformed to 0.
#'
#' 6. rank: rank replaces abundance values by their increasing ranks leaving zeros unchanged, and rrank is similar but uses relative ranks with maximum 1 (default bySample = TRUE). Average ranks are used for tied values.
#'
#' 7. standardize: scale x to zero mean and unit variance (default bySample = FALSE).
#'
#' 8. pa: scale x to presence/absence scale (0/1).
#'
#' 9. chi.square: divide by row sums and square root of column sums, and adjust for square root of matrix total (Legendre & Gallagher 2001). 
#' When used with the Euclidean distance, the distances should be similar to the Chi-square distance used in correspondence analysis. 
#' However, the results from cmdscale would still differ, since CA is a weighted ordination method (default bySample = TRUE).
#'
#' 10. hellinger: square root of method = "relative" (Legendre & Gallagher 2001).(default bySample = TRUE)
#'
#' 11. log: logarithmic transformation. Higher bases give less weight to quantities and more to presences, and logbase = Inf gives the presence/absence scaling.
#'
#' 12. alr: Additive log ratio ("alr") transformation (Aitchison 1986) reduces data skewness and compositionality bias. 
#' The transformation assumes positive values, pseudocounts can be added with the argument pseudocount. 
#' One of the rows/columns is a reference that can be given by reference (name of index). The first row/column is used by default (reference = 1). 
#' Note that this transformation drops one row or column from the transformed output data. 
#' This transformation is often used with pH and other chemistry measurenments. 
#' It is also commonly used as multinomial logistic regression. Default bySample = TRUE uses row as the reference.
#'
#' 13. clr: centered log ratio ("clr") transformation proposed by Aitchison (1986) reduces data skewness and compositionality bias. 
#' This transformation has frequent applications in microbial ecology (see e.g. Gloor et al., 2017).
#' The method can operate only with positive data; a common way to deal with zeroes is to add pseudocount,  
#  either by adding it manually to the input data, or by using the argument pseudocount as in decostand(x, method = "clr", pseudocount = 1).
#' Adding pseudocount will inevitably introduce some bias; see the rclr method for one available solution. (default bySample = TRUE)
#'
#' 14. rclr: square root of method = "relative" (Legendre & Gallagher 2001). robust clr ("rclr") is similar to regular clr (see above) but allows data that contains zeroes.
#' This method does not use pseudocounts, unlike the standard clr.
#' Robust clr divides the values by geometric mean of the observed features; zero values are kept as zeroes, and not taken into account. 
#' In high dimensional data, the geometric mean of rclr is a good approximation of the true geometric mean. (default bySample = TRUE)
#'  
#' @return EMPT object
#' @export
#'
#' @examples
#' data(MAE)
#' ## Transfer data into relative format.
#' MAE |>
#'   EMP_decostand(experiment = 'taxonomy',method = 'relative') 
#' 
#' ## Transfer data into centered log ratio ("clr") format.
#' MAE |>
#'   EMP_decostand(experiment = 'geno_ko',method = 'clr',pseudocount=0.0001)
#' 
#' ## Transfer data into logformat.
#' MAE |>
#'   EMP_decostand(experiment = 'geno_ec',method = 'log',logbase = 2) 
EMP_decostand <- function(obj,experiment,method,bySample='default',logbase =2,use_cached = TRUE,pseudocount=0.0000001,action='add',...){
  call <- match.call()
  if (inherits(obj,"MultiAssayExperiment")) {
    x <- .as.EMPT(obj,
                     experiment = experiment)
  }else if(inherits(obj,'EMPT')){
    x <- obj
    class(x) <- 'EMP_assay_data'
  }else {
    stop('Please check the input data')
  }
  if (use_cached == F) {
    memoise::forget(.EMP_decostand_m) %>% invisible()
  }
  EMPT <- .EMP_decostand_m(EMPT=x,method=method,bySample=bySample,logbase=logbase,pseudocount=pseudocount,...)
  .get.history.EMPT(EMPT) <- call
  class(EMPT) <- 'EMP_decostand'
  if (action == 'add') {
    return(EMPT)
  }else if(action == 'get') {
    return(.get.result.EMPT(EMPT))
  }else{
    warning('action should be one of add or get!')
  }
}
