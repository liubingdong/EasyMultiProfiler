#' @importFrom vegan decostand
#' @noRd
.EMP_decostand <- function(EMPT,method,bySample='default',logbase,pseudocount=pseudocount,...){
  primary <- NULL

  if (bySample != 'default' & !is.logical(bySample)) {
   stop('Paramer bySample should be default, TRUE or FALSE!')
  }

  MARGIN <- check_MARGIN(method=method,bySample=bySample)
  method_name <- check_method_name(method=method,suffix=logbase)
  message_info <- message_MARGIN(method=method,MARGIN = MARGIN)

  ## .calc_rclr unused argument pseudocount
  if (method == 'rclr') {
    assay_decostand_data <- assay(EMPT) %>% t() %>% 
    vegan::decostand(method = method, MARGIN = MARGIN,logbase=logbase,...)
  }else if(method == 'integer'){
    assay_decostand_data <- assay(EMPT) %>% t() %>% 
     round(digits = 0)
  }else{
    method_trans <- parse_log_string(method,suffix=logbase)

    if (unlist(method_trans)[1] == 'relative') {
      method_trans[[1]] <- 'total'
    }  

    if (unlist(method_trans)[1] == 'log') {
       assay_decostand_data <- assay(EMPT) %>% t()
       if (!is.na(method_trans$log_addition)) {
         assay_decostand_data <- assay_decostand_data + method_trans$log_addition
       }    
       assay_decostand_data[assay_decostand_data > 0 & !is.na(assay_decostand_data)] <- log(assay_decostand_data[assay_decostand_data > 0 & !is.na(assay_decostand_data)], 
        base = method_trans$suffix) 
    }else{
      assay_decostand_data <- assay(EMPT) %>% t() %>% 
        vegan::decostand(method = unlist(method_trans)[1], MARGIN = MARGIN,pseudocount=pseudocount,...)
    }
  }

  .get.assay.EMPT(EMPT) <- assay_decostand_data %>% as.data.frame() %>% tibble::rownames_to_column('primary') %>% tibble::as_tibble() %>% dplyr::select(primary,everything())
  .get.method.EMPT(EMPT) <- method_name
  .get.assay_name.EMPT(EMPT) <- method_name
  .get.message_info.EMPT(EMPT) <- message_info
  .get.algorithm.EMPT(EMPT) <- 'EMP_decostand'
  .get.info.EMPT(EMPT) <- 'EMP_decostand'
  return(EMPT)

}

.EMP_decostand_m <- memoise::memoise(.EMP_decostand,cache = cachem::cache_mem(max_size = 4096 * 1024^2))

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
#'   hellinger, log, alr, clr, rclr, integer
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
#' 11. log: logarithmic transformation(Here is differrent from decostand in the vegan package). The transformation only applies a logarithmic change to values greater than 0, while values equal to 0 remain unchanged. 
#' The method can be written as log2+1, which means adding 1 to all data first, followed by a base-2 logarithmic transformation. 
#' If written as log or log+1, it will operate according to the specified logarithmic base.
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
#' 15. integer: rounding of numbers using the function Round.
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
  if (is(obj,"MultiAssayExperiment")) {
    x <- .as.EMPT(obj,
                     experiment = experiment)
  }else if(is(obj,'EMPT')){
    x <- obj
    class(x) <- 'EMP_assay_data'
  }else {
    stop('Please check the input data for EMP_decostand!')
  }
  if (use_cached == FALSE) {
    memoise::forget(.EMP_decostand_m) %>% invisible()
  }
  #method <- match.arg(method,c("relative","total","max","frequency","normalize","range",
  #  "rank","standardize","pa","hellinger","log","alr","clr","rclr","integer"))
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


#' @importFrom stringr str_detect
#' @importFrom stringr str_match
parse_log_string <- function(input_string,suffix=2) {
  # 检查是否包含 log
  if (!str_detect(input_string, "^log")) {
    return(input_string)
  }
  
  # 匹配不同的情况
  match <- str_match(input_string, "^log(\\d*)(?:\\+(\\d+))?$")
  
  if (is.na(match[1, 1])) {
    stop("Please check the method format in the EMP_decostand.")
  }
  
  # 提取 log 和相关部分
  log_part <- "log"
  suffix <- ifelse(match[1, 2] == "", suffix, as.numeric(match[1, 2]))
  addition <- ifelse(is.na(match[1, 3]), NA, as.numeric(match[1, 3]))
  
  # 返回结果
  return(list(
    log = log_part,
    suffix = suffix,
    log_addition = addition
  ))
}

log_name_create<- function(parse_log_string){
  if (!is.na(parse_log_string$log_addition)) {
    log_method_name <- paste0(parse_log_string$log, parse_log_string$suffix,'+',parse_log_string$log_addition)
  }else{
    log_method_name <- paste0(parse_log_string$log, parse_log_string$suffix)
  }
  return(log_method_name)
}

check_MARGIN <- function (method,bySample) {
  method <- parse_log_string(method) |> unlist()
  if (method[1] == 'relative') {
    method[1] <- 'total'
  }
  if (bySample != 'default' & is.logical(bySample)) {
    MARGIN <- ifelse(bySample, 1, 2)
  }else {
    switch(method[1],
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
           }, 
           "integer" = {
             MARGIN <- 1
           },              
           {
             stop("MARGIN is missing!")
           }
    )
  }
  return(MARGIN)
}

message_MARGIN <- function(method,MARGIN=MARGIN,suffix=2){
  message_info <- list()
  method_trans <- parse_log_string(method,suffix=suffix)
  method_trans_unlist <- unlist(method_trans)
  if (!method_trans_unlist[1] %in% c('log','pa','integer')) {
    decostand_info <- c('Sample','Feature')
    message_info %<>% append(paste0('Standardization method proceed according to ',decostand_info[MARGIN],'!')) 
  }
  
  if (method_trans_unlist[1] == 'log') {
    method_name <- log_name_create(method_trans)
    message_info %<>% append(paste0('Standardization method proceed by ',method_name,'!')) 
  }
  
  return(message_info)
}

check_method_name <- function(method,suffix=2){
  method_trans <- parse_log_string(method,suffix=suffix)
  method_trans_unlist <- unlist(method_trans)
  
  if (method_trans_unlist[1]=="log"){
    method_name <- log_name_create(parse_log_string = method_trans)
  }else{
    method_name <- method_trans
  }
  
  if (method_trans_unlist[1] == 'total') {
    method_name  <- 'relative'
  }
  
  return(method_name)
}