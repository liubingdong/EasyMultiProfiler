.EMP_rrarefy <- function(EMPT, raresize=NULL,seed=123,only_show_depth = F, ...){
  primary <- feature <- value <- mean_value <- sd_value <- min_sd <- real_sd <- NULL
  real_mean <- mean_diff <- sd_sum <- NULL
  trim_sample_id <- NULL
  trim_feature_id <- NULL
  message_info <- list()
  assay_data <- .get.result.EMPT(EMPT) %>% tibble::column_to_rownames('primary') %>% round(digits = 0)
  total_sample_num <- nrow(assay_data)
  total_feature_num <- ncol(assay_data )
  min_raresize <- min(rowSums(assay_data))
  max_raresize <- max(rowSums(assay_data))
  if (only_show_depth ==T) {
      message_wrap('Min depth: ',min_raresize)
      message_wrap('Max depth: ',max_raresize)

  }else {
      message_info %<>% append(paste0('Min depth: ',min_raresize))
      message_info %<>% append(paste0('Max depth: ',max_raresize))

    if (is.null(raresize)) {
      raresize <- min(rowSums(assay_data))
      message_info %<>% append(paste0('Current raresize: ',raresize))
    }else if (is.numeric(raresize)) {
      if (raresize > max_raresize) {
        stop(paste0('raresize value should below '),max_raresize)
      }else {
        message_info %<>% append(paste0('Current raresize: ',raresize))
        trim_sample_id <- rownames(assay_data)[rowSums(assay_data) < raresize]
        assay_data <- assay_data[rowSums(assay_data) >= raresize, ,drop=FALSE]
        if (length(trim_sample_id) != 0) {
          message_info %<>% append(paste0(length(trim_sample_id)," of ",total_sample_num," Samples were removed, ",
                       "because their abundance were below the rarefaction: ",
                       paste(trim_sample_id,collapse = ' ')))
        }
      }
    }else
    { stop('parameter raresize should be numeric!')}


    res <- withr::with_seed(seed, vegan::rrarefy(x=assay_data, sample=raresize,...)) %>% suppressWarnings()
    remove_feature <- colSums(res)==0
    trim_feature_id <- names(remove_feature[remove_feature])

    if (length(trim_feature_id) != 0) {
      message_info %<>% append(paste0(length(trim_feature_id)," of ",total_feature_num, " Features were removed, ",
                   "because they are no longer present in any sample after rarefaction: ",
                   paste(trim_feature_id,collapse = ' ')))
    }

    rrarefy_data <- res[, !remove_feature,drop=F] %>%
      as.data.frame() %>%
      tibble::rownames_to_column('primary') %>%
      dplyr::select(primary,everything()) %>% tibble::as_tibble()

    .get.assay.EMPT(EMPT) <- rrarefy_data
    .get.info.EMPT(EMPT) <- 'EMP_rrarefy'
    .get.method.EMPT(EMPT) <- paste0('rrarefy:',raresize)
    .get.message_info.EMPT(EMPT) <- message_info
    return(EMPT)
  }
}



.EMP_rrarefy_m <- memoise::memoise(.EMP_rrarefy)

#' Rarefaction abundance or experssion Richness
#'
#' @param x Object in EMPT or MultiAssayExperiment format.
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object.
#' @param raresize An interger. Subsample size for rarefying community.
#' @param seed An interger. Set the random seed to the plot.(default:123)
#' @param only_show_depth A boolean. Whether the function only show the depth or not.
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
#' @param action A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param ... Further parameters passed to the function vegan::rrarefy.
#'
#' @return EMPT object
#' @export
#'
#' @examples
#' data(MAE)
#' ## rarefythe data according to the lowest abundance.
#' MAE |>
#'   EMP_rrarefy(experiment = 'taxonomy')
#' ## 
#' MAE |>
#'   EMP_rrarefy(experiment = 'taxonomy',only_show_depth=TRUE) # only show the depth
#' 
#' MAE |>
#'   EMP_rrarefy(experiment = 'taxonomy',raresize=1000) # Set a specific threshold

EMP_rrarefy <- function(x,experiment,raresize=NULL,seed=123,only_show_depth=FALSE,use_cached = TRUE,action = 'add',...) {
  call <- match.call()
  if (inherits(x,"MultiAssayExperiment")) {
    x <- .as.EMPT(x,
                  experiment = experiment)
  }else if(inherits(x,'EMPT')){
    x <- x
    class(x) <- 'EMP_assay_data'
  }else {
    stop('Please check the input data')
  }
  if (use_cached == F) {
    memoise::forget(.EMP_rrarefy_m) %>% invisible()
  }
  
  if (.get.assay_name.EMPT(x) != 'counts') {
     warning("If assay data was not raw counts, EMP_rrarefy may not work correctly! ")
  }

  EMPT <- .EMP_rrarefy_m(EMPT=x,raresize=raresize,seed=seed,only_show_depth=only_show_depth,...)
  # in case that .EMP_rrarefy only return depth message
  if (inherits(EMPT,'EMPT')) {
    .get.history.EMPT(EMPT) <- call
    class(EMPT) <- 'EMP_assay_data'
  }
  if (action == 'add') {
    return(EMPT)
  }else if(action == 'get') {
    return(.get.result.EMPT(EMPT))
  }else{
    warning('action should be one of add or get!')
  }

}
