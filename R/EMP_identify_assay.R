.EMP_assay_filter_default<- function(x,experiment,estimate_group=NULL,min = 0,
                                     min_ratio = 0.7,action='add') {
  call <- match.call()

  if (inherits(x,"MultiAssayExperiment")) {
    EMPT <- .as.EMPT(x,
                     experiment = experiment)
    .get.method.EMPT(EMPT) <- 'assay_filter'
    #.get.estimate_group.EMPT(EMPT) <- estimate_group
  }else if(inherits(x,'EMPT')) {
    EMPT <-x
    .get.method.EMPT(EMPT) <- 'assay_filter'
    #estimate_group <- .check_estimate_group.EMPT(x,estimate_group) ## estimate_group is necessary!
  }
  assay_name <- .get.assay_name.EMPT(EMPT)

  if (is.null(estimate_group)) {
    estimate_group <- .get.estimate_group.EMPT(EMPT)
    if(!is.null(estimate_group)){
      message('EMP_identify_assay will work according to estimate_group = ',estimate_group)
    }
  }else{
    if (estimate_group == 'none'){
      estimate_group <- NULL
    }else{
      estimate_group <- estimate_group
    }
  }

  if (!assay_name %in% c('counts','relative')) {
    stop("EMP_assay_filter only supprot counts and relative assay data!")
  }

  if (!(min >= 0 & min <= 1)) {
    stop('If EMP_assay_filter works based on default,min should be an number between 0 and 1!')
  }

  if (assay_name == 'relative') {
    data <- .get.assay.EMPT(EMPT)
  }else if(assay_name == 'counts') {
    data <- EMPT %>% EMP_decostand(method='relative',action='get')
  }


  if (min !=0) {
    #data %<>% dplyr::mutate_if(is.numeric, ~dplyr::if_else(. < !!min, 0, .))
    data[data < min] <- 0
  }

  if (!is.null(estimate_group)) {
    data <- .get.mapping.EMPT(EMPT,action='get') %>% dplyr::select(primary,!!estimate_group) %>%
      dplyr::right_join(data,by='primary')

    data %>%
      tidyr::pivot_longer(cols = c(-primary,-!!estimate_group),
                          names_to = 'feature',
                          values_to = 'abundance') %>%
      dplyr::group_by(!!dplyr::sym(estimate_group), feature) %>%
      dplyr::summarize(Prob = mean(abundance > 0),.groups='drop') %>%
      dplyr::filter(Prob >= !!min_ratio) %>% unique() %>% dplyr::pull(feature) -> id


  }else{
    message('No group or design set. Assuming all samples belong to one group.')

    data %>%
      tidyr::pivot_longer(cols = c(-primary),
                          names_to = 'feature',
                          values_to = 'abundance') %>%
      dplyr::group_by(feature) %>%
      dplyr::summarize(Prob = mean(abundance > 0),.groups='drop') %>%
      dplyr::filter(Prob >= !!min_ratio) %>% unique() %>% dplyr::pull(feature) -> id

  }

  EMPT %<>% EMP_filter(filterFeature = id,action = 'select')

  # delect EMP_filter history record,because history already contains EMP_assay_filter!
  temp_history <- .get.history.EMPT(EMPT)
  temp_history[[length(temp_history)]] <- NULL
  .get.history.EMPT(EMPT,replace = T) <- temp_history

  class(EMPT) <- 'EMP_assay_data'
  .get.info.EMPT(EMPT) <- 'EMP_assay_data'
  .get.estimate_group.EMPT(EMPT) <- estimate_group
  #.get.history.EMPT(EMPT) <- call

  if (action == 'add') {
    return(EMPT)
  }else if (action == 'get') {
    return(.get.assay.EMPT(EMPT))
  }

}

.EMP_assay_filter_bulk<- function(x,experiment,estimate_group=NULL,min = 10,min_ratio = 0.7,action='add') {
  call <- match.call()

  if (inherits(x,"MultiAssayExperiment")) {
    EMPT <- .as.EMPT(x,
                     experiment = experiment)
    .get.method.EMPT(EMPT) <- 'assay_filter'
    .get.estimate_group.EMPT(EMPT) <- estimate_group
  }else if(inherits(x,'EMPT')) {
    EMPT <-x
    .get.method.EMPT(EMPT) <- 'assay_filter'
    #estimate_group <- .check_estimate_group.EMPT(x,estimate_group) ## estimate_group is necessary!
  }
  assay_name <- .get.assay_name.EMPT(EMPT)

  if (is.null(estimate_group)) {
    estimate_group <- .get.estimate_group.EMPT(EMPT)
    if(!is.null(estimate_group)){
      message('EMP_identify_assay will work according to estimate_group = ',estimate_group)
    }
  }else{
    if (estimate_group == 'none'){
      estimate_group <- NULL
    }else{
      estimate_group <- estimate_group
    }
  }

  if (!assay_name %in% c('counts')) {
    stop("EMP_assay_filter based on edgeR only supprot counts assay data!")
  }

  if (min < 1) {
    stop('If EMP_assay_filter works based on edgeR,min should be an integer number!')
  }


  if (is.null(estimate_group)) {
    factor_of_interest <- NULL
  }else{
    factor_of_interest <- as.symbol(estimate_group)
  }

  EMPT%>%
    tidybulk::keep_abundant(factor_of_interest = !!factor_of_interest,
                            minimum_counts = min,
                            minimum_proportion = min_ratio)  %>% .get.row_info.EMPT() %>%
    dplyr::pull(feature) -> id

  EMPT %<>% EMP_filter(filterFeature = id,action = 'select')

  # delect EMP_filter history record,because history already contains EMP_assay_filter!
  temp_history <- .get.history.EMPT(EMPT)
  temp_history[[length(temp_history)]] <- NULL
  .get.history.EMPT(EMPT,replace = T) <- temp_history


  class(EMPT) <- 'EMP_assay_data'
  .get.info.EMPT(EMPT) <- 'EMP_assay_data'
  .get.estimate_group.EMPT(EMPT) <- estimate_group
  #.get.history.EMPT(EMPT) <- call

  if (action == 'add') {
    return(EMPT)
  }else if (action == 'get') {
    return(.get.assay.EMPT(EMPT))
  }
}

#' Title
#'
#' @param x wait_for_add
#' @param experiment wait_for_add
#' @param estimate_group wait_for_add
#' @param method wait_for_add
#' @param min wait_for_add
#' @param min_ratio wait_for_add
#' @param action wait_for_add
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_identify_assay <- function(x,experiment,estimate_group=NULL,
                               method=c('default','edgeR'),min=if (method == "edgeR") 10 else 0.001,
                               min_ratio = 0.7,action='add'){
  call <- match.call()
  switch(method,
         "default" = {
           EMPT <- x %>% .EMP_assay_filter_default(experiment=experiment,
                                                   estimate_group=estimate_group,
                                                   min=min,
                                                   min_ratio=min_ratio,
                                                   action='add')
         },
         "edgeR"={
           EMPT <- x %>% .EMP_assay_filter_bulk(experiment=experiment,
                                                estimate_group=estimate_group,
                                                min=min,
                                                min_ratio=min_ratio,
                                                action='add')
         },
         {
           print('method should be one of default or edgeR')
         }
  )
  if (action=='add') {
    .get.history.EMPT(EMPT) <- call
    return(EMPT)
  }else if(action=='get'){
    return(.get.assay.EMPT(EMPT))
  }else{
    warning('action should be one of add or get!')
  }
}
