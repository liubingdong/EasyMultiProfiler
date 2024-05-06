.EMP_assay_filter_default<- function(x,experiment,estimate_group=NULL,min = 0,
                                     min_ratio = 0.7,action='add') {
  primary <- feature <- abundance <- Prob <- NULL
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

    ## In case that privious estimate_group is not in the coldata,eg EMP_collapse, make the estimate_group = NULL
    if (!is.null(estimate_group)) {
      col_name <- .get.mapping.EMPT(EMPT) %>% colnames()
      if (!estimate_group %in% col_name) {
        estimate_group <- NULL
      }else {
        message('EMP_identify_assay will work according to estimate_group = ',estimate_group)
      }
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
    data <- EMPT %>% EMP_decostand(method='relative',action='get') %>% suppressMessages()
  }


  if (min !=0) {
    #data %<>% dplyr::mutate_if(is.numeric, ~dplyr::if_else(. < !!min, 0, .))
    data[data < min] <- 0
  }

  if (!is.null(estimate_group)) {
    data <- .get.mapping.EMPT(EMPT) %>% dplyr::select(primary,!!estimate_group) %>%
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

  EMPT %<>% EMP_filter(filterFeature = id,action = 'select') %>% suppressMessages()

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
  feature <- primary <- NULL
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

    ## In case that privious estimate_group is not in the coldata,eg EMP_collapse, make the estimate_group = NULL
    if (!is.null(estimate_group)) {
      col_name <- .get.mapping.EMPT(EMPT) %>% colnames()
      if (!estimate_group %in% col_name) {
        estimate_group <- NULL
      }else {
        message('EMP_identify_assay will work according to estimate_group = ',estimate_group)
      }
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

  EMPT %<>% EMP_filter(filterFeature = id,action = 'select') %>% suppressMessages()

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

#' Identify the most core experssion and abudnace from sparse data
#'
#' @param x Object in EMPT or MultiAssayExperiment format.
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object.
#' @param estimate_group A character string. Select the group name in the coldata to be calculated. When estimate_group = NULL or "none", the function will assume all samples belong to one group.
#' @param method A character string.Methods include default, edgeR. Method default is from doi: 10.3389/fgene.2021.803627. Method edgeR in from edgeR::filterByExpr.
#' @param min A number. Set the min abundance for filtering. When method='default', min means the lowest relative bundance. When method='edgeR.', min means the lowest abosulte bundance.
#' @param min_ratio Set the min ratio presence for feature.
#' @param action A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#'
#' @return EMPT object
#' @export
#'
#' @examples
#' data(MAE)
#' ## Consider the minnum relative abundance and min ratio specilally for microbial data
#' ## First, the abundance below this threshold will be converted to 0 according to the set minimum species relative abundance. Finally, the core species will be required to meet the requirement that the occurrence rate of the species group in at least one group is higher than the preset threshold, and the remaining species will Then it will be judged as "rare species" and filtered.
#' ## Note: If absolute abundance is provided as input, it will be automatically converted to relative abundance for filtering purposes during calculations. However, the output will remain in absolute abundance.
#' MAE |>
#'   EMP_assay_extract('taxonomy') |>
#'   EMP_identify_assay(estimate_group = 'Group',method = 'default',
#'                      min=0.01,min_ratio = 0.7) # consider the Group
#' MAE |>
#'   EMP_assay_extract('taxonomy') |>
#'   EMP_identify_assay(method = 'default') # consider all samples belong to one group
#' 
#' ## Consider the minnum counts abundance and min ratio specilally for microbial data
#' MAE |>
#'   EMP_assay_extract('geno_ec') |>
#'   EMP_identify_assay(method = 'edgeR',min = 10,min_ratio = 0.7) # consider all samples belong to one group
#' 
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
