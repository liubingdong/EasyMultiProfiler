.filter.EMPT <- function(EMPT,filterSample=NULL,filterFeature=NULL,action='kick') {
  primary <- feature <- NULL
  if (action == 'kick') {
    if(!is.null(filterSample)) {
      .get.mapping.EMPT(EMPT) <- .get.mapping.EMPT(EMPT) %>%
        dplyr::filter(!primary %in% !!filterSample)
    }

    if(!is.null(filterFeature)) {
      .get.row_info.EMPT(EMPT) <- .get.row_info.EMPT(EMPT) %>%
        dplyr::filter(!feature %in% !!filterFeature)
    }
  }else if(action=='select'){
    if(!is.null(filterSample)) {
      .get.mapping.EMPT(EMPT) <- .get.mapping.EMPT(EMPT) %>%
        dplyr::filter(primary %in% !!filterSample)
    }

    if(!is.null(filterFeature)) {
      .get.row_info.EMPT(EMPT) <- .get.row_info.EMPT(EMPT) %>%
        dplyr::filter(feature %in% !!filterFeature)
    }
  }else {
    warning('action must be one of kick or select!')
  }
  return(EMPT)
}


#' Filer experssion or abundance data that match a condition
#'
#' @param obj EMPT object.
#' @param sample_condition Expressions that return a logical value, and are defined in terms of the variables in coldata. If multiple expressions are included, they are combined with the &，| operator. 
#' @param feature_condition Expressions that return a logical value, and are defined in terms of the variables in rowdata. If multiple expressions are included, they are combined with the &，| operator. 
#' @param filterSample A series of character strings. Select samples in the data exactly.
#' @param filterFeature A series of character strings. Select samples in the data exactly.
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object. 
#' @param show_info A character string. Set the class of EMPT to show properly.
#' @param action A character string. You can use the filterSample and filterFeature parameters in conjunction with this. The choice is whether to keep filterSample and filterFeature (select), or simply exclude them (kick).
#'
#' @return EMPT object
#' @export
#'
#' @examples
#' data(MAE)
#' ## from MultiAssayExperiment
#' MAE |>
#'   EMP_filter(sample_condition = BMI>20 & Sex == 'M') |>
#'   EMP_summary()
#' 
#' ## from EMPT
#' MAE |>
#'   EMP_assay_extract(experiment = 'host_gene') |>
#'   EMP_identify_assay(method = 'edgeR',estimate_group = 'Group') |>
#'   EMP_diff_analysis(method='DESeq2',.formula = ~Group) |>
#'   EMP_filter(feature_condition = fdr < 0.05)
#'   
#' MAE |>
#'   EMP_assay_extract(experiment = 'taxonomy') |>
#'   EMP_alpha_analysis() |>
#'   EMP_filter(sample_condition = shannon >3 | invsimpson >19)
#' 
#' ## Precise selection
#' MAE |>
#'   EMP_assay_extract(experiment = 'taxonomy') |>
#'   EMP_alpha_analysis() |>
#'   EMP_filter(sample_condition = shannon >3 | invsimpson >19,
#'              filterSample = c('P11774','P31579'),action = 'kick') # Accurately kick samples based on satisfying sample_condition

EMP_filter <- function(obj,sample_condition,feature_condition,
                       filterSample=NULL,filterFeature=NULL,experiment=NULL,
                       show_info=NULL,action='select'){
  primary <- feature <- NULL
  call <- match.call()
  sample_condition <- dplyr::enquo(sample_condition)
  feature_condition <- dplyr::enquo(feature_condition)
  message_info <- list()
  if (inherits(obj,"MultiAssayExperiment")) {
    if (is.null(experiment)) {
      check_obj <-"MultiAssayExperiment"
    }else{
      obj %<>% .as.EMPT(experiment = experiment)  
      check_obj <-'EMPT'
      .get.method.EMPT(obj) <- 'filter'
      .get.history.EMPT(obj) <- call
      class(obj) <- 'EMP_assay_data'
      .get.info.EMPT(obj) <- 'EMP_assay_data'
    }
  }else if(inherits(obj,'EMPT')) {
    check_obj <-'EMPT'
    .get.method.EMPT(obj) <- 'filter'
    .get.history.EMPT(obj) <- call
    ## Get the initial names for recover below
    raw_rowdata_name <- .get.row_info.EMPT(obj) %>% colnames()
    raw_coldata_name <- .get.mapping.EMPT(obj) %>% colnames()

    # Temporarily merge the results in the deposit for filter below
    obj %<>% .filter.merge.EMPT() %>% suppressMessages()
  }

  total_sample_num <- obj %>% EMP_coldata_extract() %>% nrow()
  total_feature_num <- obj %>% EMP_rowdata_extract() %>% nrow()

  ## select sample
  obj %>% EMP_coldata_extract() %>%
    dplyr::filter(!!sample_condition) %>% dplyr::pull(primary) -> sample_id

  ## select feature
  obj %>% EMP_rowdata_extract() %>%
    dplyr::filter(!!feature_condition) %>% dplyr::pull(feature) -> feature_id

  if (action == 'select') {
    if (!is.null(filterSample)) {
      real_sample <- intersect(sample_id,filterSample)
    }else {
      real_sample <- sample_id
    }
    if (!is.null(filterFeature)) {
      real_feature <- intersect(feature_id,filterFeature)
    }else {
      real_feature <- feature_id
    }
  }else if(action == 'kick') {
    if (!is.null(filterSample)) {
      real_sample <- setdiff(sample_id,filterSample)
    }else {
      real_sample <- sample_id
    }
    if (!is.null(filterFeature)) {
      real_feature <- setdiff(feature_id,filterFeature)
    }else {
      real_feature <- feature_id
    }
  }else{
    stop('Paramter action should be select or kick!')
  }

  sample_filter_num<- total_sample_num - length(real_sample)
  feature_filter_num <- total_feature_num - length(real_feature)


  if (length(real_sample) == 0 | length(real_feature) == 0) {
    warning('No sample or feature meet the condition, plz reset the parameter!')
  }else {
    if (check_obj == "MultiAssayExperiment") {
      deposit <- obj[real_feature,real_sample,,drop=FALSE]
  }else if (check_obj == 'EMPT') {
        message_info %<>% append(paste0(sample_filter_num,' of ',total_sample_num,' samples were filterd out!'))
        message_info %<>% append(paste0(feature_filter_num,' of ',total_feature_num,' features were filterd out!'))

        # filter the result in the deposit
        obj %<>% .filter.deposit.EMPT(real_sample,real_feature)

        .get.message_info.EMPT(obj) <- message_info
        deposit <- obj %>% .filter.EMPT(filterSample=real_sample,filterFeature=real_feature,action='select') %>% suppressMessages() ##  Here action must be select, dont change!
        ### delete merged data used in the previous filter, and make sure the clean data
        .get.mapping.EMPT(deposit) <- .get.mapping.EMPT(deposit) %>% dplyr::select(dplyr::all_of(raw_coldata_name))
        .get.row_info.EMPT(deposit) <- .get.row_info.EMPT(deposit) %>% dplyr::select(dplyr::all_of(raw_rowdata_name))


        # delete the plot
        if (!length(deposit@plot_deposit) == 0) {
          deposit@plot_deposit <- NULL
          message("Because condtion has changed, all plot results will be removed!")
        }

        deposit@deposit_append <- NULL ### WGCNA delete
        ## select the show information
        if(!is.null(show_info)){
          class(deposit) <- show_info
          .get.info.EMPT(deposit) <- show_info
        }
        ## When result is empty, make sure the output is ok
        check_result_empty <- .get.result.EMPT(deposit) %>% suppressMessages() %>% length() == 0 | .get.result.EMPT(deposit) %>% suppressMessages() %>% is.null()
        if (check_result_empty) {
           class(deposit) <- 'EMP_assay_data'
          .get.info.EMPT(deposit) <- 'EMP_assay_data'         
        }        
  }else{
        stop('Input data should be EMP or EMPT!')
  }

  return(deposit)
  }
}




.filter.deposit.EMPT <- function(EMPT,real_sample,real_feature){
  Result <- attribute <- attribute2 <- affect_when_sample_changed <- affect_when_feature_changed <- `.` <- NULL
  primary <- feature <- NULL
  result_names <- names(EMPT@deposit)
  deposit_info <- .get.deposit_info.EMPT(EMPT) %>% dplyr::filter(Result %in% !!result_names)

  total_sample_num <- EMPT %>% EMP_coldata_extract() %>% nrow()
  total_feature_num <- EMPT %>% EMP_rowdata_extract() %>% nrow()

  check_samples <- ifelse((total_sample_num - length(real_sample)) == 0,0,1)
  check_features <- ifelse((total_feature_num - length(real_feature)) == 0,0,1)

  # for debug
  #print(check_samples)
  #print(check_features)

  for (i in result_names) {
    each_deposit_info <- deposit_info %>% dplyr::filter(Result == !!i)
    result_attribute <- each_deposit_info %>% dplyr::pull(attribute)
    result_attribute2 <- each_deposit_info %>% dplyr::pull(attribute2)
    result_source <- each_deposit_info %>% dplyr::pull(source)
    ## check the affect
    temp <- each_deposit_info %>% dplyr::select(affect_when_sample_changed,affect_when_feature_changed)
    temp[2,] <- c(check_samples,check_features)
    affect_status <- any(colSums(temp) == 2)
    ### confirm the cause for message below
    affect_cause <- c('samples','features')[colSums(temp) ==2]
    if(length(affect_cause) == 2) {
      affect_cause <- "samples and features"
    }

    if (result_attribute == 'primary') {
      if (!affect_status) {
        if (result_attribute2 == 'normal') {
          EMPT@deposit[[i]] %<>% dplyr::filter(primary %in% real_sample)
        }else if(result_attribute2 == 'diagonal'){
          # EMPT@deposit[[i]] %<>% as.data.frame() %>%
          #   dplyr::select(dplyr::all_of(!!real_sample)) %>%
          #   dplyr::filter(rownames(.) %in% !!real_sample)

          EMPT@deposit[[i]] %<>% as.data.frame() %>%
            dplyr::select(dplyr::all_of(!!real_sample))
          EMPT@deposit[[i]] <- dplyr::filter(rownames(EMPT@deposit[[i]]) %in% !!real_sample)
        }else{
          stop("The attribute2 in the deposit_info should be normal or diagonal!")
        }
      }else{
        EMPT@deposit[[i]] <- NULL
        message_wrap("If any ",affect_cause, " in the experiment have changed, the ",i,
                     " will become NULL and ",result_source, " should be re-run if needed.")
      }
    }else if(result_attribute == 'feature'){
      if (!affect_status) {
        if (result_attribute2 == 'normal') {
          EMPT@deposit[[i]] %<>% dplyr::filter(feature %in% real_feature)
        }else if(result_attribute2 == 'diagonal'){
          EMPT@deposit[[i]] %<>% as.data.frame() %>%
            dplyr::select(dplyr::all_of(!!real_feature)) %>%
            dplyr::filter(rownames(.) %in% !!real_feature)
        }else{
          stop("The attribute2 in the deposit_info should be normal or diagonal!")
        }
      }else{
        EMPT@deposit[[i]] <- NULL
        message_wrap("If any ",affect_cause, " in the experiment have changed, the ",i,
                     " will become NULL and ",result_source, " should be re-run if needed.")
      }
    }else if(result_attribute == 'all'){
      EMPT@deposit[[i]] <-NULL
      message_wrap("If any ",affect_cause, " in the experiment have changed, the ",i,
                   " will become NULL and ",result_source, " should be re-run if needed.")
    }else{
      stop("The attribute in the deposit_info should be primary, feature or all!")
    }
  }
  if(is.null(.get.info.EMPT(EMPT))) {
    class(EMPT) <- 'EMP_assay_data'
    .get.info.EMPT(EMPT) <- 'EMP_assay_data'
  }
  return(EMPT)
}


.filter.merge.EMPT <- function(EMPT){
  Result <- attribute <- attribute2 <- NULL
  result_names <- names(EMPT@deposit)
  deposit_info <- .get.deposit_info.EMPT(EMPT) %>%
    dplyr::filter(Result %in% !!result_names ) %>%
    dplyr::filter(attribute %in% c('primary','feature')) %>%
    dplyr::filter(attribute2 == 'normal')

  primary_result_name <- deposit_info %>% dplyr::filter(attribute =='primary') %>% dplyr::pull(Result)
  feature_result_name <- deposit_info %>% dplyr::filter(attribute =='feature') %>% dplyr::pull(Result)

  for (i in primary_result_name) {
    .get.mapping.EMPT(EMPT) <-  dplyr::left_join(.get.mapping.EMPT(EMPT),EMPT@deposit[[i]],by='primary')
  }
  for (j in feature_result_name) {
    .get.row_info.EMPT(EMPT) <-  dplyr::left_join(.get.row_info.EMPT(EMPT),EMPT@deposit[[j]],by='feature')
  }

  return(EMPT)
}
