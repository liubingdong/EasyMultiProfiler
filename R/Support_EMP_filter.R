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
    warning('Parameter action must be one of kick or select!')
  }
  return(EMPT)
}

.EMP_filter <- function(obj,sample_condition,feature_condition,
                       filterSample=NULL,filterFeature=NULL,experiment=NULL,
                       show_info=NULL,action='select',remove_zero_sum_sample=FALSE,remove_zero_sum_feature=FALSE,keep_result=FALSE){
  primary <- feature <- NULL
  #call <- match.call()
  sample_condition <- dplyr::enquo(sample_condition)
  feature_condition <- dplyr::enquo(feature_condition)
  message_info <- list()
  if (is(obj,"MultiAssayExperiment")) {
    if (is.null(experiment)) {
      check_obj <-"MultiAssayExperiment"
    }else{
      obj %<>% .as.EMPT(experiment = experiment)  
      check_obj <-'EMPT'
      raw_rowdata_name <- .get.row_info.EMPT(obj) %>% colnames()
      raw_coldata_name <- .get.mapping.EMPT(obj) %>% colnames()

      # Temporarily merge the results in the deposit for filter below
      obj %<>% .filter.merge.EMPT() %>% suppressMessages()      
      .get.method.EMPT(obj) <- 'filter'
      #.get.history.EMPT(obj) <- call
      class(obj) <- 'EMP_assay_data'
      .get.info.EMPT(obj) <- 'EMP_assay_data'
    }
  }else if(is(obj,'EMPT')) {
    check_obj <-'EMPT'
    .get.method.EMPT(obj) <- 'filter'
    #.get.history.EMPT(obj) <- call
    ## Get the initial names for recover below
    raw_rowdata_name <- .get.row_info.EMPT(obj) %>% colnames()
    raw_coldata_name <- .get.mapping.EMPT(obj) %>% colnames()

    # Temporarily merge the results in the deposit for filter below
    obj %<>% .filter.merge.EMPT() %>% suppressMessages()
  }

  total_sample_num <- dim(obj)[2]
  total_feature_num <- dim(obj)[1]

  ## select sample
  obj %>% EMP_coldata_extract() %>%
    dplyr::filter(!!sample_condition) %>% 
    dplyr::pull(primary) -> sample_id

  ## select feature
  obj %>% EMP_rowdata_extract() %>%
    dplyr::filter(!!feature_condition) %>% 
    dplyr::pull(feature) -> feature_id

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

  # Remove zero-sum features
  if (remove_zero_sum_sample == TRUE | remove_zero_sum_feature == TRUE) {

    if (!.get.assay_name.EMPT(obj) %in% c('counts','relative','integer')) {
      EMP_message("This assay is not count, relative abundance, or integer data. Thus, using the remove_zero_sum parameter is not recommended!",color = 31,order = 1,show='warning')
    }

    assay_data<- assay(obj)[real_feature,real_sample]
    
    check_zero_result <- check_zero(mat = assay_data,
                                  check_zero_sum_sample=remove_zero_sum_sample,
                                  check_zero_sum_feature=remove_zero_sum_feature)
      
    while (check_zero_result$flag==TRUE) {
      filter_assay_result <- filter_zero_martrx(mat = assay_data,
                                                remove_zero_sum_sample = remove_zero_sum_sample,
                                                remove_zero_sum_feature = remove_zero_sum_feature,
                                                real_sample = real_sample,
                                                real_feature = real_feature,
                                                filterSample = check_zero_result$no_zero_sample,
                                                filterFeature = check_zero_result$no_zero_feature
                                                )
      assay_data <- filter_assay_result$final_mat
      real_sample <- filter_assay_result$real_sample
      real_feature <- filter_assay_result$real_feature
      check_zero_result <- check_zero(assay_data,
                                    check_zero_sum_sample=remove_zero_sum_sample,
                                    check_zero_sum_feature=remove_zero_sum_feature)
    }
  }

  sample_filter_num<- total_sample_num - length(real_sample)
  feature_filter_num <- total_feature_num - length(real_feature)

  if (length(real_sample) == 0 | length(real_feature) == 0) {
    warning("No sample or feature meet the condition, plz reset the parameter!")
  }else {
    if (check_obj == "MultiAssayExperiment") {
      deposit <- obj[real_feature,real_sample,,drop=FALSE]
  }else if (check_obj == 'EMPT') {
        message_info %<>% append(paste0(sample_filter_num,' of ',total_sample_num,' samples were filterd out!'))
        message_info %<>% append(paste0(feature_filter_num,' of ',total_feature_num,' features were filterd out!'))

        # filter the result in the deposit
        obj %<>% .filter.deposit.EMPT(real_sample,real_feature,keep_result=keep_result) 

        .get.message_info.EMPT(obj) <- message_info
        deposit <- obj %>% .filter.EMPT(filterSample=real_sample,filterFeature=real_feature,action='select') %>% suppressMessages() ##  Here action must be select, dont change!
        ### delete merged data used in the previous filter, and make sure the clean data
        .get.mapping.EMPT(deposit) <- .get.mapping.EMPT(deposit) %>% dplyr::select(dplyr::all_of(raw_coldata_name))
        .get.row_info.EMPT(deposit) <- .get.row_info.EMPT(deposit) %>% dplyr::select(dplyr::all_of(raw_rowdata_name))


        # delete the plot
        if (!length(deposit@plot_deposit) == 0) {
          deposit@plot_deposit <- NULL
          #EMP_message("Because condtion has changed, all plot results will be eliminated.",color = 32,order = 1,show='message')
        }

        deposit@deposit_append <- NULL ### WGCNA delete
        ## select the show information
        if(!is.null(show_info)){
          class(deposit) <- show_info
          .get.info.EMPT(deposit) <- show_info
        }
        if(is(.get.result.EMPT(deposit) %>% spsUtil::quiet(),'tibble') | is(.get.result.EMPT(deposit) %>% spsUtil::quiet(),'data.frame')) {
          deposit <- deposit
        }else{
          check_result_empty <- .get.result.EMPT(deposit) %>% spsUtil::quiet() %>% is.null() ||
            .get.result.EMPT(deposit) %>% spsUtil::quiet() %>% length() == 0 ||
            all(.get.result.EMPT(deposit) %>% spsUtil::quiet()  == "No info is matched!")
          class(deposit) <- 'EMP_assay_data'
          .get.info.EMPT(deposit) <- 'EMP_assay_data' 
        }    
  }else{
        stop('Input data should be EMP or EMPT!')
  }

  return(deposit)
  }
}


.EMP_filter_m <- memoise::memoise(.EMP_filter,cache = cachem::cache_mem(max_size = 4096 * 1024^2))


.filter.deposit.EMPT <- function(EMPT,real_sample,real_feature,keep_result=FALSE){
  Result <- attribute <- attribute2 <- affect_when_sample_changed <- affect_when_feature_changed <- `.` <- method <- NULL
  primary <- feature <- NULL
  result_names <- names(EMPT@deposit)
  deposit_info <- .get.deposit_info.EMPT(EMPT) %>% dplyr::filter(Result %in% !!result_names)

  total_sample_num <- dim(EMPT)[2]
  total_feature_num <- dim(EMPT)[1]

  if (is.logical(real_sample)) {
    if (real_sample == TRUE) {
      check_samples <- 0
    }else{
      check_samples <- 1
    }
  }else{
    check_samples <- ifelse((total_sample_num - length(real_sample)) == 0,0,1)
  }

  if (is.logical(real_feature)) {
    if (real_feature == TRUE) {
      check_features <- 0
    }else{
      check_features <- 1
    }
  }else{
    check_features <- ifelse((total_feature_num - length(real_feature)) == 0,0,1)
  }

  # for debug
  #print(check_samples)
  #print(check_features)
  
  sample_affect_result <- c()
  feature_affect_result <- c()
  all_affect_result <- c()

  for (i in result_names) {
    each_deposit_info <- deposit_info %>% dplyr::filter(Result == !!i)
    
    # check the keep result name
    if (!is.logical(keep_result)) {
      if (each_deposit_info$Result %in% keep_result) {
         keep_result_real<- keep_result
      }else{
         keep_result_real <- each_deposit_info %>% dplyr::filter(source %in% keep_result) %>% dplyr::pull(Result)
      }
    }else{
      keep_result_real <- keep_result
    }

    ## Special cases in the EMP_diff_analysis
    if (i == "diff_analysis_result") {
      diff_method <- .get.result.EMPT(EMPT,info = 'EMP_diff_analysis') %>% spsUtil::quiet() %>% dplyr::pull(method) %>% unique()
      affect_diff_method <- c('edgeR_quasi_likelihood', 'edgeR_likelihood_ratio', 'edger_robust_likelihood_ratio', 'DESeq2',
                              'limma_voom',  'limma_voom_sample_weights')
      # wilcox.test may have continuity correction or exact test,so need all() to lead one match result
      if (all(diff_method %in% affect_diff_method)) {
        each_deposit_info$affect_when_sample_changed <- 1
        each_deposit_info$affect_when_feature_changed <- 1
      } 
    }
    
    if (keep_result == TRUE | i %in% keep_result_real) {
       each_deposit_info$affect_when_sample_changed <- 0
       each_deposit_info$affect_when_feature_changed <- 0
    }

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
        # info_output <- paste0("If any ",affect_cause, " have changed, the ",i,
        #              " will become NULL.\n",result_source, " should be re-run if needed.")
        # EMP_message(info_output,color = 32,order = 1,show='message')
        sample_affect_result <- append(sample_affect_result,i)
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
        # info_output <- paste0("If any ",affect_cause, " have changed, the ",i,
        #              " will become NULL.\n",result_source, " should be re-run if needed.")
        # EMP_message(info_output,color = 32,order = 1,show='message')
        feature_affect_result <- append(feature_affect_result,i)
      }
    }else if(result_attribute == 'all'){
      EMPT@deposit[[i]] <-NULL
      # info_output <- paste0("If any ",affect_cause, " have changed, the ",i,
      #              " will become NULL.\n",result_source, " should be re-run if needed.")
      # EMP_message(info_output,color = 32,order = 1,show='message')
      all_affect_result <- append(all_affect_result,i)
    }else{
      stop("The attribute in the deposit_info should be primary, feature or all!")
    }
  }

  sample_affect_analysis <- deposit_info |>
    dplyr::filter(Result %in% {{sample_affect_result}}) |>
    dplyr::pull(source) |>
    unique()
  #if (length(sample_affect_analysis) > 0) {
  #    sample_affect_analysis_info <- paste0("Due to variations in the samples, The results from ", .concat_str(sample_affect_analysis)," will be eliminated.")
  #    EMP_message(sample_affect_analysis_info,color = 32,order = 1,show='message')
  #}

  feature_affect_analysis <- deposit_info |>
    dplyr::filter(Result %in% {{feature_affect_result}}) |>
    dplyr::pull(source) |>
    unique()
  #if (length(feature_affect_analysis) > 0) {
  #    feature_affect_analysis_info <- paste0("Due to variations in the features, The results from ", .concat_str(feature_affect_analysis)," will be eliminated.")
  #    EMP_message(feature_affect_analysis_info,color = 32,order = 1,show='message')
  #}


  all_affect_analysis <- deposit_info |>
    dplyr::filter(Result %in% {{all_affect_result}}) |>
    dplyr::pull(source) |>
    unique()
  #if (length(all_affect_analysis) > 0) {
  #    all_affect_analysis_info <- paste0("Due to variations in the samples or features, The results from ", .concat_str(all_affect_analysis)," will be eliminated.")
  #    EMP_message(all_affect_analysis_info,color = 32,order = 1,show='message')
  #}

  total_assay_analysis <- c(sample_affect_analysis,feature_affect_analysis,all_affect_analysis) |> unique()
  if (length(total_assay_analysis) > 0) {
      total_assay_analysis_info <- paste0("Due to changes in the samples or features, The results from ", .concat_str(total_assay_analysis)," will be eliminated.")
      EMP_message(total_assay_analysis_info,color = 32,order = 1,show='message')
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


.concat_str <- function(str_vector){
  if (length(str_vector) > 1) {
    result <- paste(str_vector[-length(str_vector)], collapse = ", ")
    result <- paste(result, "and", str_vector[length(str_vector)], sep = " ")
  } else {
    result <- str_vector
  }
  return(result)
}



check_zero_martrix <- function(mat, MARGIN,tol = 1e-10) {
  # 处理空矩阵的情况
  deposit <- list()
  deposit$flag <- FALSE
  deposit$zero_name < NULL
  deposit$no_zero_name < NULL
  
  if (nrow(mat) == 0 | ncol(mat) == 0) {
    stop("Empty assay happend!")
  }
  
  # 计算行总和
  sums <- apply(mat, MARGIN, function(x) sum(x))
  # 使用容差比较，避免浮点数精度问题
  switch (MARGIN,
    `1` = {zero_name <- rownames(mat)[sums < tol]
         deposit$zero_name <- zero_name
         deposit$no_zero_name <- setdiff(rownames(mat),zero_name)
        },
    `2` = {zero_name <- colnames(mat)[sums < tol]
         deposit$zero_name <- zero_name
         deposit$no_zero_name <- setdiff(colnames(mat),zero_name)
    }
  )  
  if (length(deposit$zero_name) > 0) {
    deposit$flag <- TRUE
  }
  deposit$status <- MARGIN
  return(deposit)
}

filter_zero_martrx <- function(mat,real_sample,real_feature,
                               remove_zero_sum_sample=FALSE,remove_zero_sum_feature=FALSE,
                               filterSample=NULL,filterFeature=NULL){
  deposit <- list()
  assay_data<- mat[real_feature,real_sample]
  if (remove_zero_sum_sample == TRUE) {
    if (is.null(filterSample)) {
      no_zero_sample <- rownames(assay_data)[apply(assay_data, 2, function(x) sum(x) != 0)]
    }else{
      no_zero_sample <- filterSample
    }
    real_sample <- intersect(real_sample,no_zero_sample)
  }
  
  if (remove_zero_sum_feature == TRUE) {
    if (is.null(filterFeature)) {
      no_zero_feature <- rownames(assay_data)[apply(assay_data, 1, function(x) sum(x) != 0)]
    }else{
      no_zero_feature <- filterFeature
    }
    real_feature <- intersect(real_feature,no_zero_feature)
  }     
  deposit$real_sample <- real_sample
  deposit$real_feature <- real_feature
  deposit$final_mat <- mat[real_feature,real_sample]
  return(deposit)
}


check_zero <- function(mat,check_zero_sum_sample,check_zero_sum_feature){
  deposit <- list()
  check_zero_sample <- list()
  check_zero_feature <- list()
  check_zero_sample$flag <- FALSE
  check_zero_feature$flag <- FALSE

  if (check_zero_sum_sample==TRUE) {
    check_zero_sample <- check_zero_martrix(mat,2)
    deposit$no_zero_sample <- check_zero_sample$no_zero_name
  }
  if (check_zero_sum_feature==TRUE) {
    check_zero_feature <- check_zero_martrix(mat,1)
    deposit$no_zero_feature <- check_zero_feature$no_zero_name
  }

  deposit$flag <- any(check_zero_sample$flag,check_zero_feature$flag)

  return(deposit)
}




#' Filer experssion or abundance data that match a condition
#'
#' @param obj EMPT object.
#' @param sample_condition Expressions that return a logical value, and are defined in terms of the variables in coldata. If multiple expressions are included, they are combined with the &，| operator. 
#' @param feature_condition Expressions that return a logical value, and are defined in terms of the variables in rowdata. If multiple expressions are included, they are combined with the &，| operator. 
#' @param filterSample A series of character strings. Select samples in the data exactly.
#' @param filterFeature A series of character strings. Select samples in the data exactly.
#' @param remove_zero_sum_sample Remove samples with zero counts in all features.(default:FALSE)
#' @param remove_zero_sum_feature Remove features with zero counts in all samples.(default:FALSE)
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object. 
#' @param keep_result If the input is TRUE, it means to keep all analysis results,regardless of how samples and features change. If the input is a name, it means to keep the corresponding analysis results.
#' @param show_info A character string. Set the class of EMPT to show properly.
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
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
                       show_info=NULL,action='select',remove_zero_sum_sample=FALSE,remove_zero_sum_feature=FALSE,keep_result=FALSE,use_cached=TRUE){
  call <- match.call()
  deposit <- NULL
  sample_condition <- dplyr::enquo(sample_condition)
  feature_condition <- dplyr::enquo(feature_condition)
  if (use_cached == FALSE) {
    memoise::forget(.EMP_filter_m) %>% invisible()
  }  
  deposit <- .EMP_filter_m(obj=obj,sample_condition={{sample_condition}},feature_condition={{feature_condition}},
                       filterSample=filterSample,filterFeature=filterFeature,
                       remove_zero_sum_sample=remove_zero_sum_sample,remove_zero_sum_feature=remove_zero_sum_feature,
                       experiment=experiment,show_info=show_info,action=action,keep_result=keep_result)
  if (is(obj,"MultiAssayExperiment")) {
    return(deposit)
  }else if (is(obj,"EMPT")) {
    .get.history.EMPT(obj) <- call
    return(deposit)
  }
}

