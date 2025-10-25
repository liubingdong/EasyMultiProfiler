.EMP_mutate.EMPT<- function(obj,experiment='experiment',dots=dots,.by = NULL,.before = NULL,.after = NULL,mutate_by = 'sample',location='coldata',action='colwise',keep_result=FALSE){
  
  mutate_by <- match.arg(mutate_by, c("sample", "feature"))
  location <- match.arg(location, c("coldata", "rowdata","assay"))
  action <- match.arg(action, c("colwise", "rowwise"))

  if (is(obj,"MultiAssayExperiment")) {
    if (is.null(experiment)) {
      check_obj <-"MultiAssayExperiment"
    }else{
      obj %<>% .as.EMPT(experiment = experiment)  
      check_obj <-'EMPT'
      raw_rowdata <- .get.row_info.EMPT(obj)
      raw_coldata <- .get.mapping.EMPT(obj)
      
      raw_rowdata_name <- raw_rowdata %>% colnames()
      raw_coldata_name <- raw_coldata %>% colnames()
      raw_feature_name <- raw_rowdata[['feature']]
      raw_sample_name <- raw_coldata[['primary']] 
      
      # Temporarily merge the results in the deposit for filter below
      obj %<>% .filter.merge.EMPT() %>% suppressMessages()      
      .get.method.EMPT(obj) <- 'mutate'
      #.get.history.EMPT(obj) <- call
      class(obj) <- 'EMP_assay_data'      
    }
  }else if(is(obj,'EMPT')) {
    check_obj <-'EMPT'
    .get.method.EMPT(obj) <- 'mutate'
    class(obj) <- 'EMP_assay_data' 
    #.get.history.EMPT(obj) <- call
    ## Get the initial names for recover below
    raw_rowdata <- .get.row_info.EMPT(obj)
    raw_coldata <- .get.mapping.EMPT(obj)
    
    raw_rowdata_name <- raw_rowdata %>% colnames()
    raw_coldata_name <- raw_coldata %>% colnames()
    raw_feature_name <- raw_rowdata[['feature']]
    raw_sample_name <- raw_coldata[['primary']]

    # Temporarily merge the results in the deposit for filter below
    obj <- obj |> .filter.merge.EMPT() %>% suppressMessages()
  }
  
  if (mutate_by == 'sample') {
    if (location == 'coldata') {
      deposit <- obj |> .mutate.col.EMPT(raw_coldata_name = raw_coldata_name,raw_rowdata_name=raw_rowdata_name,
        #.before=.before,.after=.after,
        action = action,dots=dots,.by = .by,.before = .before,.after = .after)
      .get.info.EMPT(deposit) <- 'EMP_mutate_col'
       
    }else if (location == 'assay') {
      obj <- obj |> .mutate.assay.feature.EMPT(raw_feature_name = raw_feature_name,
        #.before=.before,.after=.after,
        raw_rowdata_name = raw_rowdata_name,action = action,dots=dots,.by = .by,.before = .before,.after = .after)
      deposit <- obj |> .filter.deposit2.EMPT(real_sample = TRUE,real_feature = raw_feature_name,keep_result=keep_result)
      .get.info.EMPT(deposit) <- 'EMP_mutate_assay'

    }else{
      stop("When mutate_by = 'sample', location only supports 'coldata' and 'assay'!")
    }
  }

  if (mutate_by == 'feature') {
    if (location == 'rowdata') {
      deposit <- obj |> .mutate.row.EMPT(raw_coldata_name = raw_coldata_name,raw_rowdata_name = raw_rowdata_name,
        #.before=.before,.after=.after,
        action=action,dots=dots,.by = .by,.before = .before,.after = .after)
      .get.info.EMPT(deposit) <- 'EMP_mutate_row'

    }else{
      stop("When mutate_by = 'feature', location only supports 'rowdata'!")
    }
  }

  
  # delete the plot
  if (!length(deposit@plot_deposit) == 0) {
    deposit@plot_deposit <- NULL
    #EMP_message("Because condtion has changed, all plot results will be eliminated.",color = 32,order = 1,show='message')
  }
  
  deposit@deposit_append <- NULL ### WGCNA delete
  
  return(deposit)

}


.mutate.assay.feature.EMPT <- function(x,raw_feature_name,raw_rowdata_name,dots,.by = NULL,.before = NULL,.after = NULL,action='colwise') {
  action <- match.arg(action, c("colwise", "rowwise"))
  assay_df <- x |> .get.assay.EMPT()
  col_df <- x |> .get.mapping.EMPT()
  assay_df <- dplyr::full_join(assay_df,col_df,by = 'primary')
  #dots <- rlang::enquos(...)

  old_feature <- colnames(assay_df)
  if (action == 'colwise') {
    assay_df <- assay_df |> 
      dplyr::mutate(!!!dots,.by = !!.by,.before = !!.before,.after = !!.after)
  }else if (action == 'rowwise') {
    assay_df <- assay_df |> 
      dplyr::rowwise() |>
      dplyr::mutate(!!!dots,.by = !!.by,.before = !!.before,.after = !!.after)
  }else {
    stop('Parameter action only allows colwise or rowwise!')
  }
  
  new_feature <- colnames(assay_df)
  add_feature <- setdiff(new_feature,old_feature)
  new_assay_df <- assay_df |> 
    #dplyr::select(dplyr::any_of(c('primary',add_feature,raw_feature_name)))
    dplyr::select(
      dplyr::all_of(
      intersect(names(assay_df), c('primary',add_feature, raw_feature_name))
      )
    )
  .get.assay.EMPT(x) <- new_assay_df
  
  row_df <- x |> .get.row_info.EMPT() |>
    dplyr::select(dplyr::any_of(!!raw_rowdata_name))
  
  .get.row_info.EMPT(x) <- row_df                            
  return(x)
}


.mutate.col.EMPT <- function(x,raw_coldata_name,raw_rowdata_name,dots,.by = NULL,.before = NULL,.after = NULL,action='colwise') {
  action <- match.arg(action, c("colwise", "rowwise"))
  assay_df <- x |> .get.assay.EMPT()
  col_df <- x |> .get.mapping.EMPT()
  col_df <- dplyr::full_join(col_df,assay_df,by = 'primary')
  #dots <- rlang::enquos(...)
  old_colname <- colnames(col_df)

  if (action == 'colwise') {
    col_df <- col_df |> 
      dplyr::mutate(!!!dots,.by = !!.by,.before = !!.before,.after = !!.after)
  }else if (action == 'rowwise') {
    col_df <- col_df |> 
      dplyr::rowwise() |>
      dplyr::mutate(!!!dots,.by = !!.by,.before = !!.before,.after = !!.after)
  }else {
    stop('Parameter action only allows colwise or rowwise!')
  }

  new_colname <- colnames(col_df)
  new_col <- setdiff(new_colname,old_colname)
  new_col_df <- col_df |> 
    #dplyr::select(dplyr::any_of(c(new_col,raw_coldata_name))) 
    dplyr::select(
      dplyr::all_of(
      intersect(names(col_df), c(new_col, raw_coldata_name))
      )
    )    
  .get.mapping.EMPT(x) <- new_col_df
  
  row_df <- x |> .get.row_info.EMPT() |>
    dplyr::select(dplyr::any_of(!!raw_rowdata_name))
  
  .get.row_info.EMPT(x) <- row_df   
  
  return(x)
}


.mutate.row.EMPT <- function(x,raw_rowdata_name,raw_coldata_name,dots,.by = NULL,.before = NULL,.after = NULL,action='colwise') {
  action <- match.arg(action, c("colwise", "rowwise"))
  row_df <- x |> EMP_rowdata_extract()
  old_colname <- colnames(row_df)
  #dots <- rlang::enquos(...)
  if (action == 'colwise') {
    row_df <- row_df |> 
      dplyr::mutate(!!!dots,.by = !!.by,.before = !!.before,.after = !!.after)
  }else if (action == 'rowwise') {
    row_df <- row_df |> 
      dplyr::rowwise() |>
      dplyr::mutate(!!!dots,.by = !!.by,.before = !!.before,.after = !!.after)
  }else {
    stop('Parameter action only allows colwise or rowwise!')
  }
  
  new_colname <- colnames(row_df)
  new_col <- setdiff(new_colname,old_colname)

  new_row_df <- row_df |> 
    #dplyr::select(dplyr::any_of(c(raw_rowdata_name,new_col)))
    dplyr::select(
      dplyr::all_of(
      intersect(names(row_df), c(new_col, raw_rowdata_name))
      )
    )     

  .get.row_info.EMPT(x) <- new_row_df
  
  col_df <- x |> .get.mapping.EMPT() |>
    dplyr::select(dplyr::any_of(!!raw_coldata_name))
  
  .get.mapping.EMPT(x) <- col_df
  
  return(x)
}


.add_NA_col <- function(x,var_name){
  col_name <- colnames(x)
  add_name <- setdiff(var_name,col_name)
  for (i in add_name) {
    x <- x |> tibble::add_column(!!i := NA)
  }
  return(x)
}

.add_NA_row <- function(x,var_name,...){
  row_name <- x[[1]]
  add_name <- setdiff(var_name,row_name)
  new_row_df <- tibble::tibble(!!colnames(x)[1] := add_name)
  x <- x |> dplyr::bind_rows(new_row_df)
  return(x)
}


.filter.deposit2.EMPT <- function(EMPT,real_sample,real_feature,keep_result=FALSE){
  Result <- attribute <- attribute2 <- affect_when_sample_changed <- affect_when_feature_changed <- `.` <- method <- NULL
  primary <- feature <- NULL
  result_names <- names(EMPT@deposit)
  deposit_info <- .get.deposit_info.EMPT(EMPT) %>% dplyr::filter(Result %in% !!result_names)

  #total_sample_num <- dim(EMPT)[2]
  #total_feature_num <- dim(EMPT)[1]
  current_sample_name <- colnames(EMPT)
  current_feature_name <- rownames(EMPT)

  if (is.logical(real_sample)) {
    if (real_sample == TRUE) {
      check_samples <- 0
      real_sample <- colnames(EMPT)
    }else{
      #check_samples <- 1
      stop('real_sample dont allow FALSE!')
    }
  }else{
    #check_samples <- ifelse((total_sample_num - length(real_sample)) == 0,0,1)
    check_samples <- ifelse(setequal(current_sample_name,real_sample),0,1)
  }

  if (is.logical(real_feature)) {
    if (real_feature == TRUE) {
      check_features <- 0
      real_feature <- rownames(EMPT)
    }else{
      #check_features <- 1
      stop('real_feature dont allow FALSE!')
    }
  }else{
    #check_features <- ifelse((total_feature_num - length(real_feature)) == 0,0,1)
    check_features <- ifelse(setequal(current_feature_name,real_feature),0,1)
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
    
    if (all(keep_result == TRUE) | i %in% keep_result_real) {
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
          #EMPT@deposit[[i]] %<>% dplyr::filter(primary %in% real_sample)
          EMPT@deposit[[i]] %<>% .add_NA_row(var_name = current_sample_name)
        }else{
          stop("The attribute2 in the deposit_info should be normal!")
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
          #EMPT@deposit[[i]] %<>% dplyr::filter(feature %in% real_feature)
          EMPT@deposit[[i]] %<>% .add_NA_row(var_name = current_feature_name)
        }else{
          stop("The attribute2 in the deposit_info should be normal!")
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

.EMP_mutate.EMPT_m <- memoise::memoise(.EMP_mutate.EMPT,cache = cachem::cache_mem(max_size = 4096 * 1024^2))



.EMP_mutate.MAE <- function(obj,dots,.by = NULL,.before = NULL,.after = NULL,
                            mutate_by='sample',location='coldata',action='colwise',
                            keep_result=FALSE,use_cached=TRUE){
  
  if (location != 'coldata' | mutate_by != 'sample') {
    stop('EMP_mutate only suppot changing the coldata for MultiAssayExperiment!')
  }

  col_df <- obj |> EMP_coldata_extract()
  
  if (action == 'colwise') {
    col_df <- col_df |> 
      dplyr::mutate(!!!dots,.by = !!.by,.before = !!.before,.after = !!.after)
  }else if (action == 'rowwise') {
    col_df <- col_df |> 
      dplyr::rowwise() |>
      dplyr::mutate(!!!dots,.by = !!.by,.before = !!.before,.after = !!.after)
  }else {
    stop('Parameter action only allows colwise or rowwise!')
  }
  
  col_df <- col_df |> 
    tibble::column_to_rownames('primary') |>
    DataFrame()
  
  obj@colData <- col_df
  
  return(obj)
}

.EMP_mutate.MAE_m<- memoise::memoise(.EMP_mutate.MAE,cache = cachem::cache_mem(max_size = 4096 * 1024^2))


#' Create, modify, and delete columns for EMPT object
#'
#' @param obj MAE or EMPT object.
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object. 
#' @param mutate_by A character string inluding sample and feature.
#' @param location A character string inluding assay, rowdata and coldata.
#' @param .by Optionally, a selection of columns to group by for just this operation, functioning as an alternative to group_by(). See also \code{\link[dplyr]{mutate}}.
#' @param .before Optionally, control where new columns should appear.
#' @param .after (the default is to add to the right hand side). See relocate() for more details.
#' @param action A character string inluding colwise and rowwise, indicating whether operations should be performed column-wise or row-wise.
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
#' @param keep_result If the input is TRUE, it means to keep all analysis results,regardless of how samples and features change. If the input is a name, it means to keep the corresponding analysis results.
#' @param ... ... Addtional parameters, see also \code{\link[dplyr]{mutate}}.
#' @return EMPT object
#' @export
#' @section Detalis:
#' Addtional parameters inherits from the function mutate from package "dplyr", but the parameters .keep is not allowed to reset.
#'
#' @examples
#' # Create, modify, and delete the features in the assay according the data and results
#' data(MAE)
#'
#' # Create, modify, and delete the features in the assay according the data and results
#' # Caculte the F/B ratio in the assay
#' MAE |>
#'   EMP_assay_extract('taxonomy') |>
#'   EMP_collapse(collapse_by = 'row',estimate_group = 'Phylum') |>
#'   EMP_mutate(F_B_ratio = Firmicutes/Bacteroidetes,
#'              .by = primary,.before = 2,
#'              mutate_by = 'sample',location ='assay') |>
#'   EMP_filter(filterFeature = 'F_B_ratio') |>
#'   EMP_boxplot(estimate_group='Group')
#' 
#' # Add the new value in the rowdata
#' # Combine the different pvalue
#' MAE |>
#'   EMP_assay_extract('host_gene') |>
#'   EMP_diff_analysis(method = 'wilcox.test',estimate_group = 'Group') |>
#'   EMP_mutate(wilcox = pvalue,
#'              wilcox_fdr = fdr,
#'              mutate_by = 'feature',location ='rowdata') |>
#'   EMP_diff_analysis(method = 'DESeq2',.formula = ~Group) |>
#'   EMP_mutate(DESeq2 = pvalue,
#'              DESeq2_fdr = fdr,
#'              mutate_by = 'feature',location ='rowdata')
#' 
#' # Add the new value in the coldata
#' # Create the new group
#' MAE |>
#'   EMP_assay_extract('host_gene',pattern = 'A1') |>
#'   EMP_mutate(Degree = dplyr::case_when(
#'     BMI < 18.5                        ~ "Lean",
#'     BMI >= 18.5 & BMI < 24            ~ "Normal",
#'     BMI >= 24 & BMI < 28              ~ "Fat",
#'     TRUE                              ~ "Need Med"
#'     ),
#'     mutate_by = 'sample',location = 'coldata',.after = Group
#'   ) |>
#'   EMP_boxplot(estimate = 'Degree')
#'
#' # Change the whole coldata for the MultiAssayExperiment
#' MAE |>
#'   EMP_mutate(Degree = dplyr::case_when(
#'     BMI < 18.5                        ~ "Lean",
#'     BMI >= 18.5 & BMI < 24            ~ "Normal",
#'     BMI >= 24 & BMI < 28              ~ "Fat",
#'     TRUE                              ~ "Need Med"),
#'     .after = Group) |>
#'   EMP_coldata_extract()


EMP_mutate <- function(obj,experiment = NULL,...,.by = NULL,.before = NULL,.after = NULL,mutate_by='sample',location='coldata',action='colwise',
                          keep_result=FALSE,use_cached=TRUE){
  call <- match.call()
  dots <- rlang::enquos(...)
  .by <- rlang::enquo(.by)
  .before <- rlang::enquo(.before)
  .after <- rlang::enquo(.after)

  deposit <- NULL

  if (use_cached == FALSE) {
    memoise::forget(.EMP_mutate.EMPT_m) %>% invisible()
    memoise::forget(.EMP_mutate.MAE_m) %>% invisible()
  }  

  if (is(obj,"MultiAssayExperiment") & !is.null(experiment)) {
    obj <- obj |> .as.EMPT(experiment=experiment)
  }

  if (is(obj,"MultiAssayExperiment")) {
    deposit <- .EMP_mutate.MAE_m(obj=obj,dots=dots,.by = .by,.before = .before,.after = .after,mutate_by=mutate_by,location=location,action=action,
                          keep_result=keep_result)    
  }else if (is(obj,"EMPT")) {
    deposit <- .EMP_mutate.EMPT_m(obj=obj,experiment=experiment,dots=dots,.by = .by,.before = .before,.after = .after,mutate_by=mutate_by,location=location,action=action,
                          keep_result=keep_result)
    .get.history.EMPT(obj) <- call
  }else {
    stop('EMP_mutate only support MultiAssayExperiment or EMPT object!')
  }
  return(deposit)

}




