#' @importFrom rlang `:=`
#' @importFrom dplyr last_col
#' @importFrom dplyr matches
#' @importFrom dplyr n_distinct
#' @importFrom SummarizedExperiment `rowData<-`
#' @importFrom SummarizedExperiment `colData<-`
.EMP_diff_analysis_tidybulk <- function(EMPT,method,.formula,p.adjust='fdr',group_level=NULL,...) {
  Estimate_group <- pvalue <- feature <- sign_group <- vs <- log2FC <- estimate_group <- fold_change <- `.` <- tiny_EMPT <- NULL
  batch_effect <- NULL
  Group_info <- as.list(.formula)[[2]]
  estimate_group <- as.list(Group_info)[[length(Group_info)]] %>% as.character()
  if (length(as.list(Group_info)) ==3) {
    batch_info <- as.list(Group_info)[[2]] %>% as.character()
    check_batch_info <- batch_info %>% as.numeric() %>% suppressWarnings()
    if (!is.na(check_batch_info)) {
      batch_info <- NULL
    }
  }else{
    batch_info <- NULL
  }

  if(length(group_level) != 2 & !is.null(group_level)){
    stop("The group_level parameter must have exactly 2 factors!")
  }

  ## reduce the cost when using tidybulk::tidybulk
  tiny_EMPT <- EMPT

  new_coldata <- colData(tiny_EMPT) 
  new_coldata <- new_coldata[,c(estimate_group,batch_info),drop=FALSE]

  new_rowdata <- rowData(tiny_EMPT)
  new_rowdata <- new_rowdata[,1,drop = FALSE]

  ## check the missing value in the group label
  if(any(is.na(new_coldata[[estimate_group]]))) {
    stop('Column ',estimate_group,' has beed deteced missing value, please check and filter them!')
  }

  rowData(tiny_EMPT) <- new_rowdata
  colData(tiny_EMPT) <- new_coldata

  melt_EMPT <- tiny_EMPT %>% tidybulk::tidybulk()

  check_group <- melt_EMPT %>% dplyr::pull(estimate_group) %>% dplyr::n_distinct() == 2
  if(!check_group ) {
    stop('For ',method,' only support supports two-category group!')
  }

  origin_group_level <- group_level

  if (!is.null(group_level)) {
    switch(method,
           "edgeR_quasi_likelihood" = {group_level <- paste0(estimate_group,group_level[1],'-',estimate_group,group_level[2])},
           "edgeR_likelihood_ratio"  = {group_level <- paste0(estimate_group,group_level[1],'-',estimate_group,group_level[2])},
           "edger_robust_likelihood_ratio"  = {group_level <- paste0(estimate_group,group_level[1],'-',estimate_group,group_level[2])},
           "DESeq2"  = {group_level <- list(c(estimate_group,group_level[1],group_level[2])) },
           "limma_voom"  = {group_level <- paste0(estimate_group,group_level[1],'-',estimate_group,group_level[2])},
           "limma_voom_sample_weights"  = {group_level <- paste0(estimate_group,group_level[1],'-',estimate_group,group_level[2])},
           {
             group_level <- group_level
           }
    )
  }

  result <- melt_EMPT %>%
    tidybulk::test_differential_abundance(
      .formula,
      method = method,action='get',contrasts=group_level,...) %>%
    suppressWarnings()  %>% suppressMessages()


  if (!is.null(group_level)) {
    #  design_raw_info <- result %>%
    #    dplyr::select(last_col()) %>% colnames() %>%
    #    strsplit(.,split = '___') %>% unlist()
    #  
    #  design_fix <- paste0('___',design_raw_info[2])
    #  
    #  result %<>%
    #    dplyr::rename_with(~stringr::str_remove(., design_fix))
    #  
    #  design_info <- design_raw_info[2]
    #  design_info %<>% gsub(estimate_group, "",.) %>% trimws() %>%
    #    gsub("-", " vs ",.)
    design_info <- paste0(origin_group_level[1],' vs ', origin_group_level[2])
    design_info_detail <- strsplit(design_info,'vs') %>% unlist() %>% trimws()
  }else {
    design_raw_info <- factor(colData(EMPT)[[estimate_group]]) %>% unique()
    design_info <- paste0(design_raw_info[2],' vs ',design_raw_info[1])
    design_info_detail <- rev(design_raw_info)
  }

  search_log2FC <- c("logFC", "log2FoldChange")
  search_pvalue <- c("p.Value", "Pvalue")

  result %<>%
    dplyr::rename_with(~ "log2FC", matches(paste(search_log2FC, collapse = "|"))) %>%
    dplyr::rename_with(~ "pvalue", matches(paste(search_pvalue, collapse = "|"))) %>%
    dplyr::select(feature,log2FC,pvalue) %>%
    dplyr::mutate(
      sign_group = dplyr::case_when(
        log2FC >= 0 ~ !!design_info_detail[1],
        log2FC < 0 ~ !!design_info_detail[2]
      )
    ) %>%
    dplyr::mutate(!!p.adjust := stats::p.adjust(pvalue,method = p.adjust)) %>%
    dplyr::mutate(Estimate_group = !!estimate_group) %>%
    dplyr::mutate(vs = !!design_info) %>%
    dplyr::mutate(method = !!method) %>%
    dplyr::mutate(fold_change = 2**log2FC) %>%
    dplyr::select(feature,Estimate_group,pvalue,!!p.adjust,sign_group,method,vs,fold_change,log2FC)

  if(method == 'DESeq2' & !is.null(batch_info)) {
    result %<>% dplyr::mutate(batch_effect=batch_info) %>%
      dplyr::select(feature,Estimate_group,batch_effect,everything())
  }

  EMPT@deposit[['diff_analysis_result']] <- result
  .get.estimate_group.EMPT(EMPT) <- estimate_group
  .get.estimate_group_info.EMPT(EMPT) <- design_info
  .get.formula.EMPT(EMPT) <- .formula
  .get.method.EMPT(EMPT) <- method
  .get.algorithm.EMPT(EMPT) <- 'diff_analysis'
  .get.info.EMPT(EMPT) <- 'EMP_diff_analysis'
  EMPT
}

.EMP_diff_analysis_tidybulk_m <- memoise::memoise(.EMP_diff_analysis_tidybulk,cache = cachem::cache_mem(max_size = 2048 * 1024^2))



#' Differential expression or abundance analysis
#'
#' @param obj Object in EMPT or MultiAssayExperiment format.
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object.
#' @param method A character string. Methods include t.test, wilcox_test, kruskal.test, oneway.test, edgeR_quasi_likelihood, edgeR_likelihood_ratio, edger_robust_likelihood_ratio, DESeq2, limma_voom, limma_voom_sample_weights
#' @param p.adjust A character string. Adjust P-values for Multiple Comparisons inluding fdr, holm, hochberg, hommel, bonferroni, BH, BY. (default:fdr)
#' @param .formula A formula representing the desired linear model. If there is more than one factor, they should be in the order factor of interest + additional factors.
#' @param estimate_group A character string. Select the group name in the coldata to be calculated.
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
#' @param action A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param group_level A series of character strings. Determine the comparison order of groups.
#' @param core A number. Select the core numbers in the parallel compute to speed up the result.
#' @param ... Further parameters passed to the function tidybulk::test_differential_abundance, or statistical function in the stats package.
#' @importFrom memoise forget
#'
#' @return EMPT object
#' @export
#'
#' @examples
#' data(MAE)
#' ## t.test or wilcox.test
#' MAE |>
#'   EMP_decostand(experiment = 'taxonomy',method = 'relative',pseudocount=0.0001) |>
#'   EMP_diff_analysis(method = 't.test',estimate_group = 'Group',p.adjust = 'fdr')
#' 
#' MAE |>
#'   EMP_assay_extract(experiment = 'taxonomy') |>
#'   EMP_diff_analysis(method = 'wilcox.test',estimate_group = 'Group',p.adjust = 'BH')
#' 
#' ## DESeq2
#' MAE |>
#'   EMP_decostand(experiment = 'geno_ec',method = 'integer') |>
#'   EMP_diff_analysis(method='DESeq2',.formula = ~Group)  
#' 
#' MAE |>
#'   EMP_decostand(experiment = 'geno_ec',method = 'integer') |>
#'   EMP_diff_analysis(method='DESeq2',.formula = ~Region+Group)  ## Eliminate the batch_effect in DESeq2
#' 
#' 
#' ### edgeR_quasi_likelihood
#' MAE |>
#'   EMP_assay_extract(experiment = 'geno_ec') |>
#'   EMP_diff_analysis(method='edgeR_quasi_likelihood',
#'                     .formula = ~0+Group,estimate_group = c('Group_B','Group_A')) ## Set the comparison order.
EMP_diff_analysis <- function(obj,experiment,.formula,
                              method = 'wilcox.test',p.adjust='fdr',estimate_group=NULL,
                              use_cached = TRUE,action='add',group_level=NULL,
                              core=NULL,...){
  call <- match.call()
  if (inherits(obj,"MultiAssayExperiment")) {
    EMPT <- .as.EMPT(obj,
                     experiment = experiment)
    .get.method.EMPT(EMPT) <- method
  }else if(inherits(obj,'EMPT')) {
    EMPT <- obj
    .get.method.EMPT(EMPT) <- method
  }
  if (use_cached == FALSE) {
    memoise::forget(.EMP_diff_analysis_m) %>% invisible()
    memoise::forget(.EMP_diff_analysis_tidybulk_m) %>% invisible()
  }

  # avoid typo for users
  if (method == 'deseq2') {
    method <- 'DESeq2'
  }

  switch(method,
         "edgeR_quasi_likelihood" = {EMPT <- .EMP_diff_analysis_tidybulk_m(EMPT = EMPT,method=method,p.adjust=p.adjust,group_level=group_level,.formula=.formula,...)},
         "edgeR_likelihood_ratio"  = {EMPT <- .EMP_diff_analysis_tidybulk_m(EMPT = EMPT,method=method,p.adjust=p.adjust,group_level=group_level,.formula=.formula,...)},
         "edger_robust_likelihood_ratio"  = {EMPT <- .EMP_diff_analysis_tidybulk_m(EMPT = EMPT,method=method,p.adjust=p.adjust,group_level=group_level,.formula=.formula,...)},
         "DESeq2"  = {EMPT <- .EMP_diff_analysis_tidybulk_m(EMPT = EMPT,method=method,p.adjust=p.adjust,group_level=group_level,.formula=.formula,...)},
         "limma_voom"  = {EMPT <- .EMP_diff_analysis_tidybulk_m(EMPT = EMPT,method=method,p.adjust=p.adjust,group_level=group_level,.formula=.formula,...)},
         "limma_voom_sample_weights"  = {EMPT <- .EMP_diff_analysis_tidybulk_m(EMPT = EMPT,method=method,p.adjust=p.adjust,group_level=group_level,.formula=.formula,...)},
         {
           EMPT <- .EMP_diff_analysis_m(EMPT = EMPT,method=method,group_level=group_level,estimate_group=estimate_group,core=core,p.adjust=p.adjust,...)
         }
  )
  .get.history.EMPT(EMPT) <- call
  class(EMPT) <- 'EMP_diff_analysis'
  if (action == 'add') {
    return(EMPT)
  }else if(action == 'get') {
    return(.get.result.EMPT(EMPT))
  }else{
    warning('action should be one of add or get!')
  }
}


.EMP_diff_analysis <- function(EMPT,method,
                               estimate_group=NULL,feature_name=NULL,
                               p.adjust='fdr',group_level=NULL,core=NULL,...){
  Estimate_group <- primary <- pvalue <- feature <- sign_group <- vs <- NULL
  message_info <- list()
  estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)

  mapping <- .get.mapping.EMPT(EMPT) %>% dplyr::select(primary,!!estimate_group)

  ## check the missing value in the group label
  if(any(is.na(mapping[[estimate_group]]))) {
    stop('Column ',estimate_group,' has beed deteced missing value, please check and filter them!')
  }

  assay_data <- .get.assay.EMPT(EMPT) %>%
              dplyr::left_join(mapping,by = 'primary') %>%
              dplyr::select(primary,!!estimate_group,everything())
  .get.assay_name.EMPT(EMPT) -> assay_name

  if (is.null(feature_name)) {
    feature_name <- assay_data %>% dplyr::select(where(is.numeric)) %>% colnames()
  }

  if (is.null(group_level)) {
    subgroup <-assay_data[[estimate_group]] %>% unique()
  }else {
    subgroup <- group_level %>% unique()
    assay_data %<>% dplyr::filter(!!dplyr::sym(estimate_group) %in% subgroup)
  }

  
  # When feature num is not many, core = 1 will be more efficient
  if (is.null(core) & length(feature_name) < 500) {
    core <- 1
  }

  diff_result <- .multi_compare(fun=method,data=assay_data,
                feature=feature_name,factorNames=estimate_group,
                subgroup=subgroup,core=core,...) %>% suppressMessages()


  diff_data <-  .get_diff_df(diff_result,feature_name,estimate_group)

  # fold change only support for the relative and counts data
  sign_fold_temp <- .get_sign_fold(assay_data,subgroup,assay_name,group_level,estimate_group)

  diff_result_brief <-  diff_data %>%
      dplyr::mutate(!!p.adjust := stats::p.adjust(pvalue,method = p.adjust)) %>%
      dplyr::full_join(sign_fold_temp, by='feature') %>%
      dplyr::select(feature,Estimate_group,pvalue,!!p.adjust,sign_group,everything())



  #EMPT@deposit[['raw']] <- diff_result
  EMPT@deposit[['diff_analysis_result']] <- diff_result_brief
  .get.method.EMPT(EMPT) <- method

  ## only counts and relative data in two groups provide vs_infomation
  check_assay <- .get.assay_name.EMPT(EMPT) %in% c('counts','relative','integer','coldata')

  check_group <- length(subgroup) == 2

  if(!check_assay | !check_group ) {
    message_info %<>% append(paste0('Only counts, relative and integer data in two groups provide vs_infomation, fold_change and log2FC!'))
  }else{
    .get.estimate_group_info.EMPT(EMPT) <- diff_result_brief %>% dplyr::pull(vs) %>% unique()
  }

  .get.estimate_group.EMPT(EMPT) <- estimate_group
  .get.message_info.EMPT(EMPT) <- message_info
  .get.algorithm.EMPT(EMPT) <- 'diff_analysis'
  .get.info.EMPT(EMPT) <- 'EMP_diff_analysis'
  EMPT
}

.EMP_diff_analysis_m <- memoise::memoise(.EMP_diff_analysis,cache = cachem::cache_mem(max_size = 2048 * 1024^2))

.get_diff_df <- function(data,feature_name,estimate_group) {
  feature <- Estimate_group <- NULL
  dfs <- lapply(data, function(x) {
    temp <- data.frame(pvalue=NA,method=NA)
    temp$pvalue <- x$p.value
    temp$method <- x$method
    return(temp)
  })
  df <- do.call(rbind, dfs) %>% tibble::as_tibble()
  df <- df %>% dplyr::mutate(feature = feature_name,Estimate_group = estimate_group) %>%
    dplyr::select(feature,Estimate_group,everything())
  return(df)
}

## This code '.multi_compare' is modified from package 'MicrobiotaProcess' for diff analysis
## Here add the parallel function

#' @importFrom spsUtil quiet
#' @importFrom snowfall sfLapply
#' @importFrom snowfall sfStop
#' @importFrom rlang new_formula
#' @importFrom parallel detectCores
#' @noRd
.multi_compare <- function(fun,
                            data,
                            feature,
                            factorNames,
                            subgroup=NULL,core,...){
  if (!is.null(subgroup)){
    data <- data[data[[factorNames]] %in% subgroup, ,drop=FALSE]
    data[[factorNames]] <- factor(data[[factorNames]], levels=subgroup)
  }
  myfun <- function(x){
    tmpformula <- rlang::new_formula(as.name(x), as.name(factorNames))
    suppressWarnings(do.call(fun,list(tmpformula,data=data,...)))
  }
  ## set the core
  if(is.null(core)){
    spsUtil::quiet(snowfall::sfInit(parallel = TRUE, cpus = parallel::detectCores() - 1),print_cat = TRUE, message = FALSE, warning = FALSE)
  }else {
    available_core <- parallel::detectCores() - 1
    if(core>0 & core <= available_core){
      spsUtil::quiet(snowfall::sfInit(parallel = TRUE, cpus = core),print_cat = TRUE, message = FALSE, warning = FALSE)
    }else{
      warning("The parameter core number is wrong,now parallel execution on ",available_core," CPUs.")
      spsUtil::quiet(snowfall::sfInit(parallel = TRUE, cpus = parallel::detectCores() - 1),print_cat = TRUE, message = FALSE, warning = FALSE)
    }
  }
  ## pass the data from envrionment
  snowfall::sfExport('factorNames','data','fun')
  snowfall::sfExport('myfun')
  result <- snowfall::sfLapply(feature,myfun)
  snowfall::sfStop()
  return(result)
}

#' @importFrom dplyr across
.get_sign_fold <- function(assay_data,subgroup,assay_name,group_level=NULL,estimate_group) {
  feature <- abundance <- vs <- `.` <- NULL
  means <- assay_data %>%
    dplyr::group_by(!!sym(estimate_group)) %>%
    dplyr::summarise(across(where(is.numeric), ~mean(., na.rm = TRUE)))

  means_long <- means %>%
    tidyr::pivot_longer(cols = -!!sym(estimate_group),
                        names_to = 'feature',
                        values_to = 'abundance'
    )


  means_long %>%
    dplyr::group_by(feature) %>%
    dplyr::top_n(1, abundance) %>%
    dplyr::rename(sign_group = !!estimate_group) %>%
    dplyr::select(-abundance) %>%
    dplyr::distinct(feature, .keep_all = TRUE) -> deposit  ## filter the same data in some extreme same feature


  if (length(subgroup) == 2 & assay_name %in% c('relative','counts','integer','coldata')) {
    if (is.null(group_level)) {
      group_level <- means[[estimate_group]] %>% as.factor() %>% rev()
    }else {
      group_level <- group_level
    }
    # if 0 value, select the 0.0001 value in other group
    means %<>% dplyr::mutate(across(everything(), ~ ifelse(. == 0, ifelse(is.na(dplyr::lag(.)), dplyr::lead(.), dplyr::lag(.)) * 0.001, .)))
    df1 <-means %>% dplyr::filter(!!sym(estimate_group) == group_level[1])
    df2 <-means %>% dplyr::filter(!!sym(estimate_group) == group_level[2])
    fold_change <- df1[,-1]/df2[,-1]
    fold_change_re <- fold_change %>% dplyr::mutate(vs = paste0(group_level[1],' vs ',group_level[2])) %>%
      tidyr::pivot_longer(cols = -vs,
                          names_to = 'feature',
                          values_to = 'fold_change') %>%
      dplyr::mutate(log2FC = log2(fold_change))


    deposit %<>% dplyr::full_join(.,fold_change_re,by = 'feature')
  }

  return(deposit)
}


