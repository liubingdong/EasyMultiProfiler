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

  try( melt_EMPT[[estimate_group]] <- droplevels(melt_EMPT[[estimate_group]]),silent = TRUE)

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
    
    if (!is.factor(melt_EMPT[[estimate_group]])) {
      melt_EMPT[[estimate_group]] <- factor(melt_EMPT[[estimate_group]])
    }

    design_raw_info <- melt_EMPT[[estimate_group]] %>% unique() %>% sort()

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

.EMP_diff_analysis_tidybulk_m <- memoise::memoise(.EMP_diff_analysis_tidybulk,cache = cachem::cache_mem(max_size = 4096 * 1024^2))



#' Differential expression or abundance analysis
#'
#' @param obj Object in EMPT or MultiAssayExperiment format.
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object.
#' @param method A character string. Methods include t.test, wilcox.test, kruskal.test, oneway.test, edgeR_quasi_likelihood, edgeR_likelihood_ratio, edger_robust_likelihood_ratio, DESeq2, limma_voom, limma_voom_sample_weights
#' @param p.adjust A character string. Adjust P-values for Multiple Comparisons inluding fdr, holm, hochberg, hommel, bonferroni, BH, BY. (default:fdr)
#' @param .formula A formula representing the desired linear model. If there is more than one factor, they should be in the order factor of interest + additional factors.
#' @param estimate_group A character string. Select the group name in the coldata to be calculated.
#' @param paired_group  A character string. Variable name corresponding to paired primary or sample.
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
#' @param action A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param group_level A series of character strings. Determine the comparison order of groups. when group_level activted,.formula should be ~0+[your interested group] for edger and limma.
#' @param core A number. Select the core numbers in the parallel compute to speed up the result.When no core is set and the number of features exceeds 2000, the default will use all CPU cores minus one.
#' @param ... Further parameters passed to the function \code{\link[tidybulk]{test_differential_abundance}} or \code{\link[stats]{t.test}}, \code{\link[stats]{wilcox.test}}, \code{\link[stats]{kruskal.test}},\code{\link[stats]{oneway.test}}.
#' @importFrom memoise forget
#' @return EMPT object
#' @export
#'
#' @section Detaild about p.adjust:
#' When using common transcriptome differential analysis methods, such as DESeq2, it’s important to note that the method internally filters genes before adjusting p-values. 
#' As a result, the adjusted p-values may differ from those obtained using unfiltered data. 
#' To ensure consistency, you can first apply gene filtering using the EMP_identify_assay function with the edgeR method. 
#' This approach ensures that the adjusted p-values remain consistent.
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
#' ## edgeR_quasi_likelihood
#' MAE |>
#'   EMP_assay_extract(experiment = 'geno_ec') |>
#'   EMP_diff_analysis(method='edgeR_quasi_likelihood',
#'                     .formula = ~0+Group,group_level = c('Group_B','Group_A')) ## Set the comparison order.
#'
#' ## Paired test
#' MAE |> 
#'   EMP_assay_extract('host_gene',
#'                     pattern = 'A1BG',pattern_ref = 'feature') |> 
#'   EMP_filter(sub_group %in% c('A','B')) |>
#'   EMP_diff_analysis(method = 't.test',
#'                     paired_group='patient', # Set the paired group
#'                     estimate_group = 'sub_group')
EMP_diff_analysis <- function(obj,experiment,.formula,
                              method = 'wilcox.test',p.adjust='fdr',estimate_group=NULL,paired_group=NULL,
                              use_cached = TRUE,action='add',group_level=NULL,
                              core='auto',...){
  call <- match.call()
  if (is(obj,"MultiAssayExperiment")) {
    EMPT <- .as.EMPT(obj,
                     experiment = experiment)
    .get.method.EMPT(EMPT) <- method
  }else if(is(obj,'EMPT')) {
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
           EMPT <- .EMP_diff_analysis_m(EMPT = EMPT,method=method,group_level=group_level,estimate_group=estimate_group,core=core,p.adjust=p.adjust,paired_group=paired_group,...)
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
                               estimate_group=NULL,feature_name=NULL,paired_group=NULL,
                               p.adjust='fdr',group_level=NULL,core=NULL,...){
  Estimate_group <- primary <- pvalue <- feature <- sign_group <- vs <- NULL
  message_info <- list()
  estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)

  mapping <- .get.mapping.EMPT(EMPT) %>% dplyr::select(primary,dplyr::any_of(c(!!estimate_group,!!paired_group)))

  ## check the missing value in the group label
  if(any(is.na(mapping[[estimate_group]]))) {
    stop('Column ',estimate_group,' has beed deteced missing value, please check and filter them!')
  }

  if (is.null(paired_group)) {
      assay_data <- .get.assay.EMPT(EMPT) %>%
              dplyr::left_join(mapping,by = 'primary') %>%
              dplyr::select(primary,!!estimate_group,everything())
      paired <- FALSE
  }else{
      ## check the missing value in the group label
      if(any(is.na(mapping[[paired_group]]))) {
        stop('Column ',paired_group,' has beed deteced missing value, please check and filter them!')
      }
      ## check the paired sample and size
      check_paired_result <- .check_group_consistency(data=mapping,group_col=estimate_group,patient_col=paired_group)
      if(check_paired_result$all_patients_equal == FALSE) {
        stop(check_paired_result$message)
      }
      assay_data <- .get.assay.EMPT(EMPT) %>%
              dplyr::left_join(mapping,by = 'primary') %>%
              dplyr::arrange(!!dplyr::sym(estimate_group),!!dplyr::sym(paired_group)) %>%
              dplyr::select(primary,!!estimate_group,everything())
      paired <- TRUE
  }
  
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
  if (core == 'auto' & length(feature_name) < 2000) {
    core <- 1
  }


  if (paired) {
    diff_result <- .multi_compare_paired(fun=method,data=assay_data,
                                    feature=feature_name,factorNames=estimate_group,
                                    subgroup=subgroup,core=core,paired=TRUE,...)  %>% suppressMessages()
  }else{
    diff_result <- .multi_compare_unpaired(fun=method,data=assay_data,
                                          feature=feature_name,factorNames=estimate_group,
                                          subgroup=subgroup,core=core,...) %>% suppressMessages()
  }

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

.EMP_diff_analysis_m <- memoise::memoise(.EMP_diff_analysis,cache = cachem::cache_mem(max_size = 4096 * 1024^2))

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

#' @importFrom spsUtil quiet
#' @importFrom snowfall sfLapply
#' @importFrom snowfall sfStop
#' @importFrom rlang new_formula
#' @importFrom parallel detectCores
#' @noRd
.multi_compare_paired <- function(fun,
                           data,
                           feature,
                           factorNames,
                           subgroup=NULL,core,paired,...){
  if (!is.null(subgroup)){
    data <- data[data[[factorNames]] %in% subgroup, ,drop=FALSE]
    data[[factorNames]] <- factor(data[[factorNames]], levels=subgroup)
  }
  
  if(length(subgroup) != 2){
    stop('Paired test only support 2 groups!')
  }

  if (core==1) {
    result <- lapply(feature,
                     function(x){
                       var1 <- data |> dplyr::filter(!!rlang::sym(factorNames) == subgroup[1]) |> dplyr::pull(x)
                       var2 <- data |> dplyr::filter(!!rlang::sym(factorNames) == subgroup[2]) |> dplyr::pull(x)
                       if (length(var1) != length(var1)) {
                        stop("Paired test need Paired data!")
                       }
                       suppressWarnings(do.call(fun,list(var1,var2,paired=TRUE,...)))})
  }else if (core== 'auto'){
    myfun <- function(x){
      var1 <- data |> dplyr::filter(!!rlang::sym(factorNames) == subgroup[1]) |> dplyr::pull(x)
      var2 <- data |> dplyr::filter(!!rlang::sym(factorNames) == subgroup[2]) |> dplyr::pull(x)
      if (length(var1) != length(var1)) {
       stop("Paired test requires equal-length pairs for feature!")
      }      
      suppressWarnings(do.call(fun,list(var1,var2,paired=TRUE,...)))
    }
    ## set the core
    spsUtil::quiet(snowfall::sfInit(parallel = TRUE, cpus = parallel::detectCores() - 1),print_cat = TRUE, message = FALSE, warning = FALSE)
    snowfall::sfExport('factorNames','data','fun')
    snowfall::sfExport('myfun')
    result <- snowfall::sfLapply(feature,myfun)
    snowfall::sfStop()
  }else{
    myfun <- function(x){
      var1 <- data |> dplyr::filter(!!rlang::sym(factorNames) == subgroup[1]) |> dplyr::pull(x)
      var2 <- data |> dplyr::filter(!!rlang::sym(factorNames) == subgroup[2]) |> dplyr::pull(x)
      if (length(var1) != length(var1)) {
       stop("Paired test requires equal-length pairs for feature!")
      }         
      suppressWarnings(do.call(fun,list(var1,var2,paired=TRUE,...)))
    }
    ## set the core
    available_core <- parallel::detectCores() - 1
    if(core>0 & core <= available_core){
      spsUtil::quiet(snowfall::sfInit(parallel = TRUE, cpus = core),print_cat = TRUE, message = FALSE, warning = FALSE)
    }else{
      warning("The parameter core number is wrong,now parallel execution on ",available_core," CPUs.")
      spsUtil::quiet(snowfall::sfInit(parallel = TRUE, cpus = parallel::detectCores() - 1),print_cat = TRUE, message = FALSE, warning = FALSE)
    }
    snowfall::sfExport('factorNames','data','fun')
    snowfall::sfExport('myfun')
    result <- snowfall::sfLapply(feature,myfun)
    snowfall::sfStop()
  }
  return(result)
}

.multi_compare_unpaired <- function(fun,
                            data,
                            feature,
                            factorNames,
                            subgroup=NULL,core,...){
   if (!is.null(subgroup)){
     data <- data[data[[factorNames]] %in% subgroup, ,drop=FALSE]
     data[[factorNames]] <- factor(data[[factorNames]], levels=subgroup)
   }
 
   if (core==1) {
     result <- lapply(feature,
            function(x){
              tmpformula <- rlang::new_formula(as.name(x), as.name(factorNames))
              suppressWarnings(do.call(fun,list(tmpformula,data=data,...)))}) 
   }else if (core== 'auto'){
     myfun <- function(x){
       tmpformula <- rlang::new_formula(as.name(x), as.name(factorNames))
       suppressWarnings(do.call(fun,list(tmpformula,data=data,...)))
     }
     ## set the core
     spsUtil::quiet(snowfall::sfInit(parallel = TRUE, cpus = parallel::detectCores() - 1),print_cat = TRUE, message = FALSE, warning = FALSE)
     snowfall::sfExport('factorNames','data','fun')
     snowfall::sfExport('myfun')
     result <- snowfall::sfLapply(feature,myfun)
     snowfall::sfStop()
   }else{
     myfun <- function(x){
       tmpformula <- rlang::new_formula(as.name(x), as.name(factorNames))
       suppressWarnings(do.call(fun,list(tmpformula,data=data,...)))
     }
     ## set the core
     available_core <- parallel::detectCores() - 1
     if(core>0 & core <= available_core){
       spsUtil::quiet(snowfall::sfInit(parallel = TRUE, cpus = core),print_cat = TRUE, message = FALSE, warning = FALSE)
     }else{
       warning("The parameter core number is wrong,now parallel execution on ",available_core," CPUs.")
       spsUtil::quiet(snowfall::sfInit(parallel = TRUE, cpus = parallel::detectCores() - 1),print_cat = TRUE, message = FALSE, warning = FALSE)
     }
     snowfall::sfExport('factorNames','data','fun')
     snowfall::sfExport('myfun')
     result <- snowfall::sfLapply(feature,myfun)
     snowfall::sfStop()
   }
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
      group_level <- means[[estimate_group]] %>% unique()
    }else {
      group_level <- group_level
    }
    # if 0 value, select the 0.0001 value in other group
    means %<>% dplyr::mutate(across(-sym(!!estimate_group), ~ ifelse(. == 0, ifelse(is.na(dplyr::lag(.)), dplyr::lead(.), dplyr::lag(.)) * 0.001, .)))
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

.check_group_consistency <- function(data, group_col, patient_col, ignore_order = TRUE) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install package 'dplyr' to use this function")
  }
  
  # Validate input columns
  if (!group_col %in% names(data)) {
    stop("Group column '", group_col, "' not found in data")
  }
  if (!patient_col %in% names(data)) {
    stop("Patient column '", patient_col, "' not found in data")
  }
  
  # Calculate group statistics
  group_stats <- data %>%
    dplyr::group_by(!!rlang::sym(group_col)) %>%
    dplyr::summarise(
      Patients = if (ignore_order) {
        paste(sort(!!rlang::sym(patient_col)), collapse = ", ")
      } else {
        paste(!!rlang::sym(patient_col), collapse = ", ")
      },
      Count = dplyr::n(),
      .groups = "drop"
    )
  
  # Check consistency
  all_patients_equal <- dplyr::n_distinct(group_stats$Patients) == 1
  all_counts_equal <- dplyr::n_distinct(group_stats$Count) == 1
  
  # Generate result message
  result_message <- if (all_patients_equal) {
    "All groups have identical patient lists and sizes"
  } else if (all_counts_equal) {
    paste0("Group sizes match but paired group differ!")
  } else {
    paste0("Both paired group and group sizes differ!")
  }
  
  # Return comprehensive results
  list(
    group_stats = as.data.frame(group_stats),
    all_patients_equal = all_patients_equal,
    all_counts_equal = all_counts_equal,
    message = result_message
  )
}
