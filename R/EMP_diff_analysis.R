#' Title
#'
#' @param EMPT wait_for_add
#' @param group_level wait_for_add
#' @param design wait_for_add
#' @param ... wait_for_add
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 DESeq
#' @importFrom DESeq2 results
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#' @importFrom tibble as_tibble
#' @importFrom stringr str_remove
#' @importFrom SummarizedExperiment colData
#' @noRd
.EMP_diff_analysis_deseq2_deprecated <- function(EMPT,group_level=NULL,design,...) {
  log2FoldChange <- log2FC <- feature <- Estimate_group <- pvalue <- padj <- sign_group <- NULL
  method <- vs <- fold_change <- batch_effect <- NULL
  # confirm group and batch info
  Group_info <- as.list(design)[[2]]
  estimate_group <- as.list(Group_info)[[length(Group_info)]] %>% as.character()
  if (length(as.list(Group_info)) ==3) {
    batch_info <- as.list(Group_info)[[2]] %>% as.character()
  }else{
    batch_info <- NULL
  }

  # prepare data for deseq2
  sub_data <- .get.result.EMPT(EMPT) %>%
    tibble::column_to_rownames('primary') %>%
    round(digits = 0) %>%
    t()
  colData <- .get.mapping.EMPT(EMPT) %>%
    tibble::column_to_rownames('primary')  %>%
    dplyr::select(!!estimate_group,!!batch_info)

  if (is.null(group_level)) {
    group_level <- colData[[estimate_group]] %>% unique()
  }

  colData[[estimate_group]] <- factor(colData[[estimate_group]],levels = rev(group_level))


  cds <- DESeq2::DESeqDataSetFromMatrix(sub_data,colData,design) %>%
    suppressMessages() %>%
    suppressWarnings()

  #fitType <- match.arg(fitType, choices=c("parametric","local","mean","glmGamPoi"))
  #test <- match.arg(test, choices=c("Wald","LRT"))
  dds <- DESeq2::DESeq(cds,...)

  res <- DESeq2::results(dds)


  # obetain vs info
  vs_info <- as.data.frame(res@elementMetadata)[2,2]


  res_df <- as.data.frame(res)  %>%
    tibble::rownames_to_column('feature') %>%
    dplyr::rename(log2FC = log2FoldChange) %>%
    dplyr::mutate(vs = vs_info) %>%
    dplyr::mutate(fold_change = 2^log2FC) %>%
    dplyr::mutate(
      sign_group = dplyr::case_when(
        log2FC > 0 ~ group_level[1],
        log2FC < 0 ~ group_level[2]
      )
    ) %>%
    dplyr::mutate(Estimate_group=estimate_group) %>%
    dplyr::mutate(method='deseq2') %>%
    dplyr::select(feature,Estimate_group,pvalue,padj,sign_group,method,vs,fold_change,log2FC,everything()) %>%
    tibble::as_tibble()

  if (!is.null(batch_info)) {
    res_df %<>%
      dplyr::mutate(batch_effect=batch_info) %>%
      dplyr::select(feature,Estimate_group,batch_effect,everything())
  }


  #EMPT@deposit[['raw']] <- res
  EMPT@deposit[['diff_analysis_result']] <- res_df
  .get.estimate_group.EMPT(EMPT) <- group_level
  .get.algorithm.EMPT(EMPT) <- 'diff_analysis'
  .get.info.EMPT(EMPT) <- 'EMP_diff_analysis'
  EMPT
}


.EMP_diff_analysis_deseq2_m_deprecated <- memoise::memoise(.EMP_diff_analysis_deseq2_deprecated)


#' @importFrom rlang `:=`
#' @importFrom dplyr last_col
#' @importFrom dplyr matches
.EMP_diff_analysis_tidybulk <- function(EMPT,method,.formula,p.adjust='fdr',group_level=NULL,...) {
  pvalue <- feature <- Estimate_group <- sign_group <- vs <- log2FC <- Estimate_group <- fold_change <- NULL
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

  result <- EMPT %>% tidybulk::tidybulk() %>%
    tidybulk::test_differential_abundance(
      .formula,
      method = method,action='get',contrasts=group_level,...) %>%
    suppressWarnings()  %>% suppressMessages()


  if (!is.null(group_level)) {
    design_raw_info <- result %>%
      dplyr::select(last_col()) %>% colnames() %>%
      strsplit(.,split = '___') %>% unlist()

    design_fix <- paste0('___',design_raw_info[2])

    result %<>%
      dplyr::rename_with(~stringr::str_remove(., design_fix))

    design_info <- design_raw_info[2]
    design_info %<>% gsub(estimate_group, "",.) %>% trimws() %>%
      gsub("-", " vs ",.)
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

.EMP_diff_analysis_tidybulk_m <- memoise::memoise(.EMP_diff_analysis_tidybulk)



#' Title
#'
#' @param x wait_for_add
#' @param experiment wait_for_add
#' @param method wait_for_add
#' @param estimate_group wait_for_add
#' @param use_cached wait_for_add
#' @param action wait_for_add
#' @param group_level wait_for_add
#' @param core wait_for_add
#' @param ... wait_for_add
#' @importFrom memoise forget
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_diff_analysis <- function(x,experiment,
                              method = 'wilcox.test',estimate_group=NULL,
                              use_cached = T,action='add',group_level=NULL,
                              core=NULL,...){
  call <- match.call()
  if (inherits(x,"MultiAssayExperiment")) {
    EMPT <- .as.EMPT(x,
                     experiment = experiment)
    .get.method.EMPT(EMPT) <- method
  }else if(inherits(x,'EMPT')) {
    EMPT <-x
    .get.method.EMPT(EMPT) <- method
  }
  if (use_cached == F) {
    memoise::forget(.EMP_diff_analysis_m) %>% invisible()
    memoise::forget(.EMP_diff_analysis_tidybulk_m) %>% invisible()
  }
  switch(method,
         "edgeR_quasi_likelihood" = {EMPT <- .EMP_diff_analysis_tidybulk_m(EMPT = EMPT,method=method,group_level=group_level,...)},
         "edgeR_likelihood_ratio"  = {EMPT <- .EMP_diff_analysis_tidybulk_m(EMPT = EMPT,method=method,group_level=group_level,...)},
         "edger_robust_likelihood_ratio"  = {EMPT <- .EMP_diff_analysis_tidybulk_m(EMPT = EMPT,method=method,group_level=group_level,...)},
         "DESeq2"  = {EMPT <- .EMP_diff_analysis_tidybulk_m(EMPT = EMPT,method=method,group_level=group_level,...)},
         "limma_voom"  = {EMPT <- .EMP_diff_analysis_tidybulk_m(EMPT = EMPT,method=method,...)},
         "limma_voom_sample_weights"  = {EMPT <- .EMP_diff_analysis_tidybulk_m(EMPT = EMPT,method=method,group_level=group_level,...)},
         {
           EMPT <- .EMP_diff_analysis_m(EMPT = EMPT,method=method,group_level=group_level,estimate_group=estimate_group,core=core,...)
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


.EMP_diff_analysis <- function(EMPT,experiment,assay_name,method,
                               estimate_group=NULL,feature_name=NULL,
                               p.adjust='fdr',group_level=NULL,core=NULL,...){
  primary <- NULL
  message_info <- list()
  estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)

  mapping <- .get.mapping.EMPT(EMPT) %>% dplyr::select(primary,!!estimate_group)

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


  diff_result <- .multi_compare(fun=method,data=assay_data,
                feature=feature_name,factorNames=estimate_group,
                subgroup=subgroup,core=core,...) %>% suppressMessages()


  diff_data <-  .get_diff_df(diff_result)

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

.EMP_diff_analysis_m <- memoise::memoise(.EMP_diff_analysis)

#' @importFrom stringr str_split
.get_diff_df <- function(data) {
  dfs <- lapply(data, function(x) {
    temp_info <- x[["data.name"]] %>% stringr::str_split(pattern = 'by') %>% unlist() %>% trimws()
    temp <- data.frame(feature = temp_info[1],Estimate_group=temp_info[2])
    temp$pvalue <- x$p.value
    temp$method <- x$method
    return(temp)
  })

  df <- do.call(rbind, dfs) %>% tibble::as_tibble()
  return(df)
}

## This code '.multi_compare' is modified from package 'MicrobiotaProcess' for diff analysis
## Here add the parallel function
#' Title
#'
#' @param fun wait_for_add
#' @param data wait_for_add
#' @param feature wait_for_add
#' @param factorNames wait_for_add
#' @param subgroup wait_for_add
#' @param core wait_for_add
#' @param ... wait_for_add
#' @importFrom spsUtil quiet
#' @importFrom snowfall sfLapply
#' @importFrom snowfall sfStop
#' @importFrom rlang new_formula
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
    spsUtil::quiet(snowfall::sfInit(parallel = TRUE, cpus = chectCores() - 1),print_cat = TRUE, message = FALSE, warning = FALSE)
  }else {
    available_core <- chectCores() - 1
    if(core>0 & core <= available_core){
      spsUtil::quiet(snowfall::sfInit(parallel = TRUE, cpus = core),print_cat = TRUE, message = FALSE, warning = FALSE)
    }else{
      warning("The parameter core number is wrong,now parallel execution on ",available_core," CPUs.")
      spsUtil::quiet(snowfall::sfInit(parallel = TRUE, cpus = chectCores() - 1),print_cat = TRUE, message = FALSE, warning = FALSE)
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
    dplyr::select(-abundance) -> deposit


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

## This code 'chectCores' is from the function 'detectCores' in the package 'parallel' for core detect
#' Title
#'
#' @param all.tests wait_for_add
#' @param logical wait_for_add
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
chectCores <- function (all.tests = FALSE, logical = TRUE){
  systems <- list(linux = "grep \"^processor\" /proc/cpuinfo 2>/dev/null | wc -l",
                  darwin = if (logical) "/usr/sbin/sysctl -n hw.logicalcpu 2>/dev/null" else "/usr/sbin/sysctl -n hw.physicalcpu 2>/dev/null",
                  solaris = if (logical) "/usr/sbin/psrinfo -v | grep 'Status of.*processor' | wc -l" else "/bin/kstat -p -m cpu_info | grep :core_id | cut -f2 | uniq | wc -l",
                  freebsd = "/sbin/sysctl -n hw.ncpu 2>/dev/null", openbsd = "/sbin/sysctl -n hw.ncpuonline 2>/dev/null")
  nm <- names(systems)
  m <- pmatch(nm, R.version$os)
  m <- nm[!is.na(m)]
  if (length(m)) {
    cmd <- systems[[m]]
    if (!is.null(a <- tryCatch(suppressWarnings(system(cmd,
                                                       TRUE)), error = function(e) NULL))) {
      a <- gsub("^ +", "", a[1])
      if (grepl("^[1-9]", a))
        return(as.integer(a))
    }
  }
  if (all.tests) {
    for (i in seq(systems)) for (cmd in systems[i]) {
      if (is.null(a <- tryCatch(suppressWarnings(system(cmd,
                                                        TRUE)), error = function(e) NULL)))
        next
      a <- gsub("^ +", "", a[1])
      if (grepl("^[1-9]", a))
        return(as.integer(a))
    }
  }
  NA_integer_
}
