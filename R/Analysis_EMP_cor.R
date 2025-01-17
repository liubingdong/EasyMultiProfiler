#' @importFrom dplyr arrange
#' @importFrom tibble column_to_rownames
#' @noRd
.EMP_cor_analysis <- function(EMP,select=NULL,method='spearman') {
  #call <- match.call()

  primary <- var1 <- NULL
  if (is.null(select)) {
    data1 <- EMP@ExperimentList[[1]] %>% EMP_assay_extract(action='get') %>%
      dplyr::arrange(primary) %>% tibble::column_to_rownames('primary') %>% suppressMessages()
    if (length(EMP@ExperimentList) == 1) {
      data2 <- EMP@ExperimentList[[1]] %>% EMP_assay_extract(action='get') %>%
        dplyr::arrange(primary) %>% tibble::column_to_rownames('primary') %>% suppressMessages()
      cor_info <- names(EMP@ExperimentList)
    }else{
      data2 <- EMP@ExperimentList[[2]] %>% EMP_assay_extract(action='get') %>%
        dplyr::arrange(primary) %>% tibble::column_to_rownames('primary') %>% suppressMessages()
      cor_info <- names(EMP@ExperimentList[1:2])
    }
  }else{
    # check the experiment num for the next step
    select_num <- select %>% unique() %>% length()
    if(select_num >3){
      stop("EMP_cor_analysis only support two experiments!")
    }else if(select_num == 1){
      select <- c(select,select)
    }

    experiment_name <- names(EMP@ExperimentList)
    if (!all(select %in% experiment_name)) {
      stop("Pararmeter select in not in the ExperimentList,please check!")
    }

    cor_info <- select

    data1 <- EMP@ExperimentList[[select[1]]] %>% EMP_assay_extract(action='get') %>%
        dplyr::arrange(primary) %>% tibble::column_to_rownames('primary') %>% suppressMessages()
    data2 <- EMP@ExperimentList[[select[2]]] %>% EMP_assay_extract(action='get') %>%
        dplyr::arrange(primary) %>% tibble::column_to_rownames('primary') %>% suppressMessages()

  }

  # filter the samnples with miss value
  data1 <- na.omit(data1)
  data2 <- na.omit(data2)

  real_sample <- intersect(rownames(data1),rownames(data2))
  data1_sample_num <- rownames(data1) %>% unique %>% length()
  data2_sample_num <- rownames(data2) %>% unique %>% length()

  data1_pre <- data1 |> dplyr::filter(rownames(data1) %in% real_sample)
  data2_pre <- data2 |> dplyr::filter(rownames(data2) %in% real_sample)
  #df.cor.p<-agricolae_correlation(x=data1,y=data2,method = method,...)

  ## check the consistent value
  data1<- data1_pre %>%
    dplyr::select_if(~ dplyr::n_distinct(.) > 1)
  data2 <- data2_pre %>%
    dplyr::select_if(~ dplyr::n_distinct(.) > 1)

  if (ncol(data1_pre) != ncol(data1) | ncol(data2_pre) != ncol(data2)) {
     warning(' The consistent value has beed detected and be removed!')
  }

  df.cor.p <- CorRcpp(x = data1,y = data2,type = method)
  names(df.cor.p) <- c('correlation','pvalue')

  df.cor.p[["correlation"]] <- round(df.cor.p[["correlation"]],2)
  df.cor.p[["pvalue"]] <- round(df.cor.p[["pvalue"]],3)

  df <- df.cor.p$correlation %>%
    as.data.frame() %>%
    tibble::rownames_to_column('var1') %>%
    tidyr::pivot_longer(cols = -var1,names_to = 'var2',values_to = 'coefficient')

  df.p <- df.cor.p$pvalue %>%
    as.data.frame() %>%
    tibble::rownames_to_column('var1') %>%
    tidyr::pivot_longer(cols = -var1,names_to = 'var2',values_to = 'pvalue')
  
  df$pvalue <- df.p$pvalue


  df.cor.p[['cor_info']] <- cor_info
  df.cor.p[["n.obs"]] <- c(data1_sample_num,data2_sample_num)
  df.cor.p[['cor_p']] <- df
  EMP@deposit[['cor_analysis_result']] <- df.cor.p
  #.get.history.EMP(EMP) <- call
  .get.info.EMP(EMP) <-'EMP_cor_analysis'
  .get.method.EMP(EMP) <- method
  class(EMP) <- 'EMP_cor_analysis'
  return(EMP)
}

.EMP_cor_analysis_m <- memoise::memoise(.EMP_cor_analysis,cache = cachem::cache_mem(max_size = 4096 * 1024^2))





.EMP_cor_analysis_multi <- function(EMP,select=NULL,method='spearman',pvalue=0.05,rvalue=0) {
  
  primary <- name <- NULL

  if (is.null(select)) {
    select <- names(EMP@ExperimentList)
  }
  
  df.cor.p <- list()
  cor_info <- select
  n.obs <- c()
  select_num <- length(select)
  experiment_num <- length(EMP@ExperimentList)
  if (select_num > experiment_num) {
    stop("Please check the parameter select is correct or not!")
  }
  
  
  data_list <- list()
  for (i in select) {
    data_list[[i]] <- EMP@ExperimentList[[i]] %>% 
      EMP_assay_extract(action='get') %>%
      dplyr::arrange(primary) %>% 
      dplyr::rename(SampleID = primary) %>% 
      as.data.frame() %>% 
      suppressMessages()
    n.obs <- append(n.obs,data_list[[i]] %>% nrow())
  }
  
  data_list <- check_duplicate_col(data_list)
    
  # contruct the relationship
  total_data <- data_list
  data_length <- length(total_data)
  relationship_stock <- list()
  count_id <- 1
  relationship_stock <- list()
  while (data_length>1) {
    if (count_id == 1 ) {
      temp <- rel_cons(data1 = total_data[[1]],data2 = total_data[[2]], 
                       pvalue = pvalue,rvalue = rvalue,cor_method = method)
    }else{
      keep_idx <- relationship_stock[[count_id-1]]$keep
      temp <- rel_cons(data1 = total_data[[1]][,c('SampleID',keep_idx)], data2 = total_data[[2]], 
                       pvalue = pvalue,rvalue = rvalue,cor_method = method)
    }
    if (length(temp$keep)>0) {
      relationship_stock[[count_id]] <- temp
      total_data <- total_data[-1]
      data_length <- data_length-1  
      count_id <- count_id+1
    }else{
      message('No proper relationship was constructed, please reset the parameters! ')
      break
    }
  }
  
  data_long_raw <- lapply(relationship_stock, function(x) {
    temp_list <- x$relation})
  data_long_raw <- do.call(rbind,data_long_raw)  
  
  # filter the isolated node
  kick_id <- kick_check(data = relationship_stock)
  relationship_stock_filter <- relationship_stock
  if (length(kick_id) != 0) {
    while (length(kick_id >0)) {
      for (i in 1:length(relationship_stock_filter)) {
        idx_temp1 <- relationship_stock_filter[[i]]$relation$target %in% kick_id
        idx_temp2 <- relationship_stock_filter[[i]]$keep %in% kick_id
        relationship_stock_filter[[i]]$relation <- relationship_stock_filter[[i]]$relation[!idx_temp1,]
        relationship_stock_filter[[i]]$keep <- relationship_stock_filter[[i]]$keep[!idx_temp2]
      }
      kick_id=kick_check(data = relationship_stock_filter)
    }
  }
  
  # combine the data
  data_long <- c()
  count_id <- length(relationship_stock_filter)
  relationship_combind <- relationship_stock_filter
  if (count_id == 1 ) {
    data_long <- relationship_stock_filter[[1]]$relation
  }else{
    while ( count_id >=2 ) {
      temp <- rbind(relationship_combind[[1]]$relation,relationship_combind[[2]]$relation)
      count_id=count_id-1
      relationship_combind=relationship_combind[-1]
      data_long <- rbind(data_long,temp)
    }
  }

  feature_info <- do.call(rbind, lapply(names(data_list), function(df_name) {
    data.frame(
      name = names(data_list[[df_name]]),
      group = df_name,
      stringsAsFactors = FALSE
    ) %>% dplyr::filter(name != 'SampleID')
  }))
  
  feature_info$group <- factor(feature_info$group,levels = select) # necessary
  
  df.cor.p[['cor_info']] <- cor_info
  df.cor.p[["n.obs"]] <- n.obs
  df.cor.p[['cor_p']] <- data_long_raw
  df.cor.p[['cor_p_filter']]<- data_long %>% dplyr::distinct()
  df.cor.p[['feature_info']]<- feature_info
  EMP@deposit[['cor_analysis_result']] <- df.cor.p
  .get.info.EMP(EMP) <-'EMP_cor_analysis'
  .get.method.EMP(EMP) <- method
  class(EMP) <- 'EMP_cor_analysis'
  return(EMP)
}

.EMP_cor_analysis_multi_m <- memoise::memoise(.EMP_cor_analysis_multi,cache = cachem::cache_mem(max_size = 4096 * 1024^2))

check_duplicate_col <- function(list_of_dfs) {
  # collect the colname
  all_column_names <- unlist(lapply(list_of_dfs, names))
  
  # check duplicate
  duplicated_column_names <- all_column_names[duplicated(all_column_names)]
  
  # 如果有重复列名，则为每个数据框的列名添加数据框的名称
  # if true , add experiment name to the colname
  if (length(duplicated_column_names) > 0) {
    for (df_name in names(list_of_dfs)) {
      df <- list_of_dfs[[df_name]]
      duplicated_cols <- names(df) %in% duplicated_column_names
      if (any(duplicated_cols)) {
        new_col_names <- paste(df_name, names(df), sep = "_")
        names(df)[duplicated_cols] <- new_col_names[duplicated_cols]
        colnames(df)[1] <- 'SampleID'
        list_of_dfs[[df_name]] <- df
        
      }
    }
  }
  
  return(list_of_dfs)
  
}

rel_cons <- function(data1,data2,pvalue=0.05,rvalue=0,cor_method='spearman'){
  
  SampleID <- value <- NULL

  deposit=list()


  # filter the samnples with miss value
  data1 <- na.omit(data1)
  data2 <- na.omit(data2)
  
  real_sample <- intersect(data1[['SampleID']],data2[['SampleID']])
  
  data1 <- data1 %>% dplyr::filter(SampleID %in% real_sample)  %>% 
    dplyr::arrange(SampleID) %>%
    tibble::column_to_rownames('SampleID')
  
  data2 <- data2 %>% dplyr::filter(SampleID %in% real_sample)  %>% 
    dplyr::arrange(SampleID) %>%
    tibble::column_to_rownames('SampleID')
  
 
  #data.corr <- agricolae_correlation(data1, data2,method = cor_method)
  data.corr <- CorRcpp(x = data1,y = data2,type = cor_method)
  names(data.corr) <- c('correlation','pvalue')

  data.corr[["correlation"]] <- round(data.corr[["correlation"]],2)
  data.corr[["pvalue"]] <- round(data.corr[["pvalue"]],2)


  occor.r <- data.corr$correlation
  occor.p <- data.corr$pvalue
  occor.r[occor.p>pvalue|abs(occor.r)<rvalue] = 0 
  
  data_long <-occor.r %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="source") %>% 
    tidyr::gather(key="target", value="value", -1) %>%
    dplyr::filter(value != 0)
  
  data_long_p <- occor.p %>%
    as.data.frame() %>% 
    tibble::rownames_to_column(var="source") %>% 
    tidyr::gather(key="target", value="pvalue", -1)
  
  data_long <- dplyr::left_join(data_long,data_long_p,by=c('source','target'))
  
  keep=unique(data_long$target)
  deposit$relation=data_long
  deposit$keep=keep
  
  return(deposit)
}

kick_check <- function(data){
  layer_length <- length(data)
  kick_id <- c()
  while( layer_length >= 2 ){
    kick_temp <- data[[1]]$keep[!data[[1]]$keep %in% unique(data[[2]]$relation$source)]
    kick_id <- append(kick_id,kick_temp)
    data <- data[-1]
    layer_length <- layer_length-1
  } 
  return(kick_id)
}

#' Correlation analysis
#'
#' @param EMP Object in EMP format.
#' @param select A character string. The experiment name in the EMP object.
#' @param method A character string. Methods include pearson (default), spearman.
#' @param rvalue A number. Set relation value forthreshold for correlation test (default:0). Only activated in sankey cor-analysis.
#' @param pvalue A number. Set pvalue forthreshold for correlation test (default:0.05). Only activated in sankey cor-analysis.
#' @param force_sankey A boolean. Whether force the cor-analysis for the sankey plot or not.
#' @param action A character string.A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
#' @rdname EMP_cor_analysis
#' @return EMPT object
#' @export
#'
#' @examples
#' data(MAE)
#' k1 <- MAE |>
#'   EMP_assay_extract('taxonomy') |>
#'   EMP_collapse(estimate_group = 'Genus',collapse_by = 'row') |>
#'   EMP_diff_analysis(method='DESeq2', .formula = ~Group) |>
#'   EMP_filter(feature_condition = pvalue<0.05)
#' 
#' k2 <- MAE |>
#'   EMP_collapse(experiment = 'untarget_metabol',na_string=c('NA','null','','-'),
#'                estimate_group = 'MS2kegg',method = 'sum',collapse_by = 'row') |>
#'   EMP_diff_analysis(method='DESeq2', .formula = ~Group) |>
#'   EMP_filter(feature_condition = pvalue<0.05 & abs(fold_change) > 1.5)
#'
#' k3 <- MAE |>
#'   EMP_assay_extract('geno_ec') |>
#'   EMP_diff_analysis(method='DESeq2', .formula = ~Group) |>
#'   EMP_filter(feature_condition = pvalue<0.05 & abs(fold_change) > 2)
#' 
#' #For two experinemnt 
#' (k1 + k2) |> EMP_cor_analysis(method = 'spearman') |>
#'   EMP_heatmap_plot() ## Visualization
#'
#' #For more experinemnt 
#' (k1 + k3 + k2) |> EMP_cor_analysis() |>
#'   EMP_sankey_plot()
EMP_cor_analysis <- function(EMP,select=NULL,method='spearman',action='add',rvalue=0,pvalue=0.05,
                             use_cached=TRUE,force_sankey=FALSE) {

  experiment_num <- NULL

  call <- match.call()
  
  if (!is(EMP,"EMP")) {
    stop("Please input the EMP format!")
  }

  if (use_cached == FALSE) {
    memoise::forget(.EMP_cor_analysis_m) %>% invisible()
    memoise::forget(.EMP_cor_analysis_multi_m) %>% invisible()
  }
  
  experiment_num <- length(EMP@ExperimentList)

  if (force_sankey == TRUE) {
    experiment_num <- 3
  }

  if (experiment_num <= 2) {
    result <- .EMP_cor_analysis_m(EMP=EMP,select=select,method=method)
  }else {
    result <- .EMP_cor_analysis_multi_m(EMP=EMP,select=select,method=method,rvalue=rvalue,pvalue=pvalue)
  }


 .get.history.EMP(result) <- call
  if (action=='get') {
    return(result@deposit[['cor_analysis_result']])
  }else if(action=='add') {
    return(result)
  }else{
    warning('action should be one of add or get!')
  }

}



