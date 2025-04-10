#' @importFrom purrr reduce
#' @importFrom bootnet estimateNetwork
#' @importFrom qgraph centralityTable
.EMP_network_analysis_EMP <- function(EMP,select=NULL,method='cor',corMethod='spearman',weighted=TRUE,
                              missing='listwise',threshold='sig',...){
  
  measure <- value <- node <- feature <- NULL

  if (is.null(select)) {
    select <- names(EMP@ExperimentList)
  }else{
    experiment_name <- names(EMP@ExperimentList)
    if (!all(select %in% experiment_name)) {
      stop("Pararmeter select in not in the ExperimentList,please check!")
    }
  }
  assay_list <- list()
  rowdata_list <- list()
  for (i in select) {
    assay_list[[i]] <- .get.assay.EMPT(EMP@ExperimentList[[i]])
    rowdata_list[[i]] <- .get.row_info.EMPT(EMP@ExperimentList[[i]]) |>
      dplyr::mutate(.experiment_source = i)
  }
  assay_df <-  reduce(assay_list, dplyr::inner_join, by = 'primary') |>
    tibble::column_to_rownames('primary')

  if (method %in% c('IsingFit','IsingSampler','graphicalVAR','GGMncv')) {
    network <- bootnet::estimateNetwork(assay_df,
                                        default = method,
                                        missing = missing,
                                        weighted=weighted,...) |> 
             spsUtil::quiet(print_cat = TRUE, message = TRUE, warning = TRUE)
  }else if (method %in% c('adalasso','mgm','relimp','TMFG','ggmModSelect','LoGo','piecewiseIsing','SVAR_lavaan')){
    stop("Due to Pararmeter confict, EasyMultiProfiler could not support ",method,', please refer to bootnet package!')
  }else{
      network <- estimateNetwork(assay_df,
                             threshold = threshold,
                             default = method,
                             corMethod = corMethod,
                             missing = missing,
                             weighted=weighted,
                             nonPositiveDefinite='continue',...) |> 
             spsUtil::quiet(print_cat = TRUE, message = TRUE, warning = TRUE)
  }
  
  centrality_info <- centralityTable(network) |>
    tidyr::pivot_wider(names_from = measure, values_from = value) |>
    dplyr::select(node,everything()) |>
    dplyr::mutate(feature = network[["labels"]]) |>
    dplyr::select(feature,everything())


  EMP@deposit[['EMP_network_analysis']][['net']] <- network
  EMP@deposit[['EMP_network_analysis']][['net_feature_info']] <- rowdata_list
  EMP@deposit[['EMP_network_analysis']][['net_centrality']] <- centrality_info

  .get.info.EMP(EMP) <-'EMP_network_analysis2'
  .get.method.EMP(EMP) <- method
  class(EMP) <- 'EMP_network_analysis2'
  return(EMP)
}

.EMP_network_analysis_EMP_m <- memoise::memoise(.EMP_network_analysis_EMP,cache = cachem::cache_mem(max_size = 4096 * 1024^2))

#' @importFrom bootnet estimateNetwork
#' @importFrom qgraph centralityTable
.EMP_network_analysis_EMPT <- function(EMPT,method='cor',corMethod='spearman',weighted=TRUE,
                                  missing='listwise',threshold='sig',coldata_to_assay=NULL,...){

  measure <- value <- node <- feature <- NULL

  if (is.null(coldata_to_assay)) {
    coldata <- NULL
    assay_df <-  .get.assay.EMPT(EMPT) |>
      tibble::column_to_rownames('primary')
    rowdata_df <- .get.row_info.EMPT(EMPT)
    
  }else{
    coldata <- EMPT %>% EMP_coldata_extract() %>%
      dplyr::select(dplyr::all_of(c('primary',!!coldata_to_assay)))
    
    assay_df <-  .get.assay.EMPT(EMPT) |> 
      dplyr::inner_join(coldata,by='primary') |>
      tibble::column_to_rownames('primary')
    
    rowdata_df <- .get.row_info.EMPT(EMPT)
    new_line <-do.call(rbind, lapply(coldata_to_assay, function(x) rep(x, ncol(rowdata_df)))) |> 
      as.data.frame()
    new_line[,-1] <- "coldata"
    colnames(new_line) <- colnames(rowdata_df)
    rowdata_df <- rbind(rowdata_df,new_line)
    
  }

  if (method %in% c('IsingFit','IsingSampler','graphicalVAR','GGMncv')) {
    network <- bootnet::estimateNetwork(assay_df,
                                        default = method,
                                        missing = missing,
                                        weighted=weighted,...) |> 
             spsUtil::quiet(print_cat = TRUE, message = TRUE, warning = TRUE)
  }else if (method %in% c('adalasso','mgm','relimp','TMFG','ggmModSelect','LoGo','piecewiseIsing','SVAR_lavaan')){
    stop("Due to Pararmeter confict, EasyMultiProfiler could not support ",method,', please refer to bootnet package!')
  }else{
      network <- estimateNetwork(assay_df,
                             threshold = threshold,
                             default = method,
                             corMethod = corMethod,
                             missing = missing,
                             weighted=weighted,
                             nonPositiveDefinite='continue',...) |> 
             spsUtil::quiet(print_cat = TRUE, message = TRUE, warning = TRUE)
  }
  
  centrality_info <- qgraph::centralityTable(network) |>
    tidyr::pivot_wider(names_from = measure, values_from = value) |>
    dplyr::select(node,everything()) |>
    dplyr::mutate(feature = network[["labels"]]) |>
    dplyr::select(feature,everything())
  
  EMPT@deposit[['net']] <- network
  EMPT@deposit[['net_feature_info']] <- rowdata_df
  EMPT@deposit[['net_centrality']] <- centrality_info
  
  .get.info.EMPT(EMPT) <-'EMP_network_analysis'
  .get.method.EMPT(EMPT) <- method
  class(EMPT) <- 'EMP_network_analysis'
  return(EMPT)
}

.EMP_network_analysis_EMPT_m <- memoise::memoise(.EMP_network_analysis_EMPT,cache = cachem::cache_mem(max_size = 4096 * 1024^2))


.network_print <- function (obj){
  if (is(obj,'EMP')) {
    network <- obj@deposit[['EMP_network_analysis']][['net']]
  }else if (is(obj,'EMPT')) {
    network <- obj@deposit[['net']]
  }else if (is(obj,'bootnetResult')) {
    network <- obj
  }else{
    stop('Input should be EMP or bootnetResult!')
  }
  sample_num <- network$data |> nrow()
  feature_num <- network$data |> ncol()
  info <- paste0("Network was constructed based on ",sample_num, 
                 ' of samples and ',feature_num, ' of features.')
  EMP_message(info,color = 32,order = 1,show='message')
}

#' Network analysis
#'
#' @param obj EMPT or EMP.
#' @param select A character string. The experiment name in the EMP object.
#' @param method A character string. Methods include cor, EBICglasso,IsingFit,IsingSampler,pcor,huge and graphicalVAR.
#' @param corMethod A character string. Methods include spearman, cor, cov, npn and cor_auto.
#' @param weighted Logical, should the analyzed network be weighted?
#' @param missing How to handle missing data? "pairwise" for pairwise deletion, "listwise" for listwise deletion, "fiml" for full-information maximum likelihood and "stop" to stop with an error.
#' @param threshold This parameter allows input of strings such as 'none', 'sig', 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', or 'fdr' for statistical calculation on edges, with the default set to 'sig' to select only statistically significant edges; alternatively, it can be a numeric value between 0 and 1, which will remove edges with an absolute value greater than this threshold.
#' @param coldata_to_assay A series of character strings. Select the column from coldata to caculate. Only for the EMPT object.
#' @param action A character string.A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
#' @param ... Further parameters passed to the function \code{\link[bootnet]{estimateNetwork}}.
#' @export
#' @rdname EMP_network_analysis
#' @section Detaild about network:
#' More detailed information are available on :
#'
#' Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating psychological networks and their accuracy: A tutorial paper. Behavior Research Methods, 50(1), 195–212. 
#'
#' Here is a brife introduction of document from bootnet package:
#'
#' 1. cor: Correlation network.
#'
#' 2. EBICglasso: Gaussian Markov random field estimation using graphical LASSO and extended Bayesian information criterion to select optimal regularization parameter. Using \code{\link[qgraph]{EBICglasso}} from the qgraph package. Calls bootnet_EBICglasso.
#'
#' 3. IsingSampler: Calls the \code{\link[IsingSampler]{EstimateIsing}} function from the IsingSampler package.
#'
#' 4. pcor: Partial correlation network (non-regularized Gaussian Markov random field), using \code{\link[corpcor]{cor2pcor}} from the corpcor package. Calls bootnet_pcor.
#'
#' 5. huge: Uses EBIC model selection of GGM networks estimated via the glasso algorithm as implemented in the huge package (as opposed to glasso and qgraph packages used in default = "EBICglasso"). Uses nonparanormal transformation in preparing the data and does not use polychoric correlations. Calls bootnet_huge.
#'
#' 6. graphicalVAR.: Estimates a graphical VAR model using the graphicalVAR package. This results in two networks which can be plotted using the 'graph' argument in the plot method. Calls bootnet_graphicalVAR.
#'
#'
#' @return EMP or EMPT object
#'
#' @examples
#'data(MAE)
#'
#' k1 <- MAE |>
#'   EMP_assay_extract('taxonomy') |>
#'   EMP_collapse(estimate_group = 'Genus',collapse_by = 'row') |>
#'   EMP_diff_analysis(method='wilcox.test', estimate_group = 'Group') |>
#'   EMP_filter(feature_condition = pvalue<0.05)
#' 
#' k2 <- MAE |>
#'   EMP_collapse(experiment = 'untarget_metabol',na_string=c('NA','null','','-'),
#'                estimate_group = 'MS2kegg',method = 'sum',collapse_by = 'row') |>
#'   EMP_diff_analysis(method='DESeq2', .formula = ~Group) |>
#'   EMP_filter(feature_condition = pvalue<0.05)
#' 
#' 
#' #For single experiment
#' k1  |>
#'   EMP_network_analysis() |>
#'   EMP_filter(feature_condition = top_detect(Closeness,1))  # filter the top node according to the network analysis
#'
#' #For muti experiment
#' (k1 + k2 ) |> 
#'   EMP_network_analysis() |> 
#'   EMP_network_plot(show = 'node') # get the node importance
EMP_network_analysis <- function(obj,select=NULL,method='cor',corMethod='spearman',weighted=TRUE,
                              missing='listwise',threshold='sig',coldata_to_assay=NULL,use_cached=TRUE,action='add',...) {

  call <- match.call()
  
  if (use_cached == FALSE) {
    memoise::forget(.EMP_network_analysis_EMP_m) %>% invisible()
    memoise::forget(.EMP_network_analysis_EMPT_m) %>% invisible()
  }
  
  if(method == 'graphicalVAR') {
    rlang::check_installed(c('BiocManager'), reason = 'for estimateNetwork().', action = install.packages)  
    rlang::check_installed(c('graphicalVAR'), reason = 'for estimateNetwork().', action = BiocManager::install)
  }
 
  if(method == 'GGMncv') {
    rlang::check_installed(c('BiocManager'), reason = 'for estimateNetwork().', action = install.packages)  
    rlang::check_installed(c('GGMncv'), reason = 'for estimateNetwork().', action = BiocManager::install)
  } 
  
  if(method == 'IsingSampler') {
    rlang::check_installed(c('BiocManager'), reason = 'for estimateNetwork().', action = install.packages)  
    rlang::check_installed(c('IsingSampler'), reason = 'for estimateNetwork().', action = BiocManager::install)
  } 

  if (is(obj,"EMP")) {
      result <- .EMP_network_analysis_EMP_m(EMP=obj,select=select,method=method,corMethod=corMethod,weighted=weighted,
                              missing=missing,threshold=threshold,...)   

      .get.history.EMP(result) <- call

      if (action=='get') {
        return(result@deposit[['EMP_network_analysis']])
      }else if(action=='add') {
        return(result)
      }else{
        warning('action should be one of add or get!')
      }
  }else if (is(obj,"EMPT")) {
      result <- .EMP_network_analysis_EMPT_m(EMPT=obj,method=method,corMethod=corMethod,weighted=weighted,
                              missing=missing,threshold=threshold,coldata_to_assay=coldata_to_assay,...) 
      .get.history.EMPT(result) <- call

      if (action=='get') {
        deposit <- list()
        deposit[['net']] <- result@deposit[['net']] 
        deposit[['net_feature_info']] <- result@deposit[['net_feature_info']] 
        deposit[['net_centrality']] <- result@deposit[['net_centrality']] 
        return(deposit)
      }else if(action=='add') {
        return(result)
      }else{
        warning('action should be one of add or get!')
      }                                
  }  
}








