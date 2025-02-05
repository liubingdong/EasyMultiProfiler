#' @importFrom purrr reduce
#' @importFrom bootnet estimateNetwork
#' @importFrom qgraph centralityTable
.EMP_network_analysis <- function(EMP,select=NULL,method='cor',corMethod='spearman',weighted=TRUE,
                              missing='pairwise',threshold=0,...){
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

  
  network <- estimateNetwork(assay_df,
                             threshold = threshold,
                             default = method,
                             corMethod = corMethod,
                             missing = missing,
                             weighted=weighted,
                             nonPositiveDefinite='continue',...) |> 
             spsUtil::quiet(print_cat = TRUE, message = TRUE, warning = TRUE)

  centrality_info <- centralityTable(network)


  EMP@deposit[['network_result']][['net']] <- network
  EMP@deposit[['network_result']][['feature_info']] <- rowdata_list
  EMP@deposit[['network_result']][['centrality']] <- centrality_info

  .get.info.EMP(EMP) <-'EMP_network_analysis'
  .get.method.EMP(EMP) <- method
  class(EMP) <- 'EMP_network_analysis'
  return(EMP)
}

.EMP_network_analysis_m <- memoise::memoise(.EMP_network_analysis,cache = cachem::cache_mem(max_size = 4096 * 1024^2))


.network_print <- function (obj){
  if (is(obj,'EMP')) {
    network <- obj@deposit[['network_result']][['net']]
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
#' @param EMP Object in EMP format.
#' @param select A character string. The experiment name in the EMP object.
#' @param method A character string. Methods include cor, EBICglasso,IsingFit,IsingSampler,pcor,huge,mgm,TMFG,LoGo,relimp,ggmModSelect and graphicalVAR.
#' @param corMethod A character string. Methods include spearman, cor, cov, npn and cor_auto.
#' @param missing How to handle missing data? "pairwise" for pairwise deletion, "listwise" for listwise deletion, "fiml" for full-information maximum likelihood and "stop" to stop with an error.
#' @param threshold Thresholding to use in partial correlation networks. Can be a fixed number to threshold all absolute edges below this value, 'locfdr' for local FDR, or any option corresponding to adjustments in \code{\link[psych]{corr.p}} ('none', 'sig', 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY' or 'fdr')). Can also be used for default = "IsingSampler" but can only be set to a logical enabling or disabling significance thresholding.
#' @param action A character string.A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
#' @param ... Further parameters passed to the function \code{\link[bootnet]{estimateNetwork}}.
#' @rdname EMP_network_analysis
#' @section Detaild about network:
#' More detailed information are available on :
#'
#' Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating psychological networks and their accuracy: A tutorial paper. Behavior Research Methods, 50(1), 195–212. 
#'
#' @return EMP object
#' @export
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
#' (k1 + k2 ) |> 
#'   EMP_network_analysis() |> 
#'   EMP_network_plot(show = 'node',layout = 'spring',
#'                    vsize = 5,
#'                    node_info = c('Phylum','MS2class'),
#'                    legend.cex=0.3,label.cex = 1,label.prop = 0.9,
#'                    fade = T,shape = 'circle',font=2)
#' (k1 + k2 ) |> 
#'   EMP_network_analysis() |> 
#'   EMP_network_plot(show = 'node') # get the node importance
EMP_network_analysis <- function(EMP,select=NULL,method='cor',corMethod='spearman',weighted=TRUE,
                              missing='pairwise',threshold=0,use_cached=TRUE,action='add',...) {


  rlang::check_installed(c('BiocManager'), reason = 'for EMP_network_analysis().', action = install.packages) 
  rlang::check_installed(c('bootnet'), reason = 'for EMP_network_analysis().', action = BiocManager::install)

  call <- match.call()
  
  if (!is(EMP,"EMP")) {
    stop("Please input the EMP format!")
  }

  if (use_cached == FALSE) {
    memoise::forget(.EMP_network_analysis_m) %>% invisible()
  }
  
  result <- .EMP_network_analysis_m(EMP=EMP,select=select,method=method,corMethod=corMethod,weighted=weighted,
                              missing=missing,threshold=threshold,...)


 .get.history.EMP(result) <- call

  if (action=='get') {
    return(result@deposit[['network_result']])
  }else if(action=='add') {
    return(result)
  }else{
    warning('action should be one of add or get!')
  }
}








