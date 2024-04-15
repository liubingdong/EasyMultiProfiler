#' @importFrom dplyr arrange
#' @importFrom tibble column_to_rownames
#' @noRd
.EMP_cor_analysis <- function(EMP,select=NULL,method='spearman',...) {
  #call <- match.call()

  var1 <- NULL
  if (is.null(select)) {
    data1 <- EMP@ExperimentList[[1]] %>% EMP_assay_extract(action='get') %>%
      dplyr::arrange('primary') %>% tibble::column_to_rownames('primary') %>% suppressMessages()
    if (length(EMP@ExperimentList) == 1) {
      data2 <- EMP@ExperimentList[[1]] %>% EMP_assay_extract(action='get') %>%
        dplyr::arrange('primary') %>% tibble::column_to_rownames('primary') %>% suppressMessages()
      cor_info <- names(EMP@ExperimentList)
    }else{
      data2 <- EMP@ExperimentList[[2]] %>% EMP_assay_extract(action='get') %>%
        dplyr::arrange('primary') %>% tibble::column_to_rownames('primary') %>% suppressMessages()
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
        dplyr::arrange('primary') %>% tibble::column_to_rownames('primary') %>% suppressMessages()
    data2 <- EMP@ExperimentList[[select[2]]] %>% EMP_assay_extract(action='get') %>%
        dplyr::arrange('primary') %>% tibble::column_to_rownames('primary') %>% suppressMessages()

  }

  real_sample <- intersect(rownames(data1),rownames(data2))
  data1_sample_num <- rownames(data1) %>% unique %>% length()
  data2_sample_num <- rownames(data2) %>% unique %>% length()

  data1 <- data1 |> dplyr::filter(rownames(data1) %in% real_sample)
  data2 <- data2 |> dplyr::filter(rownames(data2) %in% real_sample)
  df.cor.p<-agricolae_correlation(x=data1,y=data2,method = method,...)

  df <- df.cor.p$correlation %>%
    as.data.frame() %>%
    tibble::rownames_to_column('var1') %>%
    tidyr::pivot_longer(cols = -var1,names_to = 'var2',values_to = 'coefficient')

  df.p <- df.cor.p$pvalue %>%
    as.data.frame() %>%
    tibble::rownames_to_column('var1') %>%
    tidyr::pivot_longer(cols = -var1,names_to = 'var2',values_to = 'pvalue')
  df$pvalue<-round(df.p$pvalue,2)


  df.cor.p[['cor_info']] <- cor_info
  df.cor.p[["n.obs"]] <- c(data1_sample_num,data2_sample_num,length(real_sample))
  df.cor.p[['cor_p']] <- df
  EMP@deposit[['cor_analysis_result']] <- df.cor.p
  #.get.history.EMP(EMP) <- call
  .get.info.EMP(EMP) <-'EMP_cor_analysis'
  .get.method.EMP(EMP) <- method
  class(EMP) <- 'EMP_cor_analysis'
  return(EMP)
}

.EMP_cor_analysis_m <- memoise::memoise(.EMP_cor_analysis)


#' Correlation analysis
#'
#' @param EMP Object in EMP format.
#' @param select A character string. The experiment name in the EMP object.
#' @param method A character string. Methods include pearson (default), spearman and kendall.
#' @param action A character string.A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
#' @param ... Further parameters passed to the function agricolae::correlation
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
#' (k1 + k2) |> EMP_cor_analysis(method = 'spearman') |>
#'   EMP_heatmap_plot() ## Visualization
EMP_cor_analysis <- function(EMP,select=NULL,method='spearman',action='add',
                             use_cached = TRUE,...) {
  call <- match.call()
  
  if (!inherits(EMP,"EMP")) {
    stop("Please input the EMP format!")
  }

  if (use_cached == F) {
    memoise::forget(.EMP_cor_analysis_m) %>% invisible()
  }

 result <- .EMP_cor_analysis_m(EMP=EMP,select=select,method=method,...)
 .get.history.EMP(result) <- call
  if (action=='get') {
    return(result@deposit[['cor_analysis_result']])
  }else if(action=='add') {
    return(result)
  }else{
    warning('action should be one of add or get!')
  }

}





