#' @importFrom bigstatsr big_SVD
#' @importFrom vegan vegdist
#' @importFrom stats cmdscale
.EMP_dimension_analysis <- function(EMPT,method,distance=NULL,estimate_group=NULL){
  p1 <- p2 <- p3 <- R2X <- NULL
  deposit <- list()
  assay_data  <- EMPT %>%
    .get.assay.EMPT() %>% tibble::column_to_rownames('primary')
  
  switch(method,
         "pca" = {
           if(!is.null(distance)){
             message("Parameter distance in pca is euclidean!")
             distance <- 'euclidean'
           }
           sample_name <- rownames(assay_data)
           assay_data <- assay_data |> bigstatsr::as_FBM()
           pca_result <- bigstatsr::big_SVD(assay_data,k = 3)
           colnames(pca_result$u) <- c('PC1','PC2','PC3')
           dimension_reslut <- pca_result$u |> 
             tibble::as_tibble() |> 
             dplyr::mutate(primary = {{sample_name}},.before = 1) 
           axis_value <- round( (pca_result$d/sum(pca_result$d)) * 100, digits = 2)
         },
         "pcoa" = {
           if(is.null(distance)){
             stop("Parameter distance in necessary!")
           }
           assay_data_dis <-  vegan::vegdist(assay_data,method=distance) %>% as.matrix()
           
           # set the warning
           check_dim <- dim(assay_data)[1] * dim(assay_data)[2]
           if (check_dim > 1e+05) {
             message_wrap("Inputting large-scale data may lead to longer computation time. It is recommended to use other dimensionality reduction methods or filter the data first.")
           }
           
           pca_result <- stats::cmdscale(assay_data_dis,k=3,eig=T)
           colnames(pca_result$points) <- c('PCoA1','PCoA2','PCoA3')
           dimension_reslut <- pca_result$points |> 
             as.data.frame() |>
             tibble::rownames_to_column('primary') 
           axis_value <- round((pca_result$eig / sum(pca_result$eig))[1:3] * 100,digits = 2)
         },
         "pls" = {
           rlang::check_installed(c('BiocManager'), reason = 'for EMP_dimension_analysis().', action = install.packages) 
           rlang::check_installed(c('ropls'), reason = 'for EMP_dimension_analysis().', action = BiocManager::install)            
           if(!is.null(distance)){
             message("Parameter distance in pls is useless!")
             distance <- NULL
           }
           estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)
           group_info <- EMPT %>% EMP_coldata_extract() %>% dplyr::pull(!!estimate_group)
           
           plsda_raw <- ropls::opls(assay_data,
                                    group_info,
                                    orthoI=0,
                                    predI=3,
                                    scaleC = 'none',fig.pdfC='none',info.txtC='none')
           dimension_reslut <- plsda_raw@scoreMN %>% as.data.frame() %>%
             tibble::rownames_to_column('primary') %>% tibble::as_tibble() %>%
             dplyr::rename(t1 = p1,t2 =p2, t3 =p3)
           vip_score <- plsda_raw@vipVn %>% as.data.frame() %>%
             tibble::rownames_to_column('feature') %>% tibble::as_tibble()
           colnames(vip_score) <- c('feature','VIP')
           axis_value <- plsda_raw@modelDF  %>% dplyr::pull(R2X)
           axis_value <- axis_value * 100
           
         },
         "opls" = {
           rlang::check_installed(c('BiocManager'), reason = 'for EMP_dimension_analysis().', action = install.packages) 
           rlang::check_installed(c('ropls'), reason = 'for EMP_dimension_analysis().', action = BiocManager::install)       
           if(!is.null(distance)){
             message("Parameter distance in opls is useless!")
             distance <- NULL
           }
           estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)
           group_info <- EMPT %>% EMP_coldata_extract() %>% dplyr::pull(!!estimate_group)
           
           if(length(unique(group_info)) > 2){
             stop("opls model ony support two groups analysis!")
           }
           
           opsda_raw <- ropls::opls(assay_data,
                                    group_info,
                                    orthoI=3,
                                    predI=1,scaleC = 'none',fig.pdfC='none',info.txtC='none')
           
           
           if(!is.null(opsda_raw) & length(opsda_raw@modelDF) == 0) {
             message_wrap("No model was built because the first predictive component was already not significant!")
             message_wrap("Filter the input data may improve the issue!")
             stop()
           }
           
           dimension_reslut <-  opsda_raw@scoreMN %>%  #得分矩阵
             as.data.frame() %>%
             dplyr::mutate(o1 = opsda_raw@orthoScoreMN[,1]) %>%
             dplyr::rename(t1=p1) %>%
             tibble::rownames_to_column('primary')
           try(dimension_reslut %<>% dplyr::mutate(o2 = opsda_raw@orthoScoreMN[,2]),silent = T)
           
           vip_score <- opsda_raw@vipVn %>% as.data.frame() %>%
             tibble::rownames_to_column('feature') %>% tibble::as_tibble()
           colnames(vip_score) <- c('feature','VIP')
           axis_value <- opsda_raw@modelDF  %>% dplyr::pull(R2X)
           axis_value <- axis_value * 100

         },
         {
           print('method in EMP_reduce_dimension should be pca, pcoa, pls or opls!')
         }
  )
  

  EMPT@deposit[['dimension_coordinate']] <- dimension_reslut
  EMPT@deposit[['dimension_axis']] <- axis_value
  estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group) ## make sure that pca and pcoa could inherit the estimate_group in the previous step.
  .get.estimate_group.EMPT(EMPT) <- estimate_group
  try(EMPT@deposit[['dimension_VIP']] <- vip_score,silent = T)


  .get.method.EMPT(EMPT) <- method
  .get.algorithm.EMPT(EMPT) <- distance
  .get.info.EMPT(EMPT) <- 'EMP_dimension_analysis'
  EMPT
}

.EMP_dimension_analysis_m <-memoise::memoise(.EMP_dimension_analysis)

#' Dimension reduction of the abundance or experssion data
#'
#' @param x Object in EMPT or MultiAssayExperiment format.
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object.
#' @param method  A character string. Methods include pca, pcoa, pls, opls.
#' @param distance A character string. The logarithm distance used in method = "pcoa".Detailed in the vegan::vegdist.
#' @param estimate_group A character string. Select the group name in the coldata to be considerated.
#' @param action A character string.A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
#'
#' @return EMPT object
#' @export
#'
#' @examples
#' data(MAE)
#' ## PCA
#' MAE |>
#'   EMP_dimension_analysis(experiment = 'taxonomy',method = 'pca') 
#' ## Pcoa
#' MAE |>
#'   EMP_dimension_analysis(experiment = 'taxonomy',method = 'pcoa',distance = 'bray')
#' ## Pls
#' MAE |>
#'   EMP_dimension_analysis(experiment = 'untarget_metabol',
#'                          method = 'pls',estimate_group = 'Group')
#' MAE |>
#'   EMP_collapse(experiment = 'untarget_metabol',na_string=c('NA','null','','-'),
#'                estimate_group = 'MS2kegg',method = 'sum',collapse_by = 'row') |> ## Get the kegg compound
#'   EMP_dimension_analysis(method = 'pls',estimate_group = 'Group')
#' ## OPLS
#' MAE |>
#'   EMP_collapse(experiment = 'untarget_metabol',na_string=c('NA','null','','-'),
#'                estimate_group = 'MS2kegg',method = 'sum',collapse_by = 'row') |>
#'   EMP_diff_analysis(method='DESeq2',.formula = ~Group) |>
#'   EMP_filter(feature_condition = pvalue < 0.05) |>
#'   EMP_dimension_analysis(method = 'opls',estimate_group = 'Sex') |>
#'   EMP_scatterplot(estimate_group='Sex',show='p12html',ellipse=0.6) ## Visualization
EMP_dimension_analysis <- function(x,experiment,method='pcoa',distance=NULL,use_cached=TRUE,
                                   estimate_group=NULL,action='add'){
  call <- match.call()

  if (inherits(x,"MultiAssayExperiment")) {
    EMPT <- .as.EMPT(x,
                     experiment = experiment)
  }else if(inherits(x,'EMPT')) {
    EMPT <-x
  }

  if (use_cached == F) {
    memoise::forget(.EMP_dimension_analysis_m) %>% invisible()
  }

  EMPT %<>% .EMP_dimension_analysis_m(method=method,distance=distance,estimate_group=estimate_group)
  .get.history.EMPT(EMPT) <- call
  class(EMPT) <- 'EMP_dimension_analysis'


  if (action == 'add') {
    return(EMPT)
  }else if (action == 'get') {
    deposit <- list()
    deposit[['dimension_coordinate']] <- EMPT@deposit$dimension_coordinate
    deposit[['dimension_axis']] <- EMPT@deposit$dimension_axis
    try(deposit[['dimension_VIP']] <- EMPT@deposit$dimension_VIP,silent = T)
    return(deposit)
  }else{
    warning('action should be one of add or get!')
  }

}




