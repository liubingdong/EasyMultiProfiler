#' @importFrom bigstatsr big_SVD
#' @importFrom vegan vegdist
#' @importFrom stats cmdscale
.EMP_dimension_analysis <- function(EMPT,method,distance=NULL,estimate_group=NULL,umap.config=NULL,
                                    scale=NULL,bySample='default',logbase =2,pseudocount=0.0000001,...){
  p1 <- p2 <- p3 <- R2X <- NULL
  deposit <- list()

  if (!is.null(scale)) {
    EMPT_trans <- EMPT |> EMP_decostand(method=scale,bySample=bySample,logbase =logbase,pseudocount=pseudocount)
    assay_data  <- assay(EMPT_trans) %>% t()
  }else{
    assay_data  <- assay(EMPT) %>% t()
  }
  
  switch(method,
         "pca" = {
           if(!is.null(distance)){
             EMP_message("Parameter distance in pca is euclidean!",color = 32,order = 1,show='message')
           }
           distance <- 'euclidean'
           sample_name <- rownames(assay_data)
           assay_data <- assay_data |> bigstatsr::as_FBM()
           pca_result <- bigstatsr::big_SVD(assay_data,k = 3,...)
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
           assay_data_dis <-  vegan::vegdist(assay_data,method=distance,na.rm = TRUE) %>% as.matrix() %>% suppressWarnings()
           
           # set the warning
           check_dim <- dim(assay_data)[1] * dim(assay_data)[2]
           check_sample_num <- dim(assay_data)[1]
           if (check_dim > 1e+05 | check_sample_num > 500) {
            EMP_message("Large-scale data need longer computation time in PCoA.",color = 32,order = 1,show='message')
           }
           
           pca_result <- stats::cmdscale(assay_data_dis,k=3,eig=T,...)
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
             EMP_message("Parameter distance in pls is useless!",color = 32,order = 1,show='message')
             distance <- NULL
           }
           estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)
           group_info <- EMPT %>% EMP_coldata_extract() %>% dplyr::pull(!!estimate_group)
           
           ## check the missing value in the group label
           if(any(is.na(group_info))) {
             stop('Column ',estimate_group,' has beed deteced missing value, please check and filter them!')
           }

           plsda_raw <- ropls::opls(assay_data,
                                    group_info,
                                    orthoI=0,
                                    predI=3,
                                    scaleC = 'none',fig.pdfC='none',info.txtC='none',...)
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
             EMP_message("Parameter distance in opls is useless!",color = 32,order = 1,show='message')
             distance <- NULL
           }
           estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)
           group_info <- EMPT %>% EMP_coldata_extract() %>% dplyr::pull(!!estimate_group)
           
           ## check the missing value in the group label
           if(any(is.na(group_info))) {
             stop('Column ',estimate_group,' has beed deteced missing value, please check and filter them!')
           }
           
           if(length(unique(group_info)) > 2){
             stop("opls model ony support two groups analysis!")
           }
           
           opsda_raw <- ropls::opls(assay_data,
                                    group_info,
                                    orthoI=3,
                                    predI=1,scaleC = 'none',fig.pdfC='none',info.txtC='none',...)
           
           
           if(!is.null(opsda_raw) & length(opsda_raw@modelDF) == 0) {
             EMP_message("No model was built because the first predictive component was not significant!\nFilter the input data may improve the issue!",color = 31,order = 1,show='warning')
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
         "umap" = {
          rlang::check_installed(c('BiocManager'), reason = 'for EMP_dimension_analysis().', action = install.packages) 
          rlang::check_installed(c('umap'), reason = 'for EMP_dimension_analysis().', action = BiocManager::install)  

          if(!is.null(umap.config)){
            umap.config <- umap.config
          }else{
            umap.config <- umap::umap.defaults
          }

          if(!is.null(distance)){
            umap.config$metric <- distance
          }
          # set 3 axis
          umap.config$n_components <- 3
          
          sample_name <- rownames(assay_data)
          umap_result <-umap::umap(assay_data,config = umap.config,...)
          dimension_reslut <- umap_result[['layout']] %>% as.data.frame %>%
             dplyr::mutate(primary = sample_name,.before=1)
          colnames(dimension_reslut) <- c('primary',paste0('umap',1:(dim(dimension_reslut)[2]-1)))
          
          axis_value <- NULL

         },
         {
           stop('Parameter method in EMP_reduce_dimension should be pca, pcoa, pls, opls or umap!')
         }
  )
  

  EMPT@deposit[['dimension_coordinate']] <- dimension_reslut
  EMPT@deposit[['dimension_axis']] <- axis_value

  ## make sure that pca and pcoa could inherit the estimate_group in the previous step. 
  if (method %in% c('pls','opls')) {
      .get.estimate_group.EMPT(EMPT) <- estimate_group
  }
  try(EMPT@deposit[['dimension_VIP']] <- vip_score,silent = T)
  .get.method.EMPT(EMPT) <- method
  .get.algorithm.EMPT(EMPT) <- distance
  .get.info.EMPT(EMPT) <- 'EMP_dimension_analysis'
  EMPT
}

.EMP_dimension_analysis_m <-memoise::memoise(.EMP_dimension_analysis,cache = cachem::cache_mem(max_size = 4096 * 1024^2))

#' Dimension reduction of the abundance or experssion data
#'
#' @param obj Object in EMPT or MultiAssayExperiment format.
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object.
#' @param method  A character string. Methods include pca, pcoa, pls, opls, umap.
#' @param distance A character string. The logarithm distance used in method = "pcoa".Detailed in the vegan::vegdist including "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger", "aitchison", or "robust.aitchison".
#' @param estimate_group A character string. Select the group name in the coldata to be considerated.
#' @param umap.config A list. only activate in umap. More see umap::umap.
#' @param scale   A character string. The parameter works in the same way as the method in EMP_decostand.
#' @param bySample A boolean. Whether the function decostand by the sample or feature. Detaled information in the EMP_decostand.
#' @param logbase An interger. The logarithm base used in method = "log".(default=2). Detaled information in the EMP_decostand.
#' @param pseudocount A number. The logarithm pseudocount used in method = "clr" or "alr".(default=0.0000001). 
#' @param action A character string.A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
#' @param ... Additional parameters, see also \code{\link[bigstatsr]{big_SVD}} for pca, \code{\link[stats]{cmdscale}} for pcoa,\code{\link[ropls]{opls}} for pls and opls, \code{\link[umap]{umap}} for umap.
#' @section Detaild about scale:
#' When the scale parameter is enabled, data normalization only takes effect within this function and does not affect the original data in the EMPT. 
#' If you need to transform the data, consider using the EMP_EMP_decostand function in the workflow.
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
#' ## umap
#' MAE |>
#'   EMP_dimension_analysis(experiment = 'geno_ko',
#'                          method = 'umap') |>
#'   EMP_scatterplot(estimate_group='Group',show='p12html')
#' 
#' ### modify umap contig
#' umap.config <- umap::umap.defaults
#' umap.config$n_neighbors=10
#' MAE |>
#'   EMP_dimension_analysis(experiment = 'geno_ko',
#'                          method = 'umap',umap.config = umap.config) 

EMP_dimension_analysis <- function(obj,experiment,method='pcoa',distance=NULL,use_cached=TRUE,
                                  scale=NULL,bySample='default',logbase =2,pseudocount=0.0000001,
                                  estimate_group=NULL,umap.config=NULL,action='add',...){
  call <- match.call()

  if (is(obj,"MultiAssayExperiment")) {
    EMPT <- .as.EMPT(obj,
                     experiment = experiment)
  }else if(is(obj,'EMPT')) {
    EMPT <- obj
  }

  if (use_cached == FALSE) {
    memoise::forget(.EMP_dimension_analysis_m) %>% invisible()
  }

  EMPT <- EMPT |> .EMP_dimension_analysis_m(method=method,distance=distance,estimate_group=estimate_group,umap.config=umap.config,
                                            scale=scale,bySample=bySample,logbase =logbase,pseudocount=pseudocount,...)
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
    warning('Parameter action should be one of add or get!')
  }

}




