#' Title
#'
#' @param EMPT wait_for_add
#' @param method wait_for_add
#' @param distance wait_for_add
#' @param estimate_group wait_for_add
#' @importFrom vegan vegdist
#' @importFrom ape pcoa
#' @importFrom ropls opls
#' @noRd
.EMP_dimension_analysis <- function(EMPT,method,distance=NULL,estimate_group=NULL){
  deposit <- list()
  assay_data  <- EMPT %>%
    .get.assay.EMPT() %>% tibble::column_to_rownames('primary')

  switch(method,
         "pca" = {
           distance <- 'euclidean' ## PCA is one of pcoa using euclidean distance!
           assay_data_dis <-  vegan::vegdist(assay_data,method=distance) %>% as.matrix()
           pca_result <- ape::pcoa(assay_data_dis, correction = "none", rn = NULL)
           dimension_reslut <- data.frame(PC1=pca_result$vectors[,1],
                                          PC2=pca_result$vectors[,2],
                                          PC3=pca_result$vectors[,3]) %>%
             tibble::rownames_to_column('primary') %>% tibble::as_tibble()
           axis_value <- c(round(pca_result$values$Relative_eig[1]*100,digits = 2),
            round(pca_result$values$Relative_eig[2]*100,digits = 3),
            round(pca_result$values$Relative_eig[3]*100,digits = 3))

         },
         "pcoa" = {
           if(is.null(distance)){
             stop("Parameter distance in necessary!")
           }
           assay_data_dis <-  vegan::vegdist(assay_data,method=distance) %>% as.matrix()
           pca_result <- ape::pcoa(assay_data_dis, correction = "none", rn = NULL)
           dimension_reslut <- data.frame(PC1=pca_result$vectors[,1],
                                          PC2=pca_result$vectors[,2],
                                          PC3=pca_result$vectors[,3]) %>%
             tibble::rownames_to_column('primary') %>% tibble::as_tibble()
           axis_value <- c(round(pca_result$values$Relative_eig[1]*100,digits = 2),
            round(pca_result$values$Relative_eig[2]*100,digits = 3),
            round(pca_result$values$Relative_eig[3]*100,digits = 3))
         },
         "pls" = {
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
                             orthoI=NA,
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
  .get.estimate_group.EMPT(EMPT) <- estimate_group
  try(EMPT@deposit[['dimension_VIP']] <- vip_score,silent = T)


  .get.method.EMPT(EMPT) <- method
  .get.algorithm.EMPT(EMPT) <- distance
  .get.info.EMPT(EMPT) <- 'EMP_dimension_analysis'
  EMPT
}

.EMP_dimension_analysis_m <-memoise::memoise(.EMP_dimension_analysis)

#' Title
#'
#' @param x wait_for_add
#' @param experiment wait_for_add
#' @param method wait_for_add
#' @param distance wait_for_add
#' @param estimate_group wait_for_add
#' @param action wait_for_add
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_dimension_analysis <- function(x,experiment,method='pcoa',distance=NULL,
                                   estimate_group=NULL,action='add'){
  call <- match.call()

  if (inherits(x,"MultiAssayExperiment")) {
    EMPT <- .as.EMPT(x,
                     experiment = experiment)
  }else if(inherits(x,'EMPT')) {
    EMPT <-x
  }
  #estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)
  EMPT %<>% .EMP_dimension_analysis_m(method=method,distance=distance,estimate_group=estimate_group)
  .get.history.EMPT(EMPT) <- call
  class(EMPT) <- 'EMP_dimension_analysis'


  if (action == 'add') {
    return(EMPT)
  }else if (action == 'get') {
    deposit <- list()
    deposit[['dimension_coordinate']] <- EMPT@deposit$dimension_coordinate
    deposit[['dimension_axis']] <- EMPT@deposit$dimension_axis
    try(EMPT@deposit[['dimension_VIP']] <- vip_score,silent = T)
    return(deposit)
  }else{
    warning('action should be one of add or get!')
  }

}




