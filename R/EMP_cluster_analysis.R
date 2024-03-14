#' Title
#'
#' @param x wait_for_add
#' @param experiment wait_for_add
#' @param distance wait_for_add
#' @param rowdata wait_for_add
#' @param method wait_for_add
#' @param h wait_for_add
#' @param groupLabels wait_for_add
#' @param cluster wait_for_add
#' @param action wait_for_add
#' @param ... wait_for_add
#' @importFrom tibble column_to_rownames
#' @importFrom stats hclust
#' @importFrom vegan vegdist
#' @importFrom dendextend color_branches
#' @importFrom dendextend cutree
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_cluster_analysis <- function(x,experiment,distance='bray',rowdata=F,
                                 method='average',h=NULL,groupLabels=T,
                                 cluster=NULL,action='add',...) {
  colname <- primary <- NULL
  call <- match.call()
  if (inherits(x,"MultiAssayExperiment")) {

    EMPT <- .as.EMPT(x,
                     experiment = experiment)

  }else if(inherits(x,'EMPT')) {
    EMPT <- x
  }else {
    stop('Please check the input data!')
  }

  if (rowdata == F) {
    data_distance <- EMPT %>%
      .get.assay.EMPT() %>% tibble::column_to_rownames('primary') %>%
      vegan::vegdist(method=distance)
  }else if (rowdata == T) {
    data_distance <- EMPT %>%
      .get.assay.EMPT() %>% tibble::column_to_rownames('primary') %>% t() %>%
      vegan::vegdist(method=distance)
  }



  Tree = stats::hclust(as.dist(data_distance), method = method)

  dendrogram_obj <- Tree %>%
    as.dendrogram()

  if (is.null(h)) {
    dendrogram_obj %>%
      dendextend::color_branches(.,
                                 # h = 1500000,
                                 col = 'black',
                                 groupLabels = F,warn=F) %>% plot(main = "Clustering to detect outliers") %>%
      suppressWarnings()

  }else{
    dendrogram_obj %>%
      dendextend::color_branches(.,
                                 h = h,
                                 groupLabels = T,warn=F) %>% plot(main = "Clustering to detect outliers")
    abline(h = h, col = "red",lty=2)

    desired_branch <- dendextend::cutree(
      dendrogram_obj,
      h = h,
      order_clusters_as_data = F)

    if (rowdata == F) {
      cluster_result <- desired_branch  %>%
        as.data.frame() %>% tibble::rownames_to_column('primary') %>%
        dplyr::rename(cluster = '.') %>% tibble::as_tibble()

      EMPT@deposit[['sample_cluster_result']] <- cluster_result
    }else if (rowdata == T) {
      cluster_result <- desired_branch  %>%
        as.data.frame() %>% tibble::rownames_to_column('feature') %>%
        dplyr::rename(cluster = '.') %>% tibble::as_tibble()
      EMPT@deposit[['feature_cluster_result']] <- cluster_result
    }
 }

 if (action=='add') {
    EMPT@method <- paste0(distance,'+',method)
    EMPT@algorithm <- 'cluster_analysis'
    .get.info.EMPT(EMPT) <- 'EMP_cluster_analysis'
    .get.history.EMPT(EMPT) <- call
    class(EMPT) <- 'EMP_cluster_analysis'
    return(EMPT)
 }else if (action=='get') {
    return(cluster_result)
 }else{
    stop('action should be one of add or get!')
 }
}
