#' Identify clusters for sample or feature
#'
#' @param x Object in EMPT or MultiAssayExperiment format.
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object.
#' @param distance A character string.Dissimilarity index, partial match to "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger", "aitchison", or "robust.aitchison".
#' @param rowdata A boolean. Whether the function cluster the feature or not.
#' @param method The agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param h A numeric. Height at which to cut tree (passed to cutree)
#' @param groupLabels A boolean. Whether show the group label or not.
#' @param action  A character string.A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @importFrom tibble column_to_rownames
#' @importFrom stats hclust
#' @importFrom vegan vegdist
#' @importFrom graphics abline
#' @importFrom BiocManager install
#' @importFrom utils install.packages
#'
#' @return EMPT object
#' @export
#'
#' @examples
#' data(MAE)
#' ## Cluster the samples according to the assay data
#' MAE |>
#'   EMP_assay_extract(experiment = 'geno_ec') |>
#'   EMP_cluster_analysis()
#' 
#' MAE |>
#'   EMP_assay_extract(experiment = 'geno_ec') |>
#'   EMP_cluster_analysis(h=0.15) # identify the outlier samples
#' 
#' MAE |>
#'   EMP_assay_extract(experiment = 'geno_ec') |>
#'   EMP_cluster_analysis(h=0.15) |>
#'   EMP_filter(sample_condition = cluster != 1) # filter away the outlier
#' 
#' ## Cluster the samples according to the coldata
#' MAE |> 
#'   EMP_coldata_extract(action = 'add') |> # transfer the coldata to asaay
#'   EMP_impute(assay = T) |> # impute missing value
#'   EMP_cluster_analysis(method = 'ward.D2',distance='clark',h=0.2) 
#' 
#' ## Cluster the features according to the assay data
#' MAE |> 
#'   EMP_assay_extract(experiment = 'geno_ec',pattern='1.1.1.1',pattern_ref='feature') |>
#'   EMP_cluster_analysis(rowdata = T,h=0.8)
EMP_cluster_analysis <- function(x,experiment,distance='bray',rowdata=FALSE,
                                 method='average',h=NULL,groupLabels=TRUE,action='add') {
  
  # Check if package is installed, otherwise install
  if (find.package("dendextend", quiet = TRUE) %>% length() == 0) {
    message("EMP_cluster_analysis need install package dendextend!")
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install("dendextend", ask = FALSE)
  }  

  colname <- primary <- `.` <- NULL
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
