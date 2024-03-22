#' @importFrom WGCNA labels2colors
.EMP_WGCNA_cluster_analysis <- function(EMPT,powers=c(1:10, seq(from = 12, to=20, by=2)),
                                        RsquaredCut=0.85, removeFirst = FALSE, nBreaks = 10, blockSize = NULL,
                                        # corFnc = WGCNA::cor, corOptions = list(use = 'p'),
                                        networkType = "unsigned",
                                        moreNetworkConcepts = FALSE,
                                        gcInterval = NULL,
                                        TOMType = "unsigned", minModuleSize = 30,
                                        reassignThreshold = 0, mergeCutHeight = 0.25,
                                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                                        # saveTOMs = T,...) {
                                        saveTOMs = T) {
  deposit <- list()
  Var1 <- Freq <- NULL
  #enableWGCNAThreads()
  
  assay_data <- .get.assay.EMPT(EMPT) %>%
    tibble::column_to_rownames('primary') %>%
    as.matrix()
  # print("here 2_1 over")
  sft <- spsUtil::quiet(WGCNA::pickSoftThreshold(assay_data, powerVector = powers,
                          RsquaredCut=RsquaredCut, removeFirst = removeFirst, nBreaks = nBreaks, blockSize = blockSize,
                          # corFnc = corFnc, corOptions = corOptions,
                          networkType = networkType,
                          moreNetworkConcepts = moreNetworkConcepts,
                          gcInterval = gcInterval),print_cat = TRUE, message = TRUE, warning = TRUE)
  # print("here 2_2 over")
  check_best_power <- sft[['powerEstimate']]

  if (is.na(check_best_power)) {
    stop("No best soft power,please reset the RsquaredCut or other parameter!")
  }
  save(assay_data, check_best_power, TOMType, minModuleSize, reassignThreshold, mergeCutHeight, numericLabels, pamRespectsDendro, saveTOMs, file = "test.Rdata")
  # net = WGCNA::blockwiseModules(assay_data, power = check_best_power,
  net = WGCNA_blockwiseModules(assay_data, power = check_best_power,
                         TOMType = TOMType, minModuleSize = minModuleSize,
                         reassignThreshold = reassignThreshold, mergeCutHeight = mergeCutHeight,
                         numericLabels = numericLabels, pamRespectsDendro = pamRespectsDendro,
                         # saveTOMs = saveTOMs,...)
                         saveTOMs = saveTOMs)

  mergedColors = WGCNA::labels2colors(net$colors)
  # print("here 2_3 over")
  WGCNA::plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      abHeight=mergeCutHeight)

  mergedColors = WGCNA::labels2colors(net$colors)
  # print("here 2_4 over")
  WGCNA_module_elements <- table(WGCNA::labels2colors(net$colors)) %>%
                as.data.frame() %>%
                dplyr::rename(WGCNA_color=Var1,WGCNA_module_elements=Freq)
  # print("here 2_5 over")
  feature_modules <- data.frame(WGCNA_cluster=net$unmergedColors,WGCNA_color=mergedColors) %>% tibble::rownames_to_column('feature') %>%
    tibble::as_tibble() %>% dplyr::left_join(WGCNA_module_elements,by='WGCNA_color')


  #disableWGCNAThreads()

  .get.deposit_append.EMPT(EMPT,info='feature_WGCNA_cluster_result') <- net
  # print("here 2_6 over")
  EMPT@deposit[['feature_WGCNA_cluster_result']] <-feature_modules
  return(EMPT)
}

.EMP_WGCNA_cluster_analysis_m <- memoise::memoise(.EMP_WGCNA_cluster_analysis)



#' WGCNA cluster analysis
#'
#' @param x Object in EMPT or MultiAssayExperiment format.
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object.
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
#' @param powers a vector of soft thresholding powers for which the scale free topology fit indices are to be calculated.
#' @param RsquaredCut desired minimum scale free topology fitting index R2.
#' @param removeFirst should the first bin be removed from the connectivity histogram?
#' @param nBreaks number of bins in connectivity histograms.
#' @param blockSize block size into which the calculation of connectivity should be broken up. If not given, a suitable value will be calculated using function blockSize and printed if verbose>0. If R runs into memory problems, decrease this value.
#' @param networkType network type. Allowed values are (unique abbreviations of) "unsigned", "signed", "signed hybrid". See WGCNA::adjacency
#' @param moreNetworkConcepts logical: should additional network concepts be calculated? If TRUE, the function will calculate how the network density, the network heterogeneity, and the network centralization depend on the power. For the definition of these additional network concepts, see Horvath and Dong (2008). PloS Comp Biol.
#' @param gcInterval a number specifying in interval (in terms of individual genes) in which garbage collection will be performed. The actual interval will never be less than blockSize.
#' @param TOMType one of "none", "unsigned", "signed", "signed Nowick", "unsigned 2", "signed 2" and "signed Nowick 2". If "none", adjacency will be used for clustering. See WGCNA::TOMsimilarityFromExpr for details.
#' @param minModuleSize minimum module size for module detection. See WCGNA::cutreeDynamic for more details.
#' @param reassignThreshold p-value ratio threshold for reassigning genes between modules.
#' @param mergeCutHeight dendrogram cut height for module merging.
#' @param numericLabels logical: should the returned modules be labeled by colors (FALSE), or by numbers (TRUE)?
#' @param pamRespectsDendro Logical, only used when pamStage is TRUE. If TRUE, the PAM stage will respect the dendrogram in the sense an object can be PAM-assigned only to clusters that lie below it on the branch that the object is merged into. See WGCNA::cutreeDynamic for more details.
#' @param saveTOMs logical: should the consensus topological overlap matrices for each block be saved and returned?
#' @param action A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
# ' @param ... wait_for_add
#'
#' @return EMPT object
#' @export
#'
#' @examples
#' # add example
EMP_WGCNA_cluster_analysis <- function(x,experiment,use_cached=T,powers=c(1:10, seq(from = 12, to=20, by=2)),
                                        RsquaredCut=0.85, removeFirst = FALSE, nBreaks = 10, blockSize = NULL,
                                        # corFnc = WGCNA::cor, corOptions = list(use = 'p'),
                                        networkType = "unsigned",
                                        moreNetworkConcepts = FALSE,
                                        gcInterval = NULL,
                                        TOMType = "unsigned", minModuleSize = 30,
                                        reassignThreshold = 0, mergeCutHeight = 0.25,
                                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                                        # saveTOMs = T,action='add',...) {
                                        saveTOMs = TRUE,action='add') {

  call <- match.call()

  if (inherits(x,"MultiAssayExperiment")) {
    x <- .as.EMPT(x,
                     experiment = experiment)
  }else if(inherits(x,'EMPT')){
    x <- x
    class(x) <- 'EMP_assay_data'
  }else {
    stop('Please check the input data')
  }
  if (use_cached == F) {
    memoise::forget(.EMP_WGCNA_cluster_analysis_m) %>% invisible()
  }
  # print("here 1_1 over")
  EMPT <- .EMP_WGCNA_cluster_analysis_m(EMPT=x,powers=powers,
                                        RsquaredCut=RsquaredCut, removeFirst = removeFirst, nBreaks = nBreaks, blockSize = blockSize,
                                        # corFnc = corFnc, corOptions = corOptions,
                                        networkType = networkType,
                                        moreNetworkConcepts = moreNetworkConcepts,
                                        gcInterval = gcInterval,
                                        TOMType = TOMType, minModuleSize = minModuleSize,
                                        reassignThreshold = reassignThreshold, mergeCutHeight = mergeCutHeight,
                                        numericLabels = numericLabels, pamRespectsDendro = pamRespectsDendro,
                                        # saveTOMs = saveTOMs,...)
                                        saveTOMs = saveTOMs)
  # print("here 1_2 over")                                        

  .get.history.EMPT(EMPT) <- call
  .get.info.EMPT(EMPT) <- 'EMP_WGCNA_cluster_analysis'
  # print("here 1_3 over")
  class(EMPT) <- 'EMP_WGCNA_cluster_analysis'
  if (action == 'add') {
    return(EMPT)
  }else if(action == 'get') {
    return(.get.result.EMPT(EMPT))
  }else{
    warning('action should be one of add or get!')
  }

}

#' @importFrom dplyr where
.EMP_WGCNA_cor_analysis_EMPT <-function(EMPT,method='spearman',coldata_to_assay=NULL){
  var1 <- NULL
  call <- match.call()
  experiment_name <- .get.experiment.EMPT(EMPT)
  #net <- EMPT@deposit_append[['feature_WGCNA_cluster_result']]
  WGCNA_cluster_result <- .get.result.EMPT(EMPT,info='EMP_WGCNA_cluster_analysis')
  net <- WGCNA_cluster_result[['WGCNA_cluster_result']]
  if (is.null(coldata_to_assay)) {
    coldata <- EMPT %>% EMP_coldata_extract() %>%
      tibble::column_to_rownames('primary') %>%
      dplyr::select(where(is.numeric))
  }else{
    coldata <- EMPT %>% EMP_coldata_extract() %>%
      tibble::column_to_rownames('primary') %>%
      dplyr::select(dplyr::all_of(!!coldata_to_assay))
  }

  assay_data <- EMPT %>% .get.assay.EMPT() %>%
    tibble::column_to_rownames('primary') %>%
    as.matrix()

  moduleLabelsAutomatic <- net$colors
  moduleColorsAutomatic <- WGCNA::labels2colors(moduleLabelsAutomatic)
  moduleColorsWW <- moduleColorsAutomatic
  MEs0 <- WGCNA::moduleEigengenes(assay_data, moduleColorsWW)$eigengenes
  MEsWW <- orderMEs(MEs0)

  real_samples <- intersect(rownames(coldata),rownames(MEsWW))
  coldata <- coldata[real_samples,]
  assay_data <- assay_data[real_samples,]

  df.cor.p<-agricolae_correlation(x=coldata,y=MEsWW,method = method)

  df <- df.cor.p$correlation %>%
    as.data.frame() %>%
    tibble::rownames_to_column('var1') %>%
    tidyr::pivot_longer(cols = -var1,names_to = 'var2',values_to = 'coefficient')

  df.p <- df.cor.p$pvalue %>%
    as.data.frame() %>%
    tibble::rownames_to_column('var1') %>%
    tidyr::pivot_longer(cols = -var1,names_to = 'var2',values_to = 'pvalue')
  df$pvalue<-round(df.p$pvalue,2)

  data1_sample_num <- rownames(assay_data) %>% unique %>% length()
  data2_sample_num <- rownames(coldata) %>% unique %>% length()


  df.cor.p[['cor_info']] <- c(experiment_name,paste0(experiment_name,'_codata'))
  df.cor.p[["n.obs"]] <- c(data1_sample_num,data2_sample_num,length(real_samples))
  df.cor.p[['cor_p']] <- df

  EMPT@deposit_append[['WGCNA_cor_result']] <- df.cor.p
  .get.info.EMPT(EMPT) <-'EMP_WGCNA_cor_analysis'
  .get.method.EMPT(EMPT) <- method
  .get.history.EMPT(EMPT) <- call
  class(EMPT) <- 'EMP_WGCNA_cor_analysis'
  return(EMPT)
}

.EMP_WGCNA_cor_analysis_EMPT_m <- memoise::memoise(.EMP_WGCNA_cor_analysis_EMPT)

#' @importFrom WGCNA orderMEs
.EMP_WGCNA_cor_analysis_EMP <- function(EMP,select=NULL,method='spearman'){
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

  real_samples <- intersect(rownames(data1),rownames(data2))
  data1_sample_num <- rownames(data1) %>% unique %>% length()
  data2_sample_num <- rownames(data2) %>% unique %>% length()

  data1 <- dplyr::filter(rownames(data1) %in% real_samples)
  data2 <- dplyr::filter(rownames(data2) %in% real_samples)

  net <- EMP@ExperimentList[[1]]@deposit_append[['feature_WGCNA_cluster_result']]
  if(is.null(net)){
    stop("Experment 1 should run EMP_WGCNA_analysis first!")
  }

  moduleLabelsAutomatic <- net$colors
  moduleColorsAutomatic <- WGCNA::labels2colors(moduleLabelsAutomatic)
  moduleColorsWW <- moduleColorsAutomatic
  MEs0 <- WGCNA::moduleEigengenes(data1, moduleColorsWW)$eigengenes
  MEsWW <- orderMEs(MEs0)

  df.cor.p<-agricolae_correlation(x=data2,y=MEsWW,method = method)

  df <- df.cor.p$correlation %>%
    as.data.frame() %>%
    tibble::rownames_to_column('var1') %>%
    tidyr::pivot_longer(cols = -var1,names_to = 'var2',values_to = 'coefficient')

  df.p <- df.cor.p$pvalue %>%
    as.data.frame() %>%
    tibble::rownames_to_column('var1') %>%
    tidyr::pivot_longer(cols = -var1,names_to = 'var2',values_to = 'pvalue')
  df$pvalue<-round(df.p$pvalue,2)

  data1_sample_num <- rownames(data1) %>% unique %>% length()
  data2_sample_num <- rownames(data2) %>% unique %>% length()


  df.cor.p[['cor_info']] <- cor_info
  df.cor.p[["n.obs"]] <- c(data1_sample_num,data2_sample_num,length(real_samples))
  df.cor.p[['cor_p']] <- df

  EMP@deposit[['WGCNA_cor_analysis_result']] <- df.cor.p
  .get.info.EMP(EMP) <-'EMP_WGCNA_cor_analysis2'
  .get.method.EMP(EMP) <- method
  class(EMP) <- 'EMP_WGCNA_cor_analysis2'
  return(EMP)

}

.EMP_WGCNA_cor_analysis_EMP_m <- memoise::memoise(.EMP_WGCNA_cor_analysis_EMP)
