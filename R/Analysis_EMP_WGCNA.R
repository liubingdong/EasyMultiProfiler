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
                                        saveTOMs = FALSE) {
  deposit <- list()
  Var1 <- Freq <- NULL
  #enableWGCNAThreads()
  
  assay_data <- assay(EMPT) %>% t()

  sft <- spsUtil::quiet(WGCNA::pickSoftThreshold(assay_data, powerVector = powers,
                          RsquaredCut=RsquaredCut, removeFirst = removeFirst, nBreaks = nBreaks, blockSize = blockSize,
                          # corFnc = corFnc, corOptions = corOptions,
                          networkType = networkType,
                          moreNetworkConcepts = moreNetworkConcepts,
                          gcInterval = gcInterval),print_cat = TRUE, message = TRUE, warning = TRUE)

  check_best_power <- sft[['powerEstimate']]

  if (is.na(check_best_power)) {
    stop("No best soft power,please reset the RsquaredCut or other parameter!")
  }

  # WGCNA_blockwiseModules can only be applied to a 'numeric', not a 'integer'
  if (is.integer(assay_data)) {
    assay_data <- apply(assay_data, 2, as.numeric)
  }

  #save(assay_data, check_best_power, TOMType, minModuleSize, reassignThreshold, mergeCutHeight, numericLabels, pamRespectsDendro, saveTOMs, file = "test.Rdata")
  # net = WGCNA::blockwiseModules(assay_data, power = check_best_power,
  net <- WGCNA_blockwiseModules(assay_data, power = check_best_power,
                         TOMType = TOMType, minModuleSize = minModuleSize,
                         reassignThreshold = reassignThreshold, mergeCutHeight = mergeCutHeight,
                         numericLabels = numericLabels, pamRespectsDendro = pamRespectsDendro,
                         # saveTOMs = saveTOMs,...)
                         saveTOMs = saveTOMs)

  mergedColors = WGCNA::labels2colors(net$colors)

  WGCNA::plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      abHeight=mergeCutHeight)

  mergedColors = WGCNA::labels2colors(net$colors)

  WGCNA_module_elements <- table(WGCNA::labels2colors(net$colors)) %>%
                as.data.frame() %>%
                dplyr::rename(WGCNA_color=Var1,WGCNA_module_elements=Freq)

  feature_modules <- data.frame(WGCNA_cluster=net$unmergedColors,WGCNA_color=mergedColors) %>% tibble::rownames_to_column('feature') %>%
    tibble::as_tibble() %>% dplyr::left_join(WGCNA_module_elements,by='WGCNA_color')


  #disableWGCNAThreads()

  .get.deposit_append.EMPT(EMPT,info='feature_WGCNA_cluster_result') <- net
  # print("here 2_6 over")
  EMPT@deposit[['feature_WGCNA_cluster_result']] <-feature_modules
  return(EMPT)
}

.EMP_WGCNA_cluster_analysis_m <- memoise::memoise(.EMP_WGCNA_cluster_analysis,cache = cachem::cache_mem(max_size = 2048 * 1024^2))



#' WGCNA cluster analysis
#'
#' @param obj Object in EMPT or MultiAssayExperiment format.
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
#'
#' @return EMPT object
#' @export
#'
#' @examples
#' \dontrun{
#' data(MAE)
#' MAE |>
#'   EMP_assay_extract('geno_ec') |>
#'   EMP_WGCNA_cluster_analysis(RsquaredCut = 0.85)
#' 
#' MAE |>
#'   EMP_assay_extract('geno_ko') |>
#'   EMP_WGCNA_cluster_analysis(RsquaredCut = 0.8,mergeCutHeight=0.4)
#' }
EMP_WGCNA_cluster_analysis <- function(obj,experiment,use_cached=T,powers=c(1:10, seq(from = 12, to=20, by=2)),
                                        RsquaredCut=0.85, removeFirst = FALSE, nBreaks = 10, blockSize = NULL,
                                        # corFnc = WGCNA::cor, corOptions = list(use = 'p'),
                                        networkType = "unsigned",
                                        moreNetworkConcepts = FALSE,
                                        gcInterval = NULL,
                                        TOMType = "unsigned", minModuleSize = 30,
                                        reassignThreshold = 0, mergeCutHeight = 0.25,
                                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                                        # saveTOMs = T,action='add',...) {
                                        saveTOMs = FALSE,action='add') {

  call <- match.call()

  if (inherits(obj,"MultiAssayExperiment")) {
    x <- .as.EMPT(obj,
                     experiment = experiment)
  }else if(inherits(obj,'EMPT')){
    x <- obj
    class(x) <- 'EMP_assay_data'
  }else {
    stop('Please check the input data for EMP_WGCNA_cluster_analysis!')
  }
  if (use_cached == FALSE) {
    memoise::forget(.EMP_WGCNA_cluster_analysis_m) %>% invisible()
  }

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
.EMP_WGCNA_cor_analysis_EMPT <-function(obj,method='spearman',coldata_to_assay=NULL,action='add'){
  var1 <- NULL
  
  if (inherits(obj,"EMPT")) {
    EMPT <- obj
  }else{
    stop('Please check the input data!')
  }  

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

  assay_data <- assay(EMPT) %>% t()

  moduleLabelsAutomatic <- net$colors
  moduleColorsAutomatic <- WGCNA::labels2colors(moduleLabelsAutomatic)
  moduleColorsWW <- moduleColorsAutomatic
  MEs0 <- WGCNA::moduleEigengenes(assay_data, moduleColorsWW)$eigengenes
  MEsWW <- WGCNA::orderMEs(MEs0)

  # filter the samnples with miss value
  coldata <- na.omit(coldata)
  MEsWW <- na.omit(MEsWW)

  real_samples <- intersect(rownames(coldata),rownames(MEsWW))
  coldata <- coldata %>% dplyr::filter(rownames(coldata) %in% real_samples )
  MEsWW <- MEsWW %>% dplyr::filter(rownames(MEsWW) %in% real_samples ) 

  #df.cor.p<-agricolae_correlation(x=coldata,y=MEsWW,method = method)
  df.cor.p <- CorRcpp(x = coldata,y = MEsWW,type = method)
  names(df.cor.p) <- c('correlation','pvalue')

  df.cor.p[["correlation"]] <- round(df.cor.p[["correlation"]],2)
  df.cor.p[["pvalue"]] <- round(df.cor.p[["pvalue"]],2)


  df <- df.cor.p$correlation %>%
    as.data.frame() %>%
    tibble::rownames_to_column('var1') %>%
    tidyr::pivot_longer(cols = -var1,names_to = 'var2',values_to = 'coefficient')

  df.p <- df.cor.p$pvalue %>%
    as.data.frame() %>%
    tibble::rownames_to_column('var1') %>%
    tidyr::pivot_longer(cols = -var1,names_to = 'var2',values_to = 'pvalue')
 
  df$pvalue <- df.p$pvalue

  data1_sample_num <- rownames(assay_data) %>% unique %>% length()
  data2_sample_num <- rownames(coldata) %>% unique %>% length()


  df.cor.p[['cor_info']] <- c(experiment_name,paste0(experiment_name,'_coldata'))
  df.cor.p[["n.obs"]] <- c(data1_sample_num,data2_sample_num,length(real_samples))
  df.cor.p[['cor_p']] <- df
  df.cor.p[['MEsWW']] <- MEsWW
  
  if (action == 'add') {
    EMPT@deposit_append[['WGCNA_cor_result']] <- df.cor.p
    .get.info.EMPT(EMPT) <-'EMP_WGCNA_cor_analysis'
    .get.method.EMPT(EMPT) <- method
    class(EMPT) <- 'EMP_WGCNA_cor_analysis'
    return(EMPT)
  }else if(action == 'get') {
    return(df.cor.p)
  }else{
    warning('action should be one of add or get!')
  }  
}

.EMP_WGCNA_cor_analysis_EMPT_m <- memoise::memoise(.EMP_WGCNA_cor_analysis_EMPT,cache = cachem::cache_mem(max_size = 2048 * 1024^2))

#' @importFrom WGCNA orderMEs
.EMP_WGCNA_cor_analysis_EMP <- function(obj,select=NULL,method='spearman',action='add'){
  primary <- var1 <- NULL
  if (inherits(obj,"EMP")) {
    EMP <- obj
  }else{
    stop('Please check the input data!')
  }  

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

  real_samples <- intersect(rownames(data1),rownames(data2))
  data1_sample_num <- rownames(data1) %>% unique %>% length()
  data2_sample_num <- rownames(data2) %>% unique %>% length()

  data1 <- data1 %>% dplyr::filter(rownames(data1) %in% real_samples)
  data2 <- data2 %>% dplyr::filter(rownames(data2) %in% real_samples)

  net <- EMP@ExperimentList[[1]]@deposit_append[['feature_WGCNA_cluster_result']]
  if(is.null(net)){
    stop("Experment 1 should run EMP_WGCNA_analysis first!")
  }

  moduleLabelsAutomatic <- net$colors
  moduleColorsAutomatic <- WGCNA::labels2colors(moduleLabelsAutomatic)
  moduleColorsWW <- moduleColorsAutomatic
  MEs0 <- WGCNA::moduleEigengenes(data1, moduleColorsWW)$eigengenes
  MEsWW <- orderMEs(MEs0)

  #df.cor.p<-agricolae_correlation(x=data2,y=MEsWW,method = method)

  # filter the samnples with miss value
  data2 <- na.omit(data2)
  MEsWW <- na.omit(MEsWW)

  real_samples <- intersect(rownames(data2),rownames(MEsWW))
  data2 <- data2[real_samples,]
  MEsWW <- MEsWW[real_samples,]


  df.cor.p <- CorRcpp(x = data2,y = MEsWW,type = method)
  names(df.cor.p) <- c('correlation','pvalue')

  df.cor.p[["correlation"]] <- round(df.cor.p[["correlation"]],2)
  df.cor.p[["pvalue"]] <- round(df.cor.p[["pvalue"]],2)


  df <- df.cor.p$correlation %>%
    as.data.frame() %>%
    tibble::rownames_to_column('var1') %>%
    tidyr::pivot_longer(cols = -var1,names_to = 'var2',values_to = 'coefficient')

  df.p <- df.cor.p$pvalue %>%
    as.data.frame() %>%
    tibble::rownames_to_column('var1') %>%
    tidyr::pivot_longer(cols = -var1,names_to = 'var2',values_to = 'pvalue')
  
  df$pvalue <- df.p$pvalue

  data1_sample_num <- rownames(data1) %>% unique %>% length()
  data2_sample_num <- rownames(data2) %>% unique %>% length()


  df.cor.p[['cor_info']] <- cor_info
  df.cor.p[["n.obs"]] <- c(data1_sample_num,data2_sample_num,length(real_samples))
  df.cor.p[['cor_p']] <- df
  df.cor.p[['MEsWW']] <- MEsWW

  if (action == 'add') {
    EMP@deposit[['WGCNA_cor_analysis_result']] <- df.cor.p
    .get.info.EMP(EMP) <-'EMP_WGCNA_cor_analysis2'
    .get.method.EMP(EMP) <- method
    class(EMP) <- 'EMP_WGCNA_cor_analysis2'
    return(EMP)
  }else if(action == 'get') {
    return(df.cor.p)
  }else{
    warning('action should be one of add or get!')
  }  
}

.EMP_WGCNA_cor_analysis_EMP_m <- memoise::memoise(.EMP_WGCNA_cor_analysis_EMP,cache = cachem::cache_mem(max_size = 2048 * 1024^2))


#' EMP_WGCNA_cor_analysis
#'
#' @param obj EMPT or EMP object.
#' @param select A character string. The experiment name in the EMP object. Only for the EMP object.
#' @param method A character string. Methods include pearson (default), spearman.
#' @param coldata_to_assay A series of character strings. Select the column from coldata to caculate. Only for the EMPT object.
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
#' @param action A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @return EMP object
#' @export
#'
#' @examples
#' \dontrun{
#' data(MAE)
#' ## from one experiment
#' WGCNA_COR_result <- MAE |>
#'   EMP_assay_extract('geno_ec')  |> 
#'   EMP_identify_assay(method = 'edgeR',estimate_group = 'Group') |>
#'   EMP_WGCNA_cluster_analysis(RsquaredCut = 0.85,mergeCutHeight=0.4)  |>
#'   EMP_WGCNA_cor_analysis(coldata_to_assay = c('BMI','PHQ9','GAD7','HAMD','SAS','SDS'),
#'                          method='spearman',action='add') # If want the detailed result, set action = 'get'
#' 
#' ## Visualization
#' MAE |>
#'   EMP_assay_extract('geno_ec')  |> 
#'   EMP_identify_assay(method = 'edgeR',estimate_group = 'Group') |>
#'   EMP_WGCNA_cluster_analysis(RsquaredCut = 0.85,mergeCutHeight=0.4)  |>
#'   EMP_WGCNA_cor_analysis(coldata_to_assay = c('BMI','PHQ9','GAD7','HAMD','SAS','SDS'),method='spearman') |>
#'   EMP_heatmap_plot(palette = 'Spectral')
#' 
#' ## Filter the interesting module and make the enrichment analysis
#' MAE |>
#'   EMP_assay_extract('geno_ec')  |> 
#'   EMP_identify_assay(method = 'edgeR',estimate_group = 'Group') |>
#'   EMP_WGCNA_cluster_analysis(RsquaredCut = 0.85,mergeCutHeight=0.4)  |>
#'   EMP_WGCNA_cor_analysis(coldata_to_assay = c('BMI','PHQ9','GAD7','HAMD','SAS','SDS'),method='spearman') |>
#'   EMP_heatmap_plot(palette = 'Spectral') |>
#'   EMP_filter(feature_condition = WGCNA_color == 'brown' ) |> 
#'   EMP_diff_analysis(method = 'DESeq2',.formula = ~Group) |>
#'   EMP_enrich_analysis(keyType = 'ec',KEGG_Type = 'MKEGG') |>
#'   EMP_dotplot()
#' 
#' ## from two different experiments
#' k1 <- MAE |>
#'   EMP_assay_extract('geno_ec')  |> 
#'   EMP_identify_assay(method = 'edgeR',estimate_group = 'Group') |>
#'   EMP_WGCNA_cluster_analysis(RsquaredCut = 0.85,mergeCutHeight=0.4)
#' 
#' k2 <- MAE |>
#'   EMP_assay_extract('host_gene',pattern = c('A1BG','A1CF','A2MP1','AACS'),pattern_ref = 'feature')
#' 
#' (k1 + k2) |>
#'   EMP_WGCNA_cor_analysis(method='spearman') |>
#'   EMP_heatmap_plot(palette = 'Spectral') 
#' }

EMP_WGCNA_cor_analysis <- function(obj,select=NULL,method='spearman',coldata_to_assay=NULL,use_cached=TRUE,action='add'){
  deposit <- NULL
  call <- match.call()
  if (inherits(obj,"EMP")) {
    if (use_cached == FALSE) {
      memoise::forget(.EMP_WGCNA_cor_analysis_EMP_m) %>% invisible()
    }       
    deposit <- .EMP_WGCNA_cor_analysis_EMP_m(obj=obj,select=select,method=method,action=action)
    if (action=='add') {
      .get.history.EMP(deposit) <- call
    }
  }else if (inherits(obj,"EMPT")) {
    if (use_cached == FALSE) {
      memoise::forget(.EMP_WGCNA_cor_analysis_EMPT_m) %>% invisible()
    }       
    deposit <- .EMP_WGCNA_cor_analysis_EMPT_m(obj=obj,coldata_to_assay=coldata_to_assay,method=method,action=action)
    if (action=='add') {
      .get.history.EMPT(deposit) <- call    
    }    
  }else{
    stop("Please check the input data for EMP_WGCNA_cor_analysis!")
  } 
  return(deposit)
}



