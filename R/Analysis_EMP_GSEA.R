.Signal2Noise_caculate <- function(EMPT,estimate_group,group_level=NULL) {
  primary <- feature <- value <- mean_value <- sd_value <- min_sd <- real_sd <- NULL
  real_mean <- mean_diff <- sd_sum <- NULL
  coldata <- .get.mapping.EMPT(EMPT) %>% dplyr::select(primary,!!estimate_group)

  ## clean the missing value in the group label
  if(any(is.na(coldata[[estimate_group]]))) {
    stop('Column ',estimate_group,' has beed deteced missing value,please check and filter them!')
  }
  
  data  <- .get.assay.EMPT(EMPT)  %>% tidyr::pivot_longer(cols = -primary,names_to = 'feature',values_to = 'value') %>%
    dplyr::left_join(coldata,by='primary')

  ## in case that data contain more than 2 groups
  group_num <- coldata[[estimate_group]] %>% unique() %>% length()
  if(is.null(group_level)){
    if(group_num != 2){
       stop("The estimate_group parameter must have exactly 2 factors!")
    }
  }else {
    if(length(group_level) != 2){
      stop("The group_level parameter must have exactly 2 factors!")
    }else{
      if(group_level[1] == group_level[2]){
        stop("The group_level parameter must not have same factor!")
      }
      data %<>% dplyr::filter(!!dplyr::sym(estimate_group) %in% !!group_level)
    }
  }

  data_caculate <- data %>% dplyr::group_by(!!dplyr::sym(estimate_group),feature) %>% 
    dplyr::summarise(mean_value=mean(value),sd_value=sd(value),.groups = 'drop') %>%
    dplyr::mutate(min_sd = abs(mean_value) * 0.2) %>%
    dplyr::mutate(real_sd = ifelse(sd_value < min_sd,min_sd,sd_value)) %>%
    dplyr::mutate(real_mean = ifelse(real_sd == min_sd & mean_value == 0,1,mean_value))
  
  
  
  group_level <- .check_group_level.EMPT(EMPT,group_level,estimate_group)
  
  result <- data_caculate %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(mean_diff = real_mean[!!dplyr::sym(estimate_group) == group_level[1]] - real_mean[!!dplyr::sym(estimate_group) == group_level[2]],
                     sd_sum = sum(real_sd)) %>%
    dplyr::mutate(Signal2Noise = mean_diff/sd_sum) %>%
    dplyr::mutate(vs = paste0(group_level[1],' vs ',group_level[2]))
  
  return(result)
}

#' @importFrom clusterProfiler GSEA
gsea_kegg <- function(geneList, kegg.params, pvalueCutoff, seed, ...) {
  keyType <- kegg.params$keyType
  KEGG_Type <- kegg.params$KEGG_Type
  species <- kegg.params$species
  if (is.null(keyType)) {
        stop("keyType should be specified as ko, ec,cpd or entrezid!")
      }else if(!keyType %in% c('ko','ec','cpd','entrezid')){
        stop("keyType should be ko, ec ,cpd or entrezid!")
      }
  
      if(!KEGG_Type %in% c('KEGG','MKEGG')){
        stop("keyType should be KEGG or MKEGG!")
      }
  
    gson_data <- build_gson(keyType = keyType, KEGG_Type = KEGG_Type, species = species)

    enrich.data <- clusterProfiler::GSEA(geneList,gson = gson_data,pvalueCutoff=pvalueCutoff,seed=seed,...) %>% suppressWarnings()
    return(enrich.data)
}

#' @importFrom clusterProfiler gseGO
gsea_go <- function(geneList, go.params, pvalueCutoff, seed, ...) {
  keyType <- go.params$keyType |> toupper()  
  ont <- go.params$ont
  if (is.null(go.params$OrgDb)) {
    stop("Go analysis need OrgDb!")
  }else{
    OrgDb <- go.params$OrgDb
  }
 
  enrich.data <- clusterProfiler::gseGO(geneList, ont = ont, OrgDb = OrgDb,keyType=keyType, pvalueCutoff=pvalueCutoff,seed=seed,...) %>% suppressWarnings()
  return(enrich.data)
}

#' @importFrom ReactomePA gsePathway
gsea_reactome <- function(geneList, reactome.params, pvalueCutoff, seed, ...) {
  organism <- reactome.params$organism
  enrich.data <- ReactomePA::gsePathway(geneList, organism = organism, pvalueCutoff=pvalueCutoff,seed=seed,...) %>% suppressWarnings()
  return(enrich.data)
}

#' @importFrom clusterProfiler enrichWP
gsea_wikipathway <- function(geneList, wikipathway.params, ...) {
  organism <- wikipathway.params$organism
  enrich.data <- clusterProfiler::enrichWP(geneList, organism = organism, ...) %>% suppressWarnings()
  return(enrich.data)
}

#' @importFrom DOSE gseDO
gsea_do <- function(geneList, do.params, pvalueCutoff, seed, ...) {
  if (do.params$ont == 'DOLite') {
    stop("DOLite was removed in the current version.")
  }else{
    ont <- do.params$ont
  }
  organism <- do.params$organism
  enrich.data <- DOSE::gseDO(geneList, organism = organism, pvalueCutoff=pvalueCutoff,seed=seed, ...) %>% suppressWarnings()
  return(enrich.data)
}

get_enrich_data <- function(method = "kegg", geneList,
                            pvalueCutoff=1,seed=T, gson = NULL,
                            TERM2GENE = NULL,
                            TERM2NAME = NA,
                            kegg.params = list(keyType = NULL,
                                               KEGG_Type='KEGG',
                                               species = "all"),
                            go.params = list(OrgDb = NULL,
                                             keyType = "ENTREZID",
                                             ont = "MF"),
                            reactome.params = list(organism = "human"),
                            wikipathway.params = list(organism = "Homo sapiens"),
                            do.params = list(ont = "DO",
                                             organism = "hsa"),
                            ...) {
  if (!is.null(TERM2GENE)) {
    gson <- NULL
    enrich.data <- clusterProfiler::GSEA(geneList, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, pvalueCutoff=pvalueCutoff,seed=seed,...) %>% suppressWarnings()  
    return(enrich.data)  
  }
  if (!is.null(gson)) {
    enrich.data <- clusterProfiler::GSEA(geneList,gson = gson,pvalueCutoff=pvalueCutoff,seed=seed,...) %>% suppressWarnings()  
    return(enrich.data)  
  } 
  if (method == "kegg") {
    if (kegg.params$KEGG_Type == 'KEGG') {
      KEGG_info <- 'KEGG payhway'
    }else if (kegg.params$KEGG_Type == 'MKEGG') {
      KEGG_info <- 'KEGG module'
    }else {
      KEGG_info <- kegg.params$KEGG_Type
    }    
    message("KEGG analysis performed: \nkeyType: ",kegg.params$keyType,'\t KEGG_Type: ',KEGG_info,'\t species: ',kegg.params$species)    
    enrich.data <- gsea_kegg(geneList, kegg.params, pvalueCutoff = pvalueCutoff, seed = seed, ...)
  }
  if (method == "go") {
    message("Go analysis performed: \nkeyType: ",go.params$keyType,'\t ont: ',go.params$ont)    
    enrich.data <- gsea_go(geneList, go.params, pvalueCutoff = pvalueCutoff, seed = seed, ...)
  }
  if (method == "reactome") {
    reactome.params$organism <- match.arg(reactome.params$organism, c("human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly"))
    message("Reactome analysis performed: \norganism: ",reactome.params$organism)    
    enrich.data <- gsea_reactome(geneList, reactome.params, pvalueCutoff = pvalueCutoff, seed = seed, ...)
  }
  if (method == "wikipathway") {
    enrich.data <- gsea_reactome(geneList, wikipathway.params, ...)
  }
  
  if (method == "do") {
    do.params$organism <- match.arg(do.params$organism, c("hsa","mmu"))
    message("DOSE analysis performed: \nont: ",do.params$ont,"\t organism: ",do.params$organism)    
    enrich.data <- gsea_do(geneList, do.params, pvalueCutoff = pvalueCutoff, seed = seed, ...)
  } 
  return(enrich.data)
}

#' @importFrom clusterProfiler GSEA
.EMP_GSEA_analysis_Signal2Noise <- function(EMPT,estimate_group,group_level=NULL, method = "kegg", 
                                            pseudocount=0.0001,pvalueCutoff=1,threshold=NULL,seed=T, gson = NULL,
                                            TERM2GENE = NULL,
                                            TERM2NAME = NA,
                                            kegg.params = list(keyType = NULL,
                                                               KEGG_Type='KEGG',
                                                               species = "all"),
                                            go.params = list(OrgDb = NULL,
                                                             keyType = "ENTREZID",
                                                             ont = "MF"),
                                            reactome.params = list(organism = "human"),
                                            wikipathway.params = list(organism = "Homo sapiens"),
                                            do.params = list(ont = "HDO",
                                                             organism = "hsa"),
                                            ...){
    method <- tolower(method)
    method <- match.arg(method, c("kegg", "go", "reactome", "wikipathway", "do"))
    Signal2Noise <- vs <- NULL
    Signal2Noise_data <- .Signal2Noise_caculate(EMPT,estimate_group,group_level) %>%
                            dplyr::mutate(Signal2Noise = replace(Signal2Noise, Signal2Noise == 0, pseudocount)) %>%
                            tidyr::drop_na() 
                             
    
    if (!is.null(threshold)) {
      Signal2Noise_data %<>% dplyr::filter(Signal2Noise >= threshold)
    }
    
    Signal2Noise_data %<>% dplyr::arrange(dplyr::desc(Signal2Noise)) ## rank the feature
    
    # get the named vector
    geneList <- Signal2Noise_data[['Signal2Noise']]
    geneList <- setNames(geneList,Signal2Noise_data$feature)
    
    enrich.data <- get_enrich_data(method = method, geneList = geneList, 
                            pvalueCutoff=pvalueCutoff,seed=seed, gson = gson,
                            TERM2GENE = TERM2GENE,
                            TERM2NAME = TERM2NAME,
                            kegg.params = kegg.params,
                            go.params = go.params,
                            reactome.params = reactome.params,
                            wikipathway.params = wikipathway.params,
                            do.params = do.params,
                            ...)

    EMPT@deposit[['enrich_data']] <- enrich.data
    .get.estimate_group_info.EMPT(EMPT) <- Signal2Noise_data %>% dplyr::pull(var = vs) %>% unique()
    message('VS info: ',.get.estimate_group_info.EMPT(EMPT))
    message('The Signal2Noise values are arranged in descending order.')
    return(EMPT)

}

#' @importFrom dplyr desc
.EMP_GSEA_analysis_cor <- function(EMPT,estimate_group, method = "kegg", cor_method='pearson',keyType=NULL,KEGG_Type='KEGG', species = "all",
                                    pvalueCutoff=1,threshold_r=0,threshold_p=0.05,seed=T,gson = NULL, 
                                    TERM2GENE = NULL,
                                    TERM2NAME = NA,
                                    kegg.params = list(keyType = NULL,
                                                       KEGG_Type='KEGG',
                                                       species = "all"),
                                    go.params = list(OrgDb = NULL,
                                                     keyType = "ENTREZID",
                                                     ont = "MF"),
                                    reactome.params = list(organism = "human"),
                                    wikipathway.params = list(organism = "Homo sapiens"),
                                    do.params = list(ont = "HDO",
                                                     organism = "hsa"),
                                    ...){
      primary <- NULL
      if(is.null(estimate_group)){
        stop('GSEA based on correlation analysis need estimate_group parameter!')
      }

      feature_table <- .get.assay.EMPT(EMPT) %>% 
        dplyr::arrange(primary) %>% ## confirm the sample rank
        tibble::column_to_rownames('primary')

      coldata <- .get.mapping.EMPT(EMPT) %>% 
        dplyr::arrange(primary) %>% ## confirm the sample rank
        tibble::column_to_rownames('primary') %>%
        dplyr::select(!!estimate_group) 
    
      ## clean the missing value in the group label
      if(any(is.na(coldata[[estimate_group]]))) {
        stop('Column ',estimate_group,' has beed deteced missing value, please check and filter them!')
      }

      #feature_table <- na.omit(feature_table)

      real_samples <- intersect(rownames(coldata),rownames(feature_table))
      coldata <- coldata %>% dplyr::filter(rownames(coldata) %in% real_samples )
      feature_table <- feature_table %>% dplyr::filter(rownames(feature_table) %in% real_samples ) 

      #data.corr <- psych::corr.test(feature_table, coldata,method = cor_method,adjust='none')
      #data.corr <- agricolae_correlation(feature_table, coldata,method = cor_method)
      data.corr <- CorRcpp(x = feature_table,y = coldata,type = cor_method)
      names(data.corr) <- c('correlation','pvalue')

      #data.corr[["correlation"]] <- round(data.corr[["correlation"]],2)
      #data.corr[["pvalue"]] <- round(data.corr[["pvalue"]],2)
      data.r <- data.corr$correlation
      data.p <- data.corr$pvalue
   
      data.r[data.p>threshold_p|abs(data.r)<threshold_r] = 0 ## filter according to the threshold
      
      data.r %<>% as.data.frame() %>% 
        dplyr::filter(!!dplyr::sym(estimate_group) != 0) %>% ## filter the irrelevant feature
        dplyr::arrange(dplyr::desc(!!dplyr::sym(estimate_group))) 
  
    
    # get the named vector
    geneList <- data.r[[1]]
    geneList <- setNames(geneList,rownames(data.r))
    enrich.data <- get_enrich_data(method = method,geneList = geneList, 
                        pvalueCutoff=pvalueCutoff,seed=seed, gson = gson,
                        TERM2GENE = TERM2GENE,
                        TERM2NAME = TERM2NAME,
                        kegg.params = kegg.params,
                        go.params = go.params,
                        reactome.params = reactome.params,
                        wikipathway.params = wikipathway.params,
                        do.params = do.params,
                        ...)
    
    EMPT@deposit[['enrich_data']] <- enrich.data
    .get.estimate_group.EMPT(EMPT) <- estimate_group
    .get.algorithm.EMPT(EMPT) <- 'enrich_analysis'
    .get.info.EMPT(EMPT) <- 'EMP_enrich_analysis'
    message('Correlation analysis based on ',estimate_group)
    message('The correlation analysis values are arranged in descending order.')
    return(EMPT)
}

#' @importFrom dplyr desc
.EMP_GSEA_analysis_log2FC <- function(EMPT,condition,method = "kegg",pvalueCutoff=1,seed=T, gson = NULL, 
                                      TERM2GENE = NULL,
                                      TERM2NAME = NA,
                                      kegg.params = list(keyType = NULL,
                                                         KEGG_Type='KEGG',
                                                         species = "all"),
                                      go.params = list(OrgDb = NULL,
                                                       keyType = "ENTREZID",
                                                       ont = "MF"),
                                      reactome.params = list(organism = "human"),
                                      wikipathway.params = list(organism = "Homo sapiens"),
                                      do.params = list(ont = "HDO",
                                                       organism = "hsa"),
                                      ...){
  log2FC <- NULL
  data <- EMPT@deposit[['diff_analysis_result']] 

  if (is.null(data)) {
    stop('GSEA based on foldchange need EMP_diff_analysis first!')
  }else{
    data %<>% dplyr::filter({{ condition }})
  }
  
  data %<>% tidyr::drop_na() %>%
    dplyr::arrange(desc(log2FC))
  # get the named vector
  geneList <- data[['log2FC']]
  geneList <- setNames(geneList,data$feature)
  enrich.data <- get_enrich_data(method = method,geneList = geneList,
                        pvalueCutoff=pvalueCutoff,seed=seed, gson = gson,
                        TERM2GENE = TERM2GENE,
                        TERM2NAME = TERM2NAME,
                        kegg.params = kegg.params,
                        go.params = go.params,
                        reactome.params = reactome.params,
                        wikipathway.params = wikipathway.params,
                        do.params = do.params,
                        ...)
  
  EMPT@deposit[['enrich_data']] <- enrich.data
  .get.algorithm.EMPT(EMPT) <- 'enrich_analysis'
  .get.info.EMPT(EMPT) <- 'EMP_enrich_analysis'
  message('VS info: ',.get.estimate_group_info.EMPT(EMPT))
  message('The log2FC values are arranged in descending order.')
  return(EMPT)
}

#' Gene set enrichment analysis.
#'
#' @param obj Object in EMPT or MultiAssayExperiment format.
#' @param condition Expressions that return a logical value. The alogarithm condition used in method = "log2FC".eg. pvalue < 0.05
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object.
#' @param estimate_group A character string. Select the column you are interested in the coldata.
#' @param method A character string. Methods include signal2Noise, cor, log2FC.
#' @param enrich_method enrichment method, one of "kegg", "go", "reactome" and "do".
#' @param cor_method A character string including pearson, spearman. The alogarithm cor_method used in method = "cor".
#' @param group_level A series of character strings. Determine the comparison order of groups when method = "log2FC".
#' @param pseudocount A number. The alogarithm pseudocount used in method = "signal2Noise", adjust the 0 in the signal2Noise result into pseudocount value. (default:0.0001)
#' @param pvalueCutoff A character string. Adjusted pvalue cutoff on enrichment tests to report.
#' @param threshold A number. The alogarithm threshold used in method = "signal2Noise",filter out the feature below the signal2Noise threshold.
#' @param threshold_r A number. The alogarithm threshold used in method = "cor",filter out the feature below the abusolte corffcient threshold.
#' @param threshold_p A number. The alogarithm threshold used in method = "cor",filter out the feature above the cor test pavlue threshold.
#' @param seed An interger. Set a random seed to ensure the results are reproducible.
#' @param action A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param gson a GSON object, if not NULL, use it as annotation data.
#' @param TERM2GENE user input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene. Only used when gson is NULL.
#' @param TERM2NAME user input of TERM TO NAME mapping, a data.frame of 2 column with term and name. Only used when gson is NULL.
#' @param keyType For KEGG analysis, keyType include ko, ec, cpd, entrezid. For Go analysis, keyType include entrezid and symbol.
#' @param KEGG_Type A character string. KEGG_Type include KEGG and MKEGG in KEGG analysis. KEGG means KEGG pathway. MKEGG means KEGG module.
#' @param species A character string. Species includ all, hsa, mmu,...in in KEGG analysis. Supported organism listed in 'https://www.genome.jp/kegg/catalog/
#' @param OrgDb OrgDb in Go analysis.
#' @param ont For Go analysis, ont include "BP", "MF","CC", and "ALL". For DOSE analysis, ont only support "DO".
#' @param organism For Reactome analysis, organism include "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly". For DOSE analysis, organism include "hsa" and "mmu".
#' @param ... Further parameters passed to clusterProfiler::GSEA, clusterProfiler::gseGO, ReactomePA::gsePathway and DOSE::gseDO
#'
#' @return EMPT object
#' @export
#' @section Detaild about method:
#' The EMP_GSEA_analysis moudle performed based on cluserProfiler, more detailed information are available on:
#'
#' http://yulab-smu.top/biomedical-knowledge-mining-book/index.html
#' @examples
#' \dontrun{
#' data(MAE)
#' 
#' # KEGG
#' ## based on cor analysis
#' MAE |> EMP_assay_extract('geno_ko') |>
#'   EMP_GSEA_analysis(experiment = 'geno_ko',method='cor',
#'                     enrich_method = 'kegg',
#'                     keyType='ko',
#'                     estimate_group = 'BMI',cor_method = 'spearman',
#'                     threshold_r = 0.3,threshold_p = 0.05, ## filter by coe and pvalue
#'                     pvalueCutoff = 0.05)
#' 
#' ## based on diff analysis
#' MAE |>
#'   EMP_diff_analysis(experiment = 'geno_ec',method='DESeq2',.formula = ~0+Group,
#'                     group_level=c('Group_A','Group_B')) |>
#'   EMP_GSEA_analysis(method='log2FC',pvalue<0.05,enrich_method = 'kegg',
#'                     keyType = 'ec',KEGG_Type = 'KEGG',pvalueCutoff = 0.05)
#' 
#' ## based on signal2Noise
#' MAE |>
#'   EMP_GSEA_analysis(experiment = 'geno_ko',method='signal2Noise',
#'                     estimate_group = 'Group',enrich_method = 'kegg',
#'                     pvalueCutoff = 0.05,keyType = 'ko')
#' 
#' # GO
#' library(org.Hs.eg.db)
#' MAE |> 
#'   EMP_assay_extract('host_gene') |>
#'   EMP_GSEA_analysis(method='signal2Noise',ont = 'MF',keyType = 'symbol',
#'                     organism='hsa',OrgDb = org.Hs.eg.db,
#'                     estimate_group = 'Group',enrich_method = 'go',
#'                     pvalueCutoff = 1) 
#' 
#' MAE |> 
#'   EMP_assay_extract('host_gene') |>
#'   EMP_GSEA_analysis(method='cor',ont = 'MF',keyType = 'symbol',
#'                     OrgDb = org.Hs.eg.db,estimate_group = 'BMI',
#'                     cor_method = 'spearman',
#'                     enrich_method = 'go') 
#' 
#' # reactome
#' MAE |> 
#'   EMP_assay_extract('host_gene') |>
#'   EMP_feature_convert(from = 'SYMBOL',to='ENTREZID',species = 'Human') |>
#'   EMP_diff_analysis(method = 'DESeq2',.formula = ~Group) |>
#'   EMP_GSEA_analysis(method='log2FC',
#'                     organism='human',
#'                     enrich_method = 'reactome') 
#' 
#' # DO
#' MAE |> 
#'   EMP_assay_extract('host_gene') |>
#'   EMP_feature_convert(from = 'SYMBOL',to='ENTREZID',species = 'Human') |>
#'   EMP_diff_analysis(method = 'DESeq2',.formula = ~Group) |>
#'   EMP_GSEA_analysis(method='log2FC',
#'                     organism='hsa',
#'                     enrich_method = 'do') 
#' 
#' ## Visualization
#' MAE |>
#'   EMP_GSEA_analysis(experiment = 'geno_ko',method='signal2Noise',
#'                     estimate_group = 'Group',
#'                     pvalueCutoff = 0.05,keyType = 'ko') |>
#'   EMP_curveplot(geneSetID='map00680')
#'   
#' MAE |>
#'   EMP_GSEA_analysis(experiment = 'geno_ko',method='signal2Noise',
#'                     estimate_group = 'Group',
#'                     pvalueCutoff = 0.05,keyType = 'ko') |>
#'   EMP_dotplot(color='p.adjust',showCategory=10) 
#' 
#' MAE |>
#'   EMP_GSEA_analysis(experiment = 'geno_ko',method='signal2Noise',
#'                     estimate_group = 'Group',
#'                     pvalueCutoff = 0.05,keyType = 'ko') |>
#'   EMP_netplot(showCategory=5) 
#' }
EMP_GSEA_analysis <- function(obj,condition,experiment,estimate_group=NULL,method=NULL, enrich_method = "kegg", cor_method='pearson',group_level=NULL,
                               pseudocount=0.0001,pvalueCutoff=1,threshold=NULL,
                               threshold_r=0,threshold_p=0.05,seed=TRUE,action='add', gson = NULL, 
                               TERM2GENE = NULL,
                               TERM2NAME = NA,
                               KEGG_Type='KEGG',species = "all", # kegg.params
                               OrgDb = NULL, keyType = "entrezid", ont = if (enrich_method == "go") "MF" else "DO", # go.params
                               organism = if (method == "reactome") "human" else "Homo sapiens", # wikipathway.params
                               ...){
  kegg.params = list(keyType = keyType,
                   KEGG_Type = KEGG_Type,
                   species = species)
  go.params = list(OrgDb = OrgDb,
                   keyType = keyType,
                   ont = ont)
  reactome.params = list(organism = organism)
  wikipathway.params = list(organism = organism)
  do.params = list(ont = ont,
                   organism = organism)
  rlang::check_installed(c('BiocManager'), reason = 'for EMP_GSEA_analysis().', action = install.packages) 
  rlang::check_installed(c('clusterProfiler'), reason = 'for EMP_GSEA_analysis().', action = BiocManager::install)    
  
  call <- match.call()
  if (is(obj,"MultiAssayExperiment")) {
    
    EMPT <- .as.EMPT(obj,
                     experiment = experiment)
  }else if(is(obj,'EMPT')) {
    EMPT <- obj
  }
  
  switch(method,
         "signal2Noise" = {  
           EMPT %<>%  .EMP_GSEA_analysis_Signal2Noise(estimate_group,group_level,method = enrich_method, 
                                                      pseudocount,pvalueCutoff,threshold,seed,gson=gson,
                                                      TERM2GENE = TERM2GENE,
                                                      TERM2NAME = TERM2NAME,
                                                      kegg.params = kegg.params,
                                                      go.params = go.params,
                                                      reactome.params = reactome.params,
                                                      wikipathway.params = wikipathway.params,
                                                      do.params = do.params,
                                                      ...)
         },
         "cor" = {  
           if (!cor_method %in% c('pearson','spearman')) {
             stop('Parameter method in EMP_GSEA_analysis_cor should be one of pearson, spearman! ')
           }
           EMPT %<>%  .EMP_GSEA_analysis_cor(estimate_group,cor_method = cor_method,
                                             method = enrich_method, 
                                             pvalueCutoff = pvalueCutoff,threshold_r,threshold_p,seed,gson=gson,
                                             TERM2GENE = TERM2GENE,
                                             TERM2NAME = TERM2NAME,
                                             kegg.params = kegg.params,
                                             go.params = go.params,
                                             reactome.params = reactome.params,
                                             wikipathway.params = wikipathway.params,
                                             do.params = do.params,
                                             ...)
         },
         "log2FC" = {    
           EMPT %<>% .EMP_GSEA_analysis_log2FC({{condition}},method = enrich_method,pvalueCutoff = pvalueCutoff,seed = seed,gson=gson,                                           
                                                      TERM2GENE = TERM2GENE,
                                                      TERM2NAME = TERM2NAME,
                                                      kegg.params = kegg.params,
                                                      go.params = go.params,
                                                      reactome.params = reactome.params,
                                                      wikipathway.params = wikipathway.params,
                                                      do.params = do.params,
                                                      ...)
         },
         {
           stop('Parameter method in EMP_GSEA_analysis must be one of signal2Noise,cor,log2FC!')
         }
  )
  .get.history.EMPT(EMPT) <- call
  class(EMPT) <- 'EMP_enrich_analysis'
  .get.method.EMPT(EMPT) <- method
  .get.algorithm.EMPT(EMPT) <- 'enrich_analysis'
  .get.info.EMPT(EMPT) <- 'EMP_enrich_analysis'
  if (action == 'add') {
    return(EMPT)
  }else if(action == 'get') {
    return(.get.result.EMPT(EMPT))
  }else{
    warning('action should be one of add or get!')
  }  
}









