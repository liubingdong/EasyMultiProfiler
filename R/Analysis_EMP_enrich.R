#' @importFrom clusterProfiler enricher
enrich_kegg <- function(df, feature_name, kegg.params, minGSSize, maxGSSize, combineGroup, ...) {
  keyType <- kegg.params$keyType
  KEGG_Type <- kegg.params$KEGG_Type
  species <- kegg.params$species
  if (is.null(keyType)) {
    stop("keyType should be specified as ko, ec, cpd or entrezid!")
  }else if(!keyType %in% c('ko','ec','cpd','entrezid')){
    stop("keyType should be ko, ec, cpd or entrezid!")
  }

  if(!KEGG_Type %in% c('KEGG','MKEGG')){
    stop("keyType should be KEGG or MKEGG!")
  }
  gson_data <- build_gson(keyType = keyType, KEGG_Type = KEGG_Type, species = species) |> suppressMessages()

  EMP_message(paste0('KEGG database version: ',gson_data@version),color=32,order=1,show='message')
  EMP_message(paste0('Species: ',gson_data@species),color=32,order=1,show='message')

  if (combineGroup == TRUE) {
    enrich.data <- clusterProfiler::enricher(gene=feature_name, gson=gson_data,minGSSize=minGSSize, maxGSSize=maxGSSize,...) 
  }else{
    enrich.data <- clusterProfiler::compareCluster(feature~sign_group, data=df, 
          fun=clusterProfiler::enricher, gson=gson_data, 
          minGSSize=minGSSize, maxGSSize=maxGSSize,...)     
  }
  return(enrich.data)
}

#' @importFrom clusterProfiler enrichGO
enrich_go <- function(df, feature_name, go.params, minGSSize, maxGSSize, combineGroup, ...) {
  if (is.null(go.params$OrgDb)) {
    stop("Go analysis need OrgDb!")
  }else{
    OrgDb <- go.params$OrgDb
  }
  keyType <- go.params$keyType |> toupper()
  ont <- go.params$ont
  if (combineGroup == TRUE) {
    enrich.data <- clusterProfiler::enrichGO(gene=feature_name, OrgDb=OrgDb, keyType = keyType, ont = ont,
            minGSSize=minGSSize, maxGSSize=maxGSSize,...) 
  }else{
    enrich.data <- clusterProfiler::compareCluster(feature~sign_group, data=df, 
            fun=clusterProfiler::enrichGO, OrgDb=OrgDb, keyType = keyType, ont = ont,
            minGSSize=minGSSize, maxGSSize=maxGSSize,...)    
  }

  return(enrich.data)   
}

#' @importFrom ReactomePA enrichPathway
enrich_reactome <- function(df, feature_name, reactome.params, minGSSize, maxGSSize, combineGroup, ...) {
  organism <- reactome.params$organism
  if (combineGroup == TRUE) {
    enrich.data <- ReactomePA::enrichPathway(gene=feature_name, organism = organism,
          minGSSize=minGSSize, maxGSSize=maxGSSize,...) 
  }else{
    enrich.data <- clusterProfiler::compareCluster(feature~sign_group, data=df, 
            fun=ReactomePA::enrichPathway, organism = organism,
            minGSSize=minGSSize, maxGSSize=maxGSSize,...)    
  }


  return(enrich.data)   
}

enrich_wikipathway <- function(df, feature_name, wikipathway.params, combineGroup, ...) {
  organism <- wikipathway.params$organism
  if (combineGroup == TRUE) {
    enrich.data <- clusterProfiler::enrichWP(gene=feature_name, organism = organism,...) 
  }else{
    enrich.data <- clusterProfiler::enrichWP(feature~sign_group, data=df, 
            fun=ReactomePA::enrichPathway, organism = organism,...)    
  }


  return(enrich.data) 
}

#' @importFrom DOSE enrichDO
enrich_do <- function(df, feature_name, do.params, minGSSize, maxGSSize, combineGroup, ...) {
  ont <- do.params$ont
  organism <- do.params$organism
  if (combineGroup == TRUE) {
    enrich.data <- DOSE::enrichDO(gene=feature_name, ont = ont, organism = organism,
            minGSSize=minGSSize, maxGSSize=maxGSSize,...) 
  }else{
    enrich.data <- clusterProfiler::compareCluster(feature~sign_group, data=df, 
            fun=DOSE::enrichDO, ont = ont, organism = organism,
            minGSSize=minGSSize, maxGSSize=maxGSSize,...)    
  }


  return(enrich.data)  
}


#' @importFrom clusterProfiler enricher
#' @importFrom clusterProfiler compareCluster
.EMP_enrich_analysis <- function(EMPT,condition,minGSSize =1,maxGSSize =500, method = "kegg", 
  combineGroup=FALSE, gson = NULL, TERM2GENE = NULL,
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
  ...){
  method <- tolower(method)
  method <- match.arg(method, c("kegg", "go", "reactome", "do")) ## del wikipathway 
  deposit <- list()
  sign_group <- feature_name <- NULL
  diff_df <- .get.result.EMPT(EMPT,info='EMP_diff_analysis') %>% suppressMessages()
  
  if (is.null(diff_df)) {
    if (combineGroup == FALSE) {
      EMP_message("Without EMP_diff_analysis result, paramter combineGroup could not be FALSE!",color = 31,order = 1,show='message')
    }
    combineGroup <- TRUE
    feature_name <- names(EMPT)
  }else{
    condition <- dplyr::enquo(condition)
    df <- diff_df %>%
      dplyr::filter(!!condition) %>%
      tidyr::drop_na(sign_group) ## filter NA or the result will add NA group!
    feature_name <- df$feature
  }
  
  if (!is.null(TERM2GENE)) {
    gson <- NULL
    if (combineGroup == TRUE) {
      enrich.data <- clusterProfiler::enricher(gene=feature_name, TERM2GENE = TERM2GENE,
        TERM2NAME = TERM2NAME, minGSSize=minGSSize, maxGSSize=maxGSSize,...) 
    }else{
      enrich.data <- clusterProfiler::compareCluster(feature~sign_group, data=df, 
            fun=clusterProfiler::enricher, TERM2GENE = TERM2GENE,
            TERM2NAME = TERM2NAME, minGSSize=minGSSize, maxGSSize=maxGSSize,...)     
    }
    EMPT@deposit[['enrich_data']] <- enrich.data
    .get.algorithm.EMPT(EMPT) <- 'enrich_analysis'
    .get.info.EMPT(EMPT) <- 'EMP_enrich_analysis' 
    return(EMPT) 
  }

  if (!is.null(gson)) {
    if (combineGroup == TRUE) {
    enrich.data <- clusterProfiler::enricher(gene=feature_name, gson=gson,minGSSize=minGSSize, maxGSSize=maxGSSize,...) 
    } else {
      enrich.data <- clusterProfiler::compareCluster(feature~sign_group, data=df, 
            fun=clusterProfiler::enricher, gson=gson, 
            minGSSize=minGSSize, maxGSSize=maxGSSize,...)     
    }
    EMPT@deposit[['enrich_data']] <- enrich.data
    .get.algorithm.EMPT(EMPT) <- 'enrich_analysis'
    .get.info.EMPT(EMPT) <- 'EMP_enrich_analysis' 
    return(EMPT) 
  }

  if (method == "kegg") {
    if (kegg.params$KEGG_Type == 'KEGG') {
      KEGG_info <- 'KEGG payhway'
    }else if (kegg.params$KEGG_Type == 'MKEGG') {
      KEGG_info <- 'KEGG module'
    }else {
      KEGG_info <- kegg.params$KEGG_Type
    }
    info_ouput <- paste0("KEGG analysis performed: \nkeyType: ",kegg.params$keyType,'\t KEGG_Type: ',KEGG_info,'\t species: ',kegg.params$species)
    EMP_message(info_ouput,color=32,order=1,show='message')
    enrich.data <- enrich_kegg(df, feature_name = feature_name, kegg.params = kegg.params, minGSSize = minGSSize, maxGSSize = maxGSSize, combineGroup = combineGroup, ...)
  } 
  
  if (method == "go") {
    info_ouput <- paste0("Go analysis performed: \nkeyType: ",go.params$keyType,'\t ont: ',go.params$ont)
    EMP_message(info_ouput,color=32,order=1,show='message')
    enrich.data <- enrich_go(df, feature_name = feature_name, go.params = go.params, minGSSize = minGSSize, maxGSSize = maxGSSize, combineGroup = combineGroup, ...)
  } 
  if (method == "reactome") {
    rlang::check_installed(c('BiocManager'), reason = 'for enrich_reactome().', action = install.packages)  
    rlang::check_installed(c('ReactomePA'), reason = 'for enrich_reactome().', action = BiocManager::install)
    reactome.params$organism <- match.arg(reactome.params$organism, c("human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly"))
    info_ouput <- paste0("Reactome analysis performed: \norganism: ",reactome.params$organism)
    EMP_message(info_ouput,color=32,order=1,show='message')
    enrich.data <- enrich_reactome(df, feature_name = feature_name, reactome.params = reactome.params, minGSSize = minGSSize, maxGSSize = maxGSSize, combineGroup = combineGroup, ...)
  } 
  if (method == "wikipathway") {
    enrich.data <- enrich_wikipathway(df, feature_name = feature_name, wikipathway.params = wikipathway.params,combineGroup = combineGroup, ...)
  } 
  if (method == "do") {
    if (do.params$ont == 'DOLite') {
      stop("DOLite was removed in the current version.")
    }
    do.params$organism <- match.arg(do.params$organism, c("hsa","mmu"))
    info_ouput <- paste0("DOSE analysis performed: \nont: ",do.params$ont,"\t organism: ",do.params$organism)
    EMP_message(info_ouput,color=32,order=1,show='message')    
    enrich.data <- enrich_do(df, feature_name = feature_name, do.params = do.params, minGSSize=minGSSize, maxGSSize=maxGSSize, combineGroup = combineGroup, ...)
  }
   EMPT@deposit[['enrich_data']] <- enrich.data
  .get.algorithm.EMPT(EMPT) <- 'enrich_analysis'
  .get.info.EMPT(EMPT) <- 'EMP_enrich_analysis' 
  return(EMPT)  
}

#' Enrichment for EMPT object
#'
#' @param obj Object in EMPT or MultiAssayExperiment format.
#' @param condition Expressions that return a logical value according to the result of EMP_diff_analysis. eg. pvalue < 0.05
#' @param method enrichment method, one of "kegg", "go", "reactome" and "do".
#' @param TERM2GENE user input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene. Only used when gson is NULL.
#' @param TERM2NAME user input of TERM TO NAME mapping, a data.frame of 2 column with term and name. Only used when gson is NULL.
#' @param minGSSize Minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize Maximal size of genes annotated for testing.
#' @param combineGroup A boolean. Whether the function combine the enrichment or not.
#' @param action A character string.A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param gson a GSON object, if not NULL, use it as annotation data.
#' @param keyType For KEGG analysis, keyType include ko, ec, cpd, entrezid. For Go analysis, keyType include entrezid and symbol.
#' @param KEGG_Type A character string. KEGG_Type include KEGG and MKEGG in KEGG analysis. KEGG means KEGG pathway. MKEGG means KEGG module.
#' @param species A character string. Species includ all, hsa, mmu,...in in KEGG analysis. Supported organism listed in 'https://www.genome.jp/kegg/catalog/
#' @param OrgDb OrgDb in Go analysis.
#' @param ont For Go analysis, ont include "BP", "MF","CC", and "ALL". For DOSE analysis, ont only support "DO".
#' @param organism For Reactome analysis, organism include "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly". For DOSE analysis, organism include "hsa" and "mmu".
#' @param ... Further parameters passed to \code{\link[clusterProfiler]{compareCluster}},\code{\link[clusterProfiler]{enricher}}, \code{\link[clusterProfiler]{enrichGO}}, \code{\link[ReactomePA]{enrichPathway}}, \code{\link[DOSE]{enrichDO}}.
#' @return EMPT object
#' @export
#' @section Detaild about method:
#' The EMP_enrich_analysis moudle performed based on cluserProfiler, more detailed information are available on:
#'
#' http://yulab-smu.top/biomedical-knowledge-mining-book/index.html
#' @examples
#' \dontrun{
#' data(MAE)
#' # Make the enrichment after EMP_diff_analysis
#' MAE |>
#'   EMP_assay_extract(experiment = 'geno_ec') |>
#'   EMP_diff_analysis(method='DESeq2',.formula = ~Group) |>
#'   EMP_enrich_analysis(keyType ='ec',KEGG_Type = 'KEGG',pvalue<0.05,pvalueCutoff=1,species = 'all') 
#' 
#' # Make the enrichment Visualization
#' MAE |>
#'   EMP_assay_extract(experiment = 'geno_ec') |>
#'   EMP_diff_analysis(method='DESeq2',.formula = ~Group) |>
#'   EMP_enrich_analysis(keyType ='ec',KEGG_Type = 'KEGG',pvalue<0.05,pvalueCutoff=1,species = 'all') |>
#'   EMP_enrich_dotplot()
#' 
#' MAE |>
#'   EMP_assay_extract(experiment = 'geno_ec') |>
#'   EMP_diff_analysis(method='DESeq2',.formula = ~Group) |>
#'   EMP_enrich_analysis(keyType ='ec',KEGG_Type = 'KEGG',pvalue<0.05,pvalueCutoff=1,species = 'all') |>
#'   EMP_enrich_netplot()
#' 
#' # Transcriptomic data
#' ## KEGG analysis
#' MAE |>
#'   EMP_assay_extract(experiment = 'host_gene') |>
#'   EMP_feature_convert(from = 'symbol',to='entrezid',species='Human') |>
#'   EMP_diff_analysis(method = 'DESeq2',.formula = ~Group,p.adjust = 'fdr') |> 
#'   EMP_enrich_analysis(pvalue<0.05,method = 'kegg',species = 'hsa',keyType='entrezid',
#'                       pvalueCutoff=0.05) 
#' 
#' ## GO analysis
#' library(org.Hs.eg.db)
#' MAE |>
#'   EMP_assay_extract(experiment = 'host_gene') |>
#'   EMP_feature_convert(from = 'symbol',to='entrezid',species='Human') |>
#'   EMP_diff_analysis(method = 'DESeq2',.formula = ~Group,p.adjust = 'fdr') |> 
#'   EMP_enrich_analysis(pvalue<0.05,method = 'go',OrgDb=org.Hs.eg.db,ont='MF',readable=TRUE,
#'                       pvalueCutoff=0.05) |>
#'   EMP_enrich_dotplot(show=6)
#' 
#' ## DOSE analysis
#' MAE |>
#'   EMP_assay_extract(experiment = 'host_gene') |>
#'   EMP_feature_convert(from = 'symbol',to='entrezid',species='Human') |>
#'   EMP_diff_analysis(method = 'DESeq2',.formula = ~Group,p.adjust = 'fdr') |> 
#'   EMP_enrich_analysis(pvalue<0.05,method = 'do',ont="DO",organism= 'hsa',readable=TRUE) |>
#'   EMP_enrich_dotplot(show=5)
#' 
#' ## Reactome analysis
#' MAE |>
#'   EMP_assay_extract(experiment = 'host_gene') |>
#'   EMP_feature_convert(from = 'symbol',to='entrezid',species='Human') |>
#'   EMP_diff_analysis(method = 'DESeq2',.formula = ~Group,p.adjust = 'fdr') |> 
#'   EMP_enrich_analysis(pvalue<0.05,method = 'Reactome',organism= 'human',readable=TRUE) |>
#'   EMP_enrich_dotplot()
#' }
EMP_enrich_analysis <- function(obj,condition,minGSSize=1,maxGSSize=500, action='add',combineGroup=FALSE,gson=NULL,
                               method = "kegg", 
                               TERM2GENE = NULL,
                               TERM2NAME = NA,
                               KEGG_Type='KEGG',species = "all", # kegg.params
                               OrgDb = NULL, keyType = "entrezid", # go.params
                               ont = if (method == "go") "MF" else "DO", # go.params and do.params
                               organism = if (method == "reactome") "human" else "Homo sapiens", # wikipathway.params and reactome.params
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

  rlang::check_installed(c('BiocManager'), reason = 'for EMP_enrich_analysis().', action = install.packages) 
  rlang::check_installed(c('clusterProfiler'), reason = 'for EMP_enrich_analysis().', action = BiocManager::install)  

  call <- match.call()

  EMPT <- .EMP_enrich_analysis(obj,{{condition}},minGSSize=minGSSize,maxGSSize=maxGSSize, combineGroup=combineGroup,gson=gson, 
                            method = method,
                            TERM2GENE = TERM2GENE,
                            TERM2NAME = TERM2NAME,
                            kegg.params = kegg.params,
                            go.params = go.params,
                            reactome.params = reactome.params,
                            wikipathway.params = wikipathway.params,
                            do.params = do.params,
                            ...)
  .get.history.EMPT(EMPT) <- call
  class(EMPT) <- 'EMP_enrich_analysis'
  if (action == 'add') {
    return(EMPT)
  }else if(action == 'get') {
    return(.get.result.EMPT(EMPT))
  }else{
    warning('Parameter action should be one of add or get!')
  }  
}




# #' Title
# #'
# #' @param data wait_for_add
# #' @param mapID wait_for_add
# #' @param metabolite_col wait_for_add
# #' @param gene_col wait_for_add
# #' @param gradient_by wait_for_add
# #' @param gradient_col wait_for_add
# #' @importFrom dplyr pull
# #' @importFrom ggraph ggraph
# #' @importFrom ggkegg geom_node_rect
# #' @importFrom ggplot2 scale_fill_gradient
# #' @importFrom ggkegg overlay_raw_map
# #' @importFrom ggplot2 theme_void
# #' @importFrom stringr str_replace_all
# #' @importFrom ggraph geom_node_point
# #'
# #' @return xx object
# #' @export
# #'
# #' @examples
# #' # add example
# EMP_enrich_kegg <- function(data,mapID,metabolite_col = 'blue',gene_col = 'red',gradient_by = NULL,gradient_col = c('pink','steelblue')) {

#   ID <- geneID <- `.` <- name <- x <- y <- mod <- type <- graphics_name <- NULL
#   diff_node <- data$enrich.data@compareClusterResult %>%
#     dplyr::filter(ID == mapID) %>%
#     dplyr::pull(geneID) %>%  strsplit(.,'/') %>% unlist

#   mapID %<>% gsub('map','ko',.)
#   diff_df <- data$diff_df %>% dplyr::filter(ID %in% diff_node) %>% dplyr::rename(graphics_name = 'ID')
#   # where is pathway?
#   g <- pathway(mapID,group_rect_nudge=0)

#   g %<>%
#       dplyr::mutate(graphics_name = stringr::str_replace_all(name, c("ko:" = "", "cpd:" = "", "path:" = ""))) %>%
#       # where is highlight_set_nodes?
#       dplyr::mutate(mod=highlight_set_nodes(diff_node, name="graphics_name",how =F))

#   if (is.null(gradient_by)) {
#       gg <- ggraph(g, layout="manual", x=x, y=y)+
#         geom_node_rect(fill=gene_col,aes(filter=mod))+
#         geom_node_point(aes(filter=type=="compound"), color=metabolite_col, size=2)+
#         overlay_raw_map(mapID)+
#         theme_void()
#   } else {
#       g %<>%
#         dplyr::mutate(!!gradient_by := sapply(graphics_name, function(x) diff_df[[gradient_by]][str_detect(x, diff_df$graphics_name)])) %>%
#         dplyr::mutate(!!gradient_by := as.character(!!dplyr::sym(gradient_by)),
#                !!gradient_by := as.numeric(!!dplyr::sym(gradient_by))) %>% suppressWarnings()
#       gg <- ggraph(g, layout="manual", x=x, y=y)+
#         geom_node_rect(aes(filter=mod,fill = !!dplyr::sym(gradient_by)))+
#         scale_fill_gradient(low=gradient_col[1],high=gradient_col[2], name=gradient_by) +
#         # geom_node_point(aes(filter=type=="compound"), color="blue", size=2)+
#         overlay_raw_map(mapID)+
#         #ggfx::with_outer_glow(geom_node_rect(aes(filter=mod), color="white",size=2),olour="yellow",expand=5)+
#         theme_void()
#     }
#     return(gg)
# }