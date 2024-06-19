#' Title
#'
#' @param EMPT wait_for_add
#' @param condition Expressions that return a logical value according to the result of EMP_diff_analysis. eg. pvalue < 0.05
#' @param keyType A character string. Methods include ko, ec, cpd and entrezid.
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize maximal size of genes annotated for testing.
#' @param use_cache A boolean. Whether the function use the results in cache or re-compute.
#' @param ... Further parameters passed to clusterProfiler::compareCluster.
#' @importFrom clusterProfiler enricher
#' @importFrom clusterProfiler compareCluster
#' @noRd
.EMP_enrich_analysis <- function(EMPT,condition,minGSSize =1,maxGSSize =500,keyType=NULL,KEGG_Type='KEGG',species = "all",combineGroup=FALSE,...){
  deposit <- list()
  sign_group <- NULL
  condition <- dplyr::enquo(condition)
  df <- .get.result.EMPT(EMPT,info='EMP_diff_analysis') %>% suppressMessages() %>%
    dplyr::filter(!!condition) %>%
    tidyr::drop_na(sign_group) ## filter NA or the result will add NA group!
  
  if (is.null(keyType)) {
    stop("keyType should be specified as ko, ec or cpd!")
  }else if(!keyType %in% c('ko','ec','cpd','entrezid')){
    stop("keyType should be ko, ec, cpd or entrezid!")
  }

  if(!KEGG_Type %in% c('KEGG','MKEGG')){
    stop("keyType should be KEGG or MKEGG!")
  }

  gson_data <- build_gson(keyType = keyType, KEGG_Type = KEGG_Type, species = species)
  
  message('KEGG database version: ',gson_data@version)
  message('Species: ',gson_data@species)
  if (combineGroup == TRUE) {
    enrich.data <- clusterProfiler::enricher(gene=df$feature, gson=gson_data,minGSSize=minGSSize, maxGSSize=maxGSSize,...) 
  }else{
    enrich.data <- clusterProfiler::compareCluster(feature~sign_group, data=df, 
          fun=clusterProfiler::enricher, gson=gson_data, 
          minGSSize=minGSSize, maxGSSize=maxGSSize,...)     
  }

  EMPT@deposit[['enrich_data']] <- enrich.data
  .get.algorithm.EMPT(EMPT) <- 'enrich_analysis'
  .get.info.EMPT(EMPT) <- 'EMP_enrich_analysis' 
  return(EMPT) 
}

#' KEGG enrichment for EMPT object
#'
#' @param obj Object in EMPT or MultiAssayExperiment format.
#' @param condition Expressions that return a logical value according to the result of EMP_diff_analysis. eg. pvalue < 0.05
#' @param minGSSize Minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize Maximal size of genes annotated for testing.
#' @param keyType A character string. keyType include ko, ec, cpd, entrezid.
#' @param KEGG_Type A character string. KEGG_Type include KEGG and MKEGG.
#' @param species A character string. Species includ all, hsa, mmu,...Supported organism listed in 'https://www.genome.jp/kegg/catalog/org_list.html'
#' @param combineGroup A boolean. Whether the function combine the enrichment or not.
#' @param action A character string.A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param ... Further parameters passed to clusterProfiler::compareCluster or clusterProfiler::enricher.
#'
#' @return EMPT object
#' @export
#'
#' @examples
#' \dontrun{
#' data(MAE)
#' ## Make the enrichment after EMP_diff_analysis
#' MAE |>
#'   EMP_assay_extract(experiment = 'geno_ec') |>
#'   EMP_diff_analysis(method='DESeq2',.formula = ~Group) |>
#'   EMP_enrich_analysis(keyType ='ec',KEGG_Type = 'KEGG',pvalue<0.05,pvalueCutoff=1,species = 'all') 
#' 
#' ## Make the enrichment after EMP_diff_analysis and Visualization
#' MAE |>
#'   EMP_assay_extract(experiment = 'geno_ec') |>
#'   EMP_diff_analysis(method='DESeq2',.formula = ~Group) |>
#'   EMP_enrich_analysis(keyType ='ec',KEGG_Type = 'KEGG',pvalue<0.05,pvalueCutoff=1,species = 'all') |>
#'   EMP_dotplot()
#' 
#' MAE |>
#'   EMP_assay_extract(experiment = 'geno_ec') |>
#'   EMP_diff_analysis(method='DESeq2',.formula = ~Group) |>
#'   EMP_enrich_analysis(keyType ='ec',KEGG_Type = 'KEGG',pvalue<0.05,pvalueCutoff=1,species = 'all') |>
#'   EMP_netplot()
#' 
#' 
#' ## Transcriptomic data
#' MAE |>
#'   EMP_assay_extract(experiment = 'host_gene') |>
#'   EMP_feature_convert(from = 'symbol',to='entrezid',species='Human') |>
#'   EMP_diff_analysis(method = 'DESeq2',.formula = ~Group,p.adjust = 'fdr') |> 
#'   EMP_enrich_analysis(keyType ='entrezid',KEGG_Type = 'KEGG',pvalue<0.05,pvalueCutoff=0.05,species = 'hsa') |>
#'   EMP_dotplot()
#' }
EMP_enrich_analysis <- function(obj,condition,minGSSize=1,maxGSSize=500,keyType=NULL,KEGG_Type='KEGG',species = "all",action='add',combineGroup=FALSE,...){
  
  rlang::check_installed(c('BiocManager'), reason = 'for EMP_enrich_analysis().', action = install.packages) 
  rlang::check_installed(c('clusterProfiler'), reason = 'for EMP_enrich_analysis().', action = BiocManager::install)  

  call <- match.call()

  if(!is.null(.get.result.EMPT(obj,info='EMP_diff_analysis'))) {
    EMPT <- obj
  }else {
    stop('The input data be generated from EMP_diff_analysis!')
  }


  EMPT <- .EMP_enrich_analysis(EMPT,{{condition}},minGSSize=minGSSize,maxGSSize=maxGSSize,keyType=keyType,KEGG_Type=KEGG_Type,species=species,combineGroup=combineGroup,...)
  .get.history.EMPT(EMPT) <- call
  class(EMPT) <- 'EMP_enrich_analysis'
  if (action == 'add') {
    return(EMPT)
  }else if(action == 'get') {
    return(.get.result.EMPT(EMPT))
  }else{
    warning('action should be one of add or get!')
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