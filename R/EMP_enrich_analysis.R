#' Title
#'
#' @param EMPT wait_for_add
#' @param condition wait_for_add
#' @param type wait_for_add
#' @param minGSSize wait_for_add
#' @param maxGSSize wait_for_add
#' @param use_cache wait_for_add
#' @param ... wait_for_add
#' @importFrom clusterProfiler compareCluster
#' @importFrom clusterProfiler enricher
#'
#' @noRd
.EMP_enrich_analysis <- function(EMPT,condition,minGSSize =1,maxGSSize =500,keyType=NULL,KEGG_Type='KEGG',species = "all",...){
  deposit <- list()
  condition <- dplyr::enquo(condition)
  df <- .get.result.EMPT(EMPT,info='EMP_diff_analysis') %>% suppressMessages() %>%
    dplyr::filter(!!condition)
  
  if (is.null(keyType)) {
    stop("keyType should be specified as ko, ec or cpd!")
  }else if(!keyType %in% c('ko','ec','cpd')){
    stop("keyType should be ko, ec or cpd!2")
  }

  if(!KEGG_Type %in% c('KEGG','MKEGG')){
    stop("keyType should be KEGG or MKEGG!")
  }

  gason_data <- build_gson(keyType = keyType, KEGG_Type = KEGG_Type, species = species)
  
  message('KEGG database version: ',gason_data@version)
  message('Species: ',gason_data@species)

  enrich.data <- clusterProfiler::compareCluster(feature~sign_group, data=df, 
        fun=clusterProfiler::enricher, gson=gason_data, 
        minGSSize=minGSSize, maxGSSize=maxGSSize,...) 

  EMPT@deposit[['enrich_data']] <- enrich.data
  .get.algorithm.EMPT(EMPT) <- 'enrich_analysis'
  .get.info.EMPT(EMPT) <- 'EMP_enrich_analysis' 
  return(EMPT) 
}

#' Title
#'
#' @param x wait_for_add
#' @param condition filter condition
#' @param minGSSize wait_for_add
#' @param maxGSSize wait_for_add
#' @param keyType wait_for_add
#' @param KEGG_Type wait_for_add
#' @param species wait_for_add
#' @param action wait_for_add
#' @param ... wait_for_add
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_enrich_analysis <- function(x,condition,minGSSize=1,maxGSSize=500,keyType=NULL,KEGG_Type='KEGG',species = "all",action='add',...){
  call <- match.call()

  if(!is.null(.get.result.EMPT(x,info='EMP_diff_analysis'))) {
    EMPT = x
  }else {
    stop('The input data be generated from EMP_diff_analysis!')
  }


  EMPT <- .EMP_enrich_analysis(EMPT,{{condition}},minGSSize=minGSSize,maxGSSize=maxGSSize,keyType=keyType,KEGG_Type=KEGG_Type,species=species,...)
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
