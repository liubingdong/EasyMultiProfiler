
#' Multi analysis for diff reslut or enrichment
#'
#' @param EMP EMP objetc
#' @param select A character string. The experiment name in the EMP object.
#' @param method A character string including feature, diff_feature_enrich and same_feature_enrich.
#' @param combineFun A character string including enricher, ActivePathways and mitch.
#' @param pvalueCutoff A character string. Adjusted pvalue cutoff on enrichment tests to report.
#' @param combineMethod A character string including fisher, edgington and stouffer. Only actived, when combineFun == "enricher".
#' @param combineGroup A boolean. Whether the function combine the enrichment or not.
#' @param combineLevel A character string. One of "gene" and "enrichResult". Only actived, when combineFun == "enricher".
#' @param minGSSize Minimal size of genes annotated by Ontology term for testing. (default=10)
#' @param maxGSSize Maximal size of each geneSet for analyzing. default=500)
#' @param keyType A character string. Methods include ko, ec, cpd and entrezid.
#' @param KEGG_Type A character string. KEGG_Type include KEGG and MKEGG.
#' @param species A character string. Species includ all, hsa, mmu,...Supported organism listed in 'https://www.genome.jp/kegg/catalog/org_list.html'
#' @param p.adjust A character string. One of holm, hochberg, hommel, bonferroni, BH, BY, fdr, none.(defaault: fdr)
#' @param action A character string.A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param ... ...
#' @rdname EMP_multi_analysis
#' @return EMPT object
#' @export
#'
#' @examples
#' # add example
EMP_multi_analysis <- function(EMP,select=NULL,method='feature',combineFun='enricher',
                                 pvalueCutoff=0.05,combineMethod='fisher',combineGroup=FALSE,
                                 combineLevel='gene',p.adjust='fdr',minGSSize=10,maxGSSize=500,
                                 keyType = "ko", KEGG_Type = "KEGG", species = "all",action='add',...) {
  experment_num <- deposit <- ExperimentList <- experment_num <- feature <- pvalue <- deposit_enrichment <- gson_data <- input_list_each <- NULL
  
  if (!inherits(EMP,"EMP")) {
    stop("Please input the EMP format!")
  }

  call <- match.call()
  if (method == 'feature') {
    input_list <- list()

    if (is.null(select)) {
        ExperimentList <- .get.ExperimentList.EMP(EMP) 
    }else if (!is.null(select)) {
        ExperimentList <- .get.ExperimentList.EMP(EMP)[select]
    }

    experment_num <- ExperimentList %>% length()
    for (i in seq_len(experment_num)) {
      input_list[[i]] <- ExperimentList[[i]] %>% 
        EMP_result(info='diff_analysis_result') %>% 
        dplyr::select(feature,pvalue) %>% as.data.frame() %>% ## as.data.frame is necessary or wrong in .EMP_multi_feature!
        suppressMessages()
    }     

    deposit_df <- .EMP_multi_feature(list_df=input_list,combineMethod=combineMethod,p.adjust=p.adjust)

    if (action=='get') {
        return(deposit_df)
    }else if(action=='add') {
        EMP@deposit[['multi_same_df']] <- deposit_df
        .get.history.EMP(EMP) <- call
        .get.info.EMP(EMP) <-'EMP_multi_same_df'
        .get.method.EMP(EMP) <- method
        class(EMP) <- 'EMP_multi_same_df'
        return(EMP)       
    }else{
        warning('action should be one of add or get!')
    }

  }else if (method == 'diff_feature_enrich') {
    input_list <- list()

    if (is.null(select)) {
      ExperimentList <- .get.ExperimentList.EMP(EMP) 
    }else if (!is.null(select)) {
      ExperimentList <- .get.ExperimentList.EMP(EMP)[select]
    }

    experment_num <- ExperimentList %>% length()
    for (i in seq_len(experment_num)) {
      input_list_each <- ExperimentList[[i]] %>% 
        EMP_result(info='enrich_data') %>% suppressMessages()

      if (is.null(input_list_each)) {
        experiment_name <- names(ExperimentList[[i]])
        stop("For ",experiment_name,", it needs run EMP_enrich_analysis first!")
      }
      input_list[[i]] <- input_list_each
    } 
    
    deposit_enrichment <- .diff_enrich(list_compareCluster=input_list,combineMethod=combineMethod,p.adjust=p.adjust)

    if (action=='get') {
        if (combineGroup==FALSE) {
            return(deposit_enrichment@compareClusterResult %>% tibble::as_tibble() %>% print(n=Inf))
        }else if (combineGroup==TRUE) {
            return(deposit_enrichment@result %>% tibble::as_tibble() %>% print(n=Inf))
        }
    }else if(action=='add') {
        EMP@deposit[['multi_diff_enrich']] <- deposit_enrichment
        .get.history.EMP(EMP) <- call
        .get.info.EMP(EMP) <-'EMP_multi_diff_enrich'
        .get.method.EMP(EMP) <- method
        class(EMP) <- 'EMP_multi_diff_enrich'
        return(EMP)       
    }else{
        warning('action should be one of add or get!')
    }

  }else if (method == 'same_feature_enrich') {
    gson_data <- build_gson(keyType = keyType,KEGG_Type = KEGG_Type,species = species)
    deposit_enrichment <- .same_enrich(obj=EMP,combineFun=combineFun,pvalueCutoff=pvalueCutoff,combineMethod=combineMethod,combineGroup=combineGroup,
                         gson=gson_data,combineLevel=combineLevel,p.adjust=p.adjust,minGSSize=minGSSize,maxGSSize=maxGSSize,...) %>% suppressWarnings()

    if (action=='get') {
        if (combineGroup==FALSE) {
            return(deposit_enrichment@compareClusterResult %>% tibble::as_tibble() %>% print(n=Inf))
        }else if (combineGroup==TRUE) {
            return(deposit_enrichment@result %>% tibble::as_tibble() %>% print(n=Inf))
        }
    }else if(action=='add') {
        EMP@deposit[['multi_same_enrich']] <- deposit_enrichment
        .get.history.EMP(EMP) <- call
        .get.info.EMP(EMP) <-'EMP_multi_same_enrich'
        .get.method.EMP(EMP) <- method
        class(EMP) <- 'EMP_multi_same_enrich'
    return(EMP)  
    
  }else{
    stop('Parameter type in EMP_multi_analysis only allow feature and enrich')
  }  

 }

}





#' @importFrom metap sumlog
#' @importFrom metap sump
#' @importFrom metap sumz
#' @noRd
.EMP_multi_feature <- function(list_df, combineMethod = "fisher", 
                           p.adjust = "fdr") {
    genes <- list_df |> lapply(function(x) x[, 1]) |> do.call(c, args = _) |> unique()   
    gene_df <- matrix(NA, nrow = length(genes), ncol = length(list_df))
    for (i in seq_len(ncol(gene_df))) {
        gene_df[, i] <- list_df[[i]][match(genes, list_df[[i]][, 1]), 2]
    }
    
    pvalue <- rep(0, nrow(gene_df))
    for (i in seq_len(nrow(gene_df))) {
        pp <- as.numeric(gene_df[i, ])
        delete <- which(is.na(pp))
        if (length(delete) > 0) pp <- pp[-delete]
        if (length(pp) == 1) {
            pvalue[i] <- pp
        } else {
            if (combineMethod == "fisher") pvalue[i] <- sumlog(pp)$p %>% suppressWarnings() # for all NA of pvalue
            if (combineMethod == "edgington") pvalue[i] <- sump(pp)$p %>% suppressWarnings()
            if (combineMethod == "stouffer") pvalue[i] <- sumz(pp) %>% suppressWarnings()
        }

    }
     p.adj <- p.adjust(pvalue, method=p.adjust)
     deposit <- tibble::tibble(feature = genes, pvalue = pvalue, {{p.adjust}} := p.adj)

}




.diff_enrich <- function(list_compareCluster,combineMethod='fisher',p.adjust='fdr') {
    result <- list_compareCluster[[1]]  
    compareClusterResult_list <- lapply(list_compareCluster, function(x) x@compareClusterResult)
    for (i in seq_len(length(compareClusterResult_list))) {
        compareClusterResult_list[[i]]$source <- paste0("source", i)
    }
    compareClusterResult_df <- do.call(rbind, compareClusterResult_list)
    
    # 按照cluster分开，然后将不同source的合并
    compareClusterResult_df_list <- split(compareClusterResult_df, compareClusterResult_df$Cluster)
    compareClusterResult_df_list <- lapply(compareClusterResult_df_list, combine_compareClusterResult,combineMethod=combineMethod,p.adjust=p.adjust)
    compare_result <- do.call(rbind, compareClusterResult_df_list)
    compare_result <- compare_result[!duplicated(compare_result$ID), 1:(ncol(compare_result)-1)]
    # 将同样的 ID的通路合并
    geneClusters_list <- lapply(list_compareCluster, function(x) x@geneClusters)
    geneClusters <- combine_geneClusters(geneClusters_list)
    names(geneClusters) <- names(list_compareCluster[[1]]@geneClusters)
    group_length <- lapply(geneClusters, length) |> unlist()

    compare_result$GeneRatio <- paste(compare_result$Count, group_length[as.character(compare_result$Cluster)], sep = "/")
    
    result@compareClusterResult <- compare_result
    
    result@geneClusters <- geneClusters
    result
}




#' @importFrom metap sumlog
#' @importFrom metap sump
#' @importFrom metap sumz
combine_pvalue_for_diff_enrich <- function(x, combineMethod = "fisher") {
    if (length(x) == 1) { 
       return(x)
    } else {
        if (combineMethod == "fisher") x <- sumlog(x)$p %>% suppressWarnings() # for all NA of pvalue
        if (combineMethod == "edgington") x <- sump(x)$p %>% suppressWarnings() # for all NA of pvalue
        if (combineMethod == "stouffer") x <- sumz(x)  %>% suppressWarnings() # for all NA of pvalue
        return(x)
    }
}

combine_geneID <- function(x) {
    strsplit(x, "/") |> do.call(c, args = _) |> unique() |> paste(collapse = "/")
}

combine_Count <- function(x) {
    strsplit(x, "/") |> do.call(c, args = _) |> unique() |> length()
}

combine_compareClusterResult <- function(x, combineMethod = "fisher", 
                           p.adjust = "fdr"){
   pvalue <- split(x$pvalue, x$ID) |> lapply(combine_pvalue_for_diff_enrich,combineMethod=combineMethod) |> unlist()
   geneID <- split(x$geneID, x$ID) |> lapply(combine_geneID) |> unlist()
   Count <- split(x$geneID, x$ID) |> lapply(combine_Count) |> unlist()
   # BgRatio 暂时不计算
   x$pvalue <- pvalue[x$ID]
   x$geneID <- geneID[x$ID]
   x$Count <- Count[x$ID]
   x$p.adjust <- p.adjust(x$pvalue, method=p.adjust) 
   return(x)   
}
 
combine_geneClusters <- function(geneClusters_list) {
    geneClusters_result <- vector("list", length(geneClusters_list[[1]]))
    for (i in seq_len(length(geneClusters_result))) {
        aa <- lapply(geneClusters_list, function(x) x[[i]])
        geneClusters_result[[i]] <- do.call(c, aa) |> unique()
    }
    return(geneClusters_result)
}

.same_enrich <- function(obj,combineFun='enricher',pvalueCutoff=0.05,combineMethod='fisher',combineGroup=FALSE,
                         gson,combineLevel='gene',p.adjust='fdr',minGSSize=10,maxGSSize=500,...) {
  
  
  experiment_name <- experment_list <- feature <- pvalue <- sign_group <- NULL
  ## avoid typo
  if (combineFun == 'activepathways') {
    combineFun <- 'ActivePathways'
  }
  
  pvalue_list <- list()
  experment_list <- .get.ExperimentList.EMP(obj)
  experiment_name <- names(experment_list)
  if (combineGroup == TRUE) {
    for (i in experiment_name) {
      pvalue_list[[i]] <- EMP_result(experment_list[[i]],info='diff_analysis_result') %>%
        dplyr::select(feature,pvalue) %>%
        tidyr::drop_na(pvalue) %>% 
        dplyr::rename({{i}} := pvalue)
    }
    merge_p_df <- purrr::reduce(pvalue_list, function(x, y) dplyr::full_join(x, y, by = "feature")) %>%
      tibble::column_to_rownames('feature')
    merge_enrich <- multiEnrichment(multiGene = merge_p_df, method = combineFun, pvalueCutoff = pvalueCutoff,combineMethod=combineMethod,
                                    TERM2GENE = gson@gsid2gene,TERM2NAME = gson@gsid2name,combineLevel=combineLevel,pAdjustMethod=p.adjust,
                                    minGSSize=minGSSize,maxGSSize=maxGSSize,...)
    return(merge_enrich)
  }else if (combineGroup == FALSE) {
    for (i in experiment_name) {
      pvalue_list[[i]] <- EMP_result(experment_list[[i]],info='diff_analysis_result') %>%
        dplyr::select(feature,pvalue,sign_group) %>%
        tidyr::drop_na(pvalue) %>% 
        dplyr::rename({{i}} := pvalue)
    }
    merge_p_df <- purrr::reduce(pvalue_list, function(x, y) dplyr::full_join(x, y, by = c("feature","sign_group"))) 
    
    group_name <- merge_p_df  %>% dplyr::pull(sign_group) %>%  unique()
    
    merge_p_df_list <- list()
    
    for (i in group_name) {
      merge_p_df_list[[i]] <- merge_p_df %>% 
        dplyr::filter(sign_group=={{i}}) %>% 
        dplyr::select(-sign_group) %>%
        tibble::column_to_rownames('feature')
    }
    
    enrichment_list <- list()
    for (i in group_name) {
      enrichment_list[[i]] <- multiEnrichment(multiGene = merge_p_df_list[[i]], method = combineFun, pvalueCutoff = pvalueCutoff,
                                              TERM2GENE = gson@gsid2gene,TERM2NAME = gson@gsid2name,
                                              combineLevel=combineLevel,combineMethod=combineMethod,pAdjustMethod=p.adjust,...) 
    }
    
    merge_enrich <-  clusterProfiler::merge_result(enrichment_list)
    return(merge_enrich)
  }else{
    stop("combineGroup must be TRUE or FALSE!")
  }  
}



