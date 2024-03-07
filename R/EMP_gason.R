#' @importFrom utils getFromNamespace
kegg_rest <- getFromNamespace("kegg_rest", "clusterProfiler")

#########
# 原版的仅支持pathway，不支持module，也不支持单物种。
gson_cpd2 <- function(KEGG_Type = "KEGG", species = "all") {
    if (KEGG_Type == "KEGG") {
        target <- "pathway"
    } else {
        target <- "module"
    }
    kegg_rest <- getFromNamespace("kegg_rest", "clusterProfiler")
    url_k1 <- paste0("https://rest.kegg.jp/link/cpd/", target, collapse = "")
    k1 <- kegg_rest(url_k1)
    k1[, 1]  <- gsub("[^:]+:", "", k1[, 1])
    k1[, 2]  <-  gsub("[^:]+:", "",  k1[, 2])
    if (KEGG_Type == "KEGG") {
        k1 <- k1[grep("map", k1[, 1]),]
    }
    url_k2 <- paste0("https://rest.kegg.jp/list/", target, collapse = "")
    k2  <- kegg_rest(url_k2)
    k2[, 1]  <- gsub("path:", "", k2[, 1])
    if (species != "all") {
        url_pathway2gene <- paste0("https://rest.kegg.jp/link/", species, "/", target, collapse="")
        pathway2gene <- clusterProfiler:::kegg_rest(url_pathway2gene)
        if (KEGG_Type == "KEGG") {
            pathways <- gsub(paste0("path:", species), "map", pathway2gene[, 1])
        } else {
            pathways <- gsub(paste0("md:", species, "_"), "", pathway2gene[, 1])
        }
        
        k1 <- k1[k1[, 1] %in% pathways, ]
    }
    
    gsid2gene <- setNames(k1, c("gsid", "gene"))
    gsid2name <- setNames(k2, c("gsid", "name"))
    y <- yulab.utils::yread("https://rest.kegg.jp/info/cpd")
    version <- sub("\\w+\\s+", "", y[grep('Release', y)])
    gson::gson(
        gsid2gene = gsid2gene,
        gsid2name = gsid2name,
        species = species,
        gsname = "KEGG",
        version = version,
        keytype = "kegg_compound",
        accessed_date = as.character(Sys.Date())
    )
}


# clusterProfiler有个同名函数gson_KEGG，因此改为gson_KEGG2
gson_KEGG2 <- function(keyType = "ko", KEGG_Type = "KEGG", species = "all") {
    keyType <- match.arg(keyType, c("ko", "ec"))
    if (KEGG_Type == "KEGG") {
        target <- "pathway"
    } else {
        target <- "module"
    }
    kegg_rest <- getFromNamespace("kegg_rest", "clusterProfiler")
    if (species == "all") {
        url_k1 <- paste0("https://rest.kegg.jp/link/", keyType, "/", target, collapse = "")
        k1 <- kegg_rest(url_k1)
    } else {
        url_pathway2gene <- paste0("https://rest.kegg.jp/link/", species, "/", target, collapse="")
        pathway2gene <- clusterProfiler:::kegg_rest(url_pathway2gene)
        if (KEGG_Type == "KEGG") {
            pathway2gene[, 1] <- gsub(paste0("path:", species), "map", pathway2gene[, 1])
        } else {
            pathway2gene[, 1] <- gsub(paste0("md:", species, "_"), "", pathway2gene[, 1])
        }
        url_id2gene <- paste0("https://rest.kegg.jp/link/", species, "/", keyType, collapse="")
        id2gene <- clusterProfiler:::kegg_rest(url_id2gene)     
        k1 <- merge(id2gene, pathway2gene, by = "to")[, c(3, 2)]
        k1[, 2]  <-  gsub("[^:]+:", "",  k1[, 2])
    }
   
    k1[, 1]  <- gsub("[^:]+:", "", k1[, 1])
    k1[, 2]  <-  gsub("[^:]+:", "",  k1[, 2])
    if (KEGG_Type == "KEGG") {
        k1 <- k1[grep("map", k1[, 1]),]
    }
    url_k2 <- paste0("https://rest.kegg.jp/list/", target, collapse = "")
    k2  <- kegg_rest(url_k2)
    k2[, 1]  <- gsub("path:", "", k2[, 1])
    gsid2gene <- setNames(k1, c("gsid", "gene"))
    gsid2name <- setNames(k2, c("gsid", "name"))
    y <- yulab.utils::yread("https://rest.kegg.jp/info/ko")
    version <- sub("\\w+\\s+", "", y[grep('Release', y)])
    if (keyType == "ko") {
        keytype <- "kegg_orthology"
    } else {
        keytype <- "kegg_enzyme"
    }
    gson::gson(
        gsid2gene = gsid2gene,
        gsid2name = gsid2name,
        species = species,
        gsname = "KEGG",
        version = version,
        keytype = keytype,
        accessed_date = as.character(Sys.Date())
    )
}

# 主函数
build_gson <- function(keyType = "ko", KEGG_Type = "KEGG", species = "hsa") {
    keyType <-  tolower(keyType) |> 
        match.arg(choices = c("ko", "ec", "cpd", "entrezid"))
    KEGG_Type <- match.arg(KEGG_Type, c("KEGG", "MKEGG"))

    if (keyType == "ko") {
        return(gson_KEGG2(keyType = "ko", KEGG_Type = KEGG_Type, species = species))
    }

    if (keyType == "ec") {
        return(gson_KEGG2(keyType = "ec", KEGG_Type = KEGG_Type, species = species))
    }
    
    if (keyType == "cpd") {
        return(gson_cpd2(KEGG_Type = KEGG_Type, species = species))
    }
    
    if (keyType == "entrezid") {
        if (species == "all") {
            stop("When keyType is 'entrezid', a specific species must be specified, not 'all'.")
        } else {
            return(clusterProfiler::gson_KEGG(species = species, KEGG_Type = KEGG_Type, keyType = "kegg"))
        }
        
    }

  
}

