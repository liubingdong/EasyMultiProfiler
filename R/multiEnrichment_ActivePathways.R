#' ActivePathways method
#'
#' @param multiGene a data.frame of multi-omics gene difference analysis results (pvalue).
#' Each row is a gene, and each column represents an omics dataset. 
#' @param pvalueCutoff Cutoff value of pvalue.
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", 
#' "bonferroni", "BH", "BY", "fdr", "none"
#' @param universe background genes
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize maximal size of each geneSet for analyzing
#' @param qvalueCutoff Cutoff of qvalue.
#' @param TERM2GENE user input annotation of TERM TO GENE mapping, 
#' a data.frame of 2 column with term and gene
#' @param TERM2NAME user input of TERM TO NAME mapping, 
#' a data.frame of 2 column with term and name
#' @param combineMethod The method of combining pvalues, one of 
#' "fisher", and "Brown".
#' @param ... Other parameters.
#' @noRd
ActivePathways_method <- function(multiGene,
                                  pvalueCutoff = 0.05,
                                  pAdjustMethod = "BH",
                                  universe = NULL,
                                  minGSSize=10,
                                  maxGSSize=500,
                                  qvalueCutoff = 0.2,
                                  TERM2GENE,
                                  TERM2NAME = NULL,
                                  combineMethod = "fisher",
                                  ...) {
    rlang::check_installed(c('BiocManager'), reason = 'for ActivePathways_method().', action = install.packages)  
    rlang::check_installed(c('ActivePathways'), reason = 'for ActivePathways_method().', action = BiocManager::install)
    
    TERM2GENE$term <- as.character(TERM2GENE[, 1])
    TERM2GENE$gene <- as.character(TERM2GENE[, 2])
    geneSets <- split(TERM2GENE$gene, TERM2GENE$term)
    idx <- get_geneSet_index(geneSets, minGSSize, maxGSSize)
    if (sum(idx) == 0) {
        msg <- paste("No gene set have size >", minGSSize, "...")
        message(msg)
        message("--> return NULL...")
        return (NULL)
    }
    gmts <- split(TERM2GENE, TERM2GENE$term)
    # geneSet_size <- sapply(gmts, nrow)
    gmts <- gmts[idx]    
    gmts2 <- lapply(gmts, function(x) { list(id=x$term[1], name=x$term[1], genes=x$gene)})
    if (!is.null(TERM2NAME)) {
      rownames(TERM2NAME) <- TERM2NAME[, 1]
      termName <- TERM2NAME[names(gmts), 2]
      for (k in seq(length(gmts2))) {
          gmts2[[k]]$name <- termName[k]
      }
    }   
    class(gmts2) <- 'GMT'
    scores <- as.matrix(multiGene)
    # scores can not contain missing values
    scores[is.na(scores)] <- 1
    background <- TERM2GENE$gene |> as.character() |> unique()
    if(!is.null(universe)) {
        if (is.character(universe)) {
            universe <- intersect(background, universe)
            background <- universe
        } else {
            message("`universe` is not in character and will be ignored...")
        }
    }

    if (!combineMethod %in% c("Brown", "fisher")) {
        message(paste("combineMethod can only be 'fisher' or 'Brown'",
            "in ActivePathways method.",
            "'combineMethod' parameter will be set to 'Brown'."))
    }

    if (combineMethod == "fisher") {
        merge.method <- "Fisher"
    } else {
        merge.method <- "Brown"
    }
    em <- ActivePathways::ActivePathways(scores = scores, gmt = gmts2,
        significant = 1, background = background, correction_method = "none", 
        merge_method = merge.method, ...)

    Count <- vapply(em$overlap, length, FUN.VALUE = 1)
    nInputGene <- length(intersect(rownames(multiGene), TERM2GENE$gene))
    GeneRatio <- apply(data.frame(a=Count, b=nInputGene), 1, function(x)
                       paste(x[1], "/", x[2], sep="", collapse="")
                       )
    nTermGene <- vapply(gmts[em$term_id], nrow, FUN.VALUE = 1)
    N <- length(unique(TERM2GENE$gene))
    BgRatio <- apply(data.frame(a=nTermGene, b=N), 1, function(x)
                     paste(x[1], "/", x[2], sep="", collapse="")
                     )
    p.adj <- p.adjust(em$adjusted_p_val, method=pAdjustMethod)
    qobj <- tryCatch(qvalue(p=em$adjusted_p_val, lambda=0.05, pi0.method="bootstrap"),
        error=function(e) NULL)
    if (is.null(qobj)) {
    qvalues <- p.adj
    } else {
        qvalues <- qobj$qvalue
    }
    geneID <- vapply(em$overlap, function(i) paste(i, collapse="/"), FUN.VALUE = "1")
    result <- data.frame(ID          = em$term_id,
                         Description = em$term_name,
                         GeneRatio   = GeneRatio,
                         BgRatio     = BgRatio,
                         pvalue      = em$adjusted_p_val,
                         p.adjust    = p.adj,
                         qvalue      = qvalues,
                         geneID      = geneID,
                         Count       = Count)
    result <- result[result$pvalue < pvalueCutoff, ]
    result <- result[result$qvalue < qvalueCutoff, ]
    result <- result[order(result$pvalue), ]
    x <- new("enrichResult",
             result         = result,
             pvalueCutoff   = pvalueCutoff,
             pAdjustMethod  = pAdjustMethod,
             qvalueCutoff   = qvalueCutoff,
             gene           = rownames(multiGene),
             universe       = background,
             geneSets       = geneSets,
             organism       = "UNKNOWN",
             keytype        = "UNKNOWN",
             ontology       = "UNKNOWN",
             readable       = FALSE
            )
    return(x)
}
