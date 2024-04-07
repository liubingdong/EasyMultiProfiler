#' mitch method
#'
#' @param multiGene a data.frame of multi-omics gene difference analysis results (-log10 pvalue * sign(logFC)).
#' Each row is a gene, and each column represents an omics dataset. 
#' @param TERM2GENE user input annotation of TERM TO GENE mapping,
#' a data.frame of 2 column with term and gene
#' @param TERM2NAME user input of TERM TO NAME mapping, 
#' a data.frame of 2 column with term and name
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param pvalueCutoff Cutoff value of pvalue.
#' @param qvalueCutoff Cutoff of qvalue.
#' @param ... Other parameters.
#' @noRd
mitch_method <- function(multiGene, TERM2GENE, TERM2NAME = NULL, minGSSize, 
                         pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.25, ...) {
    PATHID2EXTID <- split(as.character(TERM2GENE[,2]), as.character(TERM2GENE[,1]))
    results <- mitch::mitch_calc(x = multiGene, genesets = PATHID2EXTID,
        minsetsize = minGSSize, ...)
    enrichment_result <- results$enrichment_result
    if (is.null(TERM2NAME)) {
        enrichment_result$Description <- enrichment_result$set     
    } else {
        rownames(TERM2NAME) <- TERM2NAME[, 1]
        enrichment_result$Description <- TERM2NAME[enrichment_result$set, 2]
    }

    # 先把multiGene保持不变
    geneSets <- split(TERM2GENE[, 2], TERM2GENE[, 1])
    gmts <- split(TERM2GENE, TERM2GENE[, 1])
    sets <- names(geneSets)
    Count <- lapply(geneSets, function(x) intersect(x, rownames(multiGene)) |> length()) |>
        unlist()
    names(Count) <- names(geneSets)
    enrichment_result$Count <- Count[enrichment_result$set]
    # Count <- vapply(em$overlap, length, FUN.VALUE = 1)
    nInputGene <- length(intersect(rownames(multiGene), TERM2GENE$gene))
    GeneRatio <- apply(data.frame(a=Count, b=nInputGene), 1, function(x)
                       paste(x[1], "/", x[2], sep="", collapse="")
                       )
    names(GeneRatio) <- names(Count)
    enrichment_result$GeneRatio <- GeneRatio[enrichment_result$set]

    nTermGene <- vapply(gmts[enrichment_result$set], nrow, FUN.VALUE = 1)
    N <- length(unique(TERM2GENE[, 2]))
    BgRatio <- apply(data.frame(a=nTermGene, b=N), 1, function(x)
                     paste(x[1], "/", x[2], sep="", collapse="")
                     )
    enrichment_result$BgRatio <- BgRatio
    p.adj <- enrichment_result$p.adjustMANOVA
    qobj <- tryCatch(qvalue(p=enrichment_result$p.adjustMANOVA, lambda=0.05, pi0.method="bootstrap"),
        error=function(e) NULL)
    qvalue <- p.adj
    if (!is.null(qobj)) {
        qvalue <- qobj$qvalue
    }
    overlap <- lapply(results$input_genesets, function(x) {
        intersect(x, rownames(multiGene))
    })
    geneID <- vapply(overlap, function(i) paste(i, collapse="/"), FUN.VALUE = "1")
    names(geneID) <- names(results$input_genesets)
    enrichment_result$geneID <- geneID[enrichment_result$set]
    



    result <- data.frame(ID          = enrichment_result$set,
                         Description = enrichment_result$Description,
                         GeneRatio   = enrichment_result$GeneRatio,
                         BgRatio     = enrichment_result$BgRatio,
                         pvalue      = enrichment_result$p.adjustMANOVA,
                         p.adjust    = p.adj,
                         qvalue      = qvalue,
                         geneID      = enrichment_result$geneID,
                         Count       = enrichment_result$Count)
    result <- result[result$pvalue < pvalueCutoff, ]
    result <- result[result$qvalue < qvalueCutoff, ]
    result <- result[order(result$pvalue), ]
    background <- unique(TERM2GENE[, 2])
    x <- new("enrichResult",
             result         = result,
             pvalueCutoff   = pvalueCutoff,
             pAdjustMethod  = "none",
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
