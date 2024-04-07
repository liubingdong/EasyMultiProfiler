#' Multi-omics enrichment analysis
#'
#' @importFrom stats p.adjust
#' @importFrom stats setNames
#' @importFrom utils head
#' @importFrom utils tail
#' @importFrom qvalue qvalue
#' @importFrom methods new
#' @param multiGene a data.frame of multi-omics gene difference analysis results (pvalue).
#' Each row is a gene, and each column represents an omics dataset. 
#' @param network network
#' @param method enrichment analysis method, one of "enricher"(the default)
#' "GSEA", "mitch", "ActivePathways", and "multiNetEnrich".
#' @param cutoff Pvalue threshold of differentially expressed genes.
#' @param pvalueCutoff Cutoff value of pvalue.
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", 
#' "bonferroni", "BH", "BY", "fdr", "none"
#' @param combineLevel one of "gene" and "enrichResult"
#' @param universe background genes
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize maximal size of each geneSet for analyzing
#' @param qvalueCutoff Cutoff of qvalue.
#' @param TERM2GENE user input annotation of TERM TO GENE mapping, 
#' a data.frame of 2 column with term and gene
#' @param TERM2NAME user input of TERM TO NAME mapping, 
#' a data.frame of 2 column with term and name
#' @param combineMethod The method of combining pvalues, one of 
#' "fisher", "edgington", "stouffer" and "Brown"(only used in ActivePathways method).
#' @param stoufferWeights weights of stouffer combine method.
#' @param ... Other parameters.
#' @export
multiEnrichment <- function(multiGene,
                            method = "enricher",
                            network = NULL,
                            cutoff = 0.05,
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            universe = NULL,
                            minGSSize=10,
                            maxGSSize=500,
                            qvalueCutoff = 0.2,
                            TERM2GENE,
                            TERM2NAME = NULL,
                            combineMethod = "fisher",
                            stoufferWeights = NULL,
                            combineLevel = "enrichResult",
                            ...) {

    method <- match.arg(method,
        c("enricher", "GSEA", "mitch", "ActivePathways", "multiNetEnrich"))
    if (method == "mitch") {
        em <- mitch_method(multiGene, TERM2GENE, TERM2NAME = TERM2NAME, 
                           minGSSize = minGSSize, pvalueCutoff = pvalueCutoff, 
                           pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff, ...)

    }

    if (method == "ActivePathways") {
        em <- ActivePathways_method(multiGene = multiGene,
                                    pvalueCutoff = pvalueCutoff,
                                    pAdjustMethod = pAdjustMethod,
                                    universe = universe,
                                    minGSSize = minGSSize,
                                    maxGSSize = maxGSSize,
                                    qvalueCutoff = qvalueCutoff,
                                    TERM2GENE = TERM2GENE,
                                    TERM2NAME = TERM2NAME,
                                    combineMethod = combineMethod,
                                    ...)
    }

    if (method == "enricher") {
        em <- multi_enricher(multiGene = multiGene,
                             cutoff = cutoff,
                             pvalueCutoff = pvalueCutoff,
                             pAdjustMethod = pAdjustMethod,
                             universe = universe,
                             minGSSize = minGSSize,
                             maxGSSize = maxGSSize,
                             qvalueCutoff = qvalueCutoff,
                             TERM2GENE = TERM2GENE,
                             TERM2NAME = TERM2NAME,
                             combineMethod = combineMethod,
                             stoufferWeights = stoufferWeights,
                             combineLevel = combineLevel,
                             ...)
    }

    if (method == "GSEA") {
        em <- multi_GSEA(multiGene = multiGene,
                         pvalueCutoff = pvalueCutoff,
                         pAdjustMethod = pAdjustMethod,
                         universe = universe,
                         minGSSize = minGSSize,
                         maxGSSize = maxGSSize,
                         qvalueCutoff = qvalueCutoff,
                         TERM2GENE = TERM2GENE,
                         TERM2NAME = TERM2NAME,
                         combineMethod = combineMethod,
                         stoufferWeights = stoufferWeights,
                         combineLevel = combineLevel,
                         ...)
    }
    
    if (method == "multiNetEnrich") {
        em <- multiNetEnrich(multiGene, network = network, TERM2GENE = TERM2GENE,
                             TERM2NAME = TERM2NAME, 
                             pvalueCutoff = pvalueCutoff, cutoff = cutoff,
                             pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff,
                             combineMethod = combineMethod,
                             stoufferWeights = stoufferWeights, 
                             combineLevel = combineLevel,
                             ...)
    }
    return(em)
}


