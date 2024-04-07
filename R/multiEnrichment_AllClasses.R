##' Class "multiEnrichResult"
##' This class represents the result of multi-omics ORA enrichment analysis.
##'
##'
##' @name multiEnrichResult-class
##' @aliases multiEnrichResult-class
##'   show,multiEnrichResult-method summary,multiEnrichResult-method
##'
##' @docType class
##' @slot result multi-omics enrichment analysis.
##' @slot geneID enriched genes.
##' @slot Count number of enriched genes.
##' @slot GeneRatio Ratio of the number of genes enriched 
##' to the number of genes input.
##' @slot pvalueCutoff pvalueCutoff.
##' @slot pAdjustMethod pvalue adjust method.
##' @slot qvalueCutoff qvalueCutoff.
##' @slot organism organism.
##' @slot ontology biological ontology.
##' @slot gene input gene IDs.
##' @slot keytype Gene ID type.
##' @slot universe background gene.
##' @slot gene2Symbol mapping gene to Symbol.
##' @slot geneSets gene sets.
##' @slot readable logical flag of gene ID in symbol or not.
##' @slot termsim Similarity between term.
##' @slot method method of calculating the similarity between nodes.
##' @exportClass multiEnrichResult
##' @author Erqiang Hu 
##' @seealso \code{\link{multiEnrichment}}
##' @keywords classes
setClass("multiEnrichResult",
         representation=representation(
             result         = "data.frame",
             geneID         = "list",
             Count          = "list",
             GeneRatio      = "list",
             pvalueCutoff   = "numeric",
             pAdjustMethod  = "character",
             qvalueCutoff   = "numeric",
             organism       = "character",
             ontology       = "character",
             gene           = "list",
             keytype        = "character",
             universe       = "character",
             gene2Symbol    = "character",
             geneSets       = "list",
             readable       = "logical",
             termsim        = "matrix",
             method         = "character"
             ),
         prototype=prototype(readable = FALSE)
)



##' Class "multiGseaResult"
##' This class represents the result of multi-omics GSEA analysis
##'
##'
##' @name multiGseaResult-class
##' @aliases gseahResult-class
##'   show,multiGseaResult-method summary,multiGseaResult-method
##'
##' @docType class
##' @slot result multi-omics GSEA anaysis.
##' @slot enrichmentScore enrichmentScore.
##' @slot NES Enrichment score after normalization.
##' @slot rank rank
##' @slot leading_edge leading_edge.
##' @slot core_enrichment core enrichment genes.
##' @slot organism organism.
##' @slot setType setType.
##' @slot geneSets geneSets.
##' @slot geneList order rank geneList.
##' @slot keytype ID type of gene.
##' @slot permScores permutation scores.
##' @slot params parameters.
##' @slot gene2Symbol gene ID to Symbol.
##' @slot readable whether convert gene ID to symbol.
##' @exportClass multiGseaResult
##' @author Erqiang Hu
##' @keywords classes
setClass("multiGseaResult",
         representation   = representation(
             result          = "data.frame",
             enrichmentScore = "list",
             NES             = "list",
             rank            = "list",
             leading_edge    = "list",
             core_enrichment = "list",
             geneList        = "list",
             setSize         = "list",
             organism        = "character",
             setType         = "character",
             geneSets        = "list",
             keytype         = "character",
             permScores      = "list",
             params          = "list",
             gene2Symbol     = "character",
             readable        = "logical",
             termsim         = "matrix",
             method          = "character"
         )
)
