get_geneSet_index <- getFromNamespace("get_geneSet_index", "DOSE")
# get_enriched <- getFromNamespace("get_enriched", "DOSE")

#' Title
#'
#' @param object enrichment result.
#' @noRd
get_enriched2 <- function(object) {

    Over <- object@result
    if (class(object) == "multiEnrichResult") {
        pvalueCutoff <- object@pvalueCutoff
        qvalueCutoff <- object@qvalueCutoff
    } else {
        pvalueCutoff <- object@params$pvalueCutoff
        qvalueCutoff <- object@params$qvalueCutoff
    }

    if (length(pvalueCutoff) != 0) {
        ## if groupGO result, numeric(0)
        Over <- Over[ Over$pvalue <= pvalueCutoff, ]
        Over <- Over[ Over$p.adjust <= pvalueCutoff, ]
    }


    if (length(qvalueCutoff) != 0) {
        if (! any(is.na(Over$qvalue))) {
            if (length(qvalueCutoff) > 0)
                Over <- Over[ Over$qvalue <= qvalueCutoff, ]
        }
    }

    object@result <- Over
    return(object)
}
