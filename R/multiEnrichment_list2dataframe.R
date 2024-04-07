#' Transform the result data of gene differential expression from list to dataframe
#'
#' @param genelists Results of gene difference analysis.
#' @param method The method of merging, one of "intersection" and "union"
#' @export
#'
#' @examples
#' \dontrun{
#'     library(multiGSEA)
#'     data(transcriptome)
#'     data(proteome)
#'     data(metabolome) 
#'     class(transcriptome) <- class(proteome) <- class(metabolome) <- "data.frame"
#'     transcriptome2 <- transcriptome$pValue
#'     names(transcriptome2) <- transcriptome$Symbol
#'     proteome2 <- proteome$pValue
#'     names(proteome2) <- proteome$Symbol
#'     proteome2 <- proteome2[-which(names(proteome2)== "")]
#'     genelists <- list(transcriptome2 = transcriptome2, proteome2 = proteome2)
#'     geneDf <- list2dataframe(genelists)
#' }
list2dataframe <- function(genelists, method = "intersection") {
    method <- match.arg(method, c("intersection", "union"))
    genes <- unique(unlist(lapply(genelists, names)))
    geneDf <- matrix(NA, nrow  = length(genes), ncol = length(genelists))
    rownames(geneDf) <- genes
    for (i in seq_len(length(genelists))) {
        geneDf[names(genelists[[i]]), i] <- genelists[[i]]
    }
    if (method == "intersection") {
        geneDf <- stats::na.omit(geneDf)
    }
    return(geneDf)
}
