#' geneID generic
#'
#' @param x enrichResult object
#' @return 'geneID' return the 'geneID' column of the enriched result which can be converted to data.frame via 'as.data.frame'
#' @export
geneID <- function(x) {
   UseMethod("geneID", x)
}

#' geneInCategory generic
#'
#' @param x enrichResult
#' @return 'geneInCategory' return a list of genes, by spliting the input gene vector to enriched functional categories
#' @export
geneInCategory <- function(x) {
   UseMethod("geneInCategory", x)
}

