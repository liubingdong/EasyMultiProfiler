##' @method as.data.frame multiEnrichResult
##' @export
as.data.frame.multiEnrichResult <- function(x, ...) {
    as.data.frame(x@result, ...)
}

##' @method as.data.frame multiGseaResult
##' @export
as.data.frame.multiGseaResult <- function(x, ...) {
    as.data.frame(x@result, ...)
}

##' @method geneID multiEnrichResult
##' @export
geneID.multiEnrichResult <- function(x) as.character(x@geneID)

##' @method geneID multiGseaResult
##' @export
geneID.multiGseaResult <- function(x) as.character(x@core_enrichment)


##' @method geneInCategory multiEnrichResult
##' @export
##' @importFrom stats setNames
geneInCategory.multiEnrichResult <- function(x) {
    setNames(strsplit(geneID(x), "/", fixed=TRUE), rownames(x@result))
}
    

##' @method geneInCategory multiGseaResult
##' @export
geneInCategory.multiGseaResult <- function(x) {
    setNames(strsplit(geneID(x), "/", fixed=TRUE), rownames(x@result))
}
    


##' @method [ multiEnrichResult
##' @export
`[.multiEnrichResult` <- function(x, i, j, asis = FALSE, ...) {
    y <- x@result[i, j, ...]
    if (!asis)
        return(y)
    x@result <- y
    return(x)
}

##' @method [ multiGseaResult
##' @export
`[.multiGseaResult` <- function(x, i, j, asis = FALSE, ...) {
    y <- x@result[i, j, ...]
    if (!asis)
        return(y)
    x@result <- y
    return(x)
}


##' @method $ multiEnrichResult
##' @export
`$.multiEnrichResult` <-  function(x, name) {
    x@result[, name]
}

##' @method $ multiGseaResult
##' @export
`$.multiGseaResult` <- function(x, name) {
    x@result[, name]
}



##' @method [[ multiEnrichResult
##' @export
`[[.multiEnrichResult` <- function(x, i) {
    gc <- geneInCategory(x)
    if (!i %in% names(gc))
        stop("input term not found...")
    gc[[i]]
}


##' @method [[ multiGseaResult
##' @export
`[[.multiGseaResult` <- function(x, i) {
    gc <- geneInCategory(x)
    if (!i %in% names(gc))
        stop("input term not found...")
    gc[[i]]
}


##' @importFrom utils head
##' @method head multiEnrichResult
##' @export
head.multiEnrichResult <- function(x, n=6L, ...) {
    head(x@result, n, ...)
}

##' @method head multiGseaResult
##' @export
head.multiGseaResult <- function(x, n=6L, ...) {
    head(x@result, n, ...)
}

##' @importFrom utils tail
##' @method tail multiEnrichResult
##' @export
tail.multiEnrichResult <- function(x, n=6L, ...) {
    tail(x@result, n, ...)
}

##' @method tail multiGseaResult
##' @export
tail.multiGseaResult <- function(x, n=6L, ...) {
    tail(x@result, n, ...)
}

##' @method dim multiEnrichResult
##' @export
dim.multiEnrichResult <- function(x) {
    dim(x@result)
}

##' @method dim multiGseaResult
##' @export
dim.multiGseaResult <- function(x) {
    dim(x@result)
}


