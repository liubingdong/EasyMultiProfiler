#' @import tidySummarizedExperiment
mutate.SummarizedExperiment <- getFromNamespace("mutate.SummarizedExperiment", "tidySummarizedExperiment")
#' Title
#'
#' @param x wait_for_add
#' @param experiment wait_for_add
#' @param ... wait_for_add
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_mutate_deprecated <- function(x,experiment=NULL,...) {
  if (inherits(x,"MultiAssayExperiment")) {
    EMPT <- .as.EMPT(x,
                     experiment = experiment)
  }else if(inherits(x,'EMPT')) {
    EMPT <-x
  }

  EMPT %<>% mutate.SummarizedExperiment(...) %>% suppressMessages()
 return(EMPT)
}

filter.SummarizedExperiment <- getFromNamespace("filter.SummarizedExperiment", "tidySummarizedExperiment")

#' Title
#'
#' @param x wait_for_add
#' @param condition wait_for_add
#' @param experiment wait_for_add
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_filter_deprecated <- function(x,condition,experiment=NULL) {
   condition <- dplyr::enquo(condition)
   if (inherits(x,"MultiAssayExperiment")) {
    EMPT <- .as.EMPT(x,
                     experiment = experiment)
  }else if(inherits(x,'EMPT')) {
    EMPT <-x
  }

  EMPT %<>%  filter.SummarizedExperiment(!!condition)
  return(EMPT)
}
