#' Title
#'
#' @param EMPT wait_for_add
#' @param geneSetID wait_for_add
#' @param show wait_for_add
#' @param ... wait_for_add
#' @importFrom enrichplot gseaplot2
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_curveplot_enrich_default <- function(EMPT,geneSetID,show='pic',...) {
  enrich_plot <- list()
  if(is.null(geneSetID)){
    stop('Please select a result generated from EMP_GSEA_analysis!')

  }
  p <- .get.result.EMPT(EMPT) %>%
    enrichplot::gseaplot2(geneSetID,...)

  enrich_plot[['pic']] <- p
  #enrich_plot[['html']] <- plotly::ggplotly(p)

  EMPT@deposit[['enrich_analysis_plot']] <- enrich_plot
  .get.plot_specific.EMPT(EMPT) <- show
  .get.algorithm.EMPT(EMPT) <- 'enrich_analysis_plot'
  .get.info.EMPT(EMPT) <- 'EMP_enrich_analysis_curveplot'
  class(EMPT) <- 'EMP_enrich_analysis_curveplot'
  EMPT
}

#' Title
#'
#' @param EMPT wait_for_add
#' @param plot_category wait_for_add
#' @param seed wait_for_add
#' @param ... wait_for_add
#' @importFrom withr with_seed
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_curveplot_enrich <- function(EMPT,plot_category = 1,seed =123,...) {
  call <- match.call()
  .get.plot_category.EMPT(EMPT) <- plot_category
  .get.history.EMPT(EMPT) <- call

  switch(.get.plot_category.EMPT(EMPT),
         "1" = {
           withr::with_seed(seed,EMP_curveplot_enrich_default(EMPT,...))
         },
         "2" = {
           withr::with_seed(seed,EMP_curveplot_enrich_2(EMPT,...))
         }

  )

}

.show_EMP_curveplot_enrich <- function(obj,plot) {
  switch(plot,
         "pic" = print(obj@deposit$enrich_analysis_plot$pic),
         "html" = print(obj@deposit$enrich_analysis_plot$html)
  )
}

