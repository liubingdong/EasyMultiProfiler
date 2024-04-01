#' Netplot for enrichment result
#'
#' @param EMPT EMPT object
#' @param plot_category wait_for_add
#' @param seed An interger. Set the random seed to the plot.(default:123)
#' @param layout Layout of the map, e.g. 'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' or 'lgl'.
#' @param showCategory A number or a vector of terms. If it is a number, the first n terms will be displayed. If it is a vector of terms, the selected terms will be displayed.
#' @param ... Further parameters passed to enrichplot::cnetplot.
#' @rdname EMP_netplot
#' @return EMPT object
#' @export
#'
#' @examples
#' # add example
EMP_netplot_enrich <- function(EMPT,plot_category = 1,seed =123,layout='kk',showCategory=5,...) {
  call <- match.call()
  .get.plot_category.EMPT(EMPT) <- plot_category
  .get.history.EMPT(EMPT) <- call

  switch(.get.plot_category.EMPT(EMPT),
         "1" = {
           withr::with_seed(seed,EMP_netplot_enrich_default(EMPT,...))
         },
         "2" = {
           # where is EMP_netplot_enrich_2?
           # withr::with_seed(seed,EMP_netplot_enrich_2(EMPT,...))
         }

  )

}

#' @importFrom enrichplot cnetplot

EMP_netplot_enrich_default <- function(EMPT,show='pic',layout='kk',showCategory=5,...) {
  enrich_plot <- list()
  p <- .get.result.EMPT(EMPT) %>%
    cnetplot(layout=layout,showCategory=showCategory,...)

  enrich_plot[['pic']] <- p
  #enrich_plot[['html']] <- plotly::ggplotly(p)

  .get.plot_deposit.EMPT(EMPT,info = 'enrich_analysis_netplot') <- enrich_plot
  .get.plot_specific.EMPT(EMPT) <- show
  .get.algorithm.EMPT(EMPT) <- 'enrich_analysis_netplot'
  .get.info.EMPT(EMPT) <- 'EMP_enrich_analysis_netplot'
  class(EMPT) <- 'EMP_enrich_analysis_netplot'
  EMPT
}

.show_EMP_netplot_enrich <- function(obj,plot) {
  result <- .get.plot_deposit.EMPT(obj,info = 'enrich_analysis_netplot')
  switch(plot,
         "pic" = print(result$pic),
         "html" = print(result$html)
  )
}

