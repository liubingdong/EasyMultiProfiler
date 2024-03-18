#' Title
#'
#' @param EMPT wait_for_add
#' @param plot_category wait_for_add
#' @param seed wait_for_add
#' @param ... wait_for_add
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_netplot_enrich <- function(EMPT,plot_category = 1,seed =123,...) {
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

#' Title
#'
#' @param EMPT wait_for_add
#' @param show wait_for_add
#' @param ... wait_for_add
#' @importFrom enrichplot cnetplot
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_netplot_enrich_default <- function(EMPT,show='pic',...) {
  enrich_plot <- list()
  p <- .get.result.EMPT(EMPT) %>%
    enrichplot::cnetplot(...)

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

