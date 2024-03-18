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
EMP_dotplot_enrich <- function(EMPT,plot_category = 1,seed =123,...) {
  call <- match.call()
  .get.plot_category.EMPT(EMPT) <- plot_category
  .get.history.EMPT(EMPT) <- call

  switch(.get.plot_category.EMPT(EMPT),
         "1" = {
           withr::with_seed(seed,EMP_dotplot_enrich_default(EMPT,...))
         },
         "2" = {
           # where is EMP_dotplot_enrich_2?
           # withr::with_seed(seed,EMP_dotplot_enrich_2(EMPT,...))
         }

  )

}

#' Title
#'
#' @param EMPT wait_for_add
#' @param show wait_for_add
#' @param ... wait_for_add
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_dotplot_enrich_default <- function(EMPT,show='pic',...) {
  enrich_plot <- list()
  p <- .get.result.EMPT(EMPT) %>%
    enrichplot::dotplot(...)

  enrich_plot[['pic']] <- p
  enrich_plot[['html']] <- plotly::ggplotly(p)

  .get.plot_deposit.EMPT(EMPT,info = 'enrich_analysis_dotplot') <- enrich_plot
  .get.plot_specific.EMPT(EMPT) <- show
  .get.algorithm.EMPT(EMPT) <- 'enrich_analysis_dotplot'
  .get.info.EMPT(EMPT) <- 'EMP_enrich_analysis_dotplot'
  class(EMPT) <- 'EMP_enrich_analysis_dotplot'
  EMPT
}

.show_EMP_dotplot_enrich <- function(obj,plot) {
  result <- .get.plot_deposit.EMPT(obj,info = 'enrich_analysis_dotplot')
  switch(plot,
         "pic" = print(result$pic),
         "html" = print(result$html)
  )
}

