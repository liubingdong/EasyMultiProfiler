#' @param EMPT EMPT object
#' @param plot_category An interger.More plot style.(under constrution)
#' @param seed An interger. Set the random seed to the plot.(default:123)
#' @param color A character string. Variable that used to color enriched terms, e.g. 'pvalue', 'p.adjust' or 'qvalue'
#' @param showCategory A number or a list of terms. If it is a number, the first n terms will be displayed. If it is a list of terms, the selected terms will be displayed.
#' @param ... Further parameters passed to enrichplot::dotplot.
#' @rdname EMP_dotplot
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_dotplot_enrich <- function(EMPT,plot_category = 1,seed =123,color='p.adjust',showCategory=10,...) {
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


#' @importFrom enrichplot dotplot
EMP_dotplot_enrich_default <- function(EMPT,show='pic',color='p.adjust',showCategory=10,...) {
  enrich_plot <- list()
  p <- .get.result.EMPT(EMPT) %>% 
    dotplot(color=color,showCategory=showCategory,...)

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

