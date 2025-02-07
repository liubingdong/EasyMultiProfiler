
#' @importFrom enrichplot gseaplot2
EMP_curveplot_enrich_default <- function(EMPT,geneSetID,show='pic',...) {
  enrich_plot <- list()
  if(is.null(geneSetID)){
    stop('Please select a result generated from EMP_GSEA_analysis!')

  }
  p <- .get.result.EMPT(EMPT) %>%
    gseaplot2(geneSetID,...)

  enrich_plot[['pic']] <- p[[1]] + p[[2]] + p[[3]] +  patchwork::plot_layout(ncol = 1,nrow = 3,heights = c(4.5,1.5,3))
  #enrich_plot[['html']] <- plotly::ggplotly(p)

  #EMPT@deposit[['gsea_curve_plot']] <- enrich_plot
  .get.plot_deposit.EMPT(EMPT,info = 'gsea_curve_plot') <- enrich_plot
  .get.plot_specific.EMPT(EMPT) <- show
  .get.algorithm.EMPT(EMPT) <- 'gsea_curve_plot'
  .get.info.EMPT(EMPT) <- 'EMP_enrich_analysis_curveplot'
  class(EMPT) <- 'EMP_enrich_analysis_curveplot'
  EMPT
}


#' @param obj Object in EMPT format.
#' @param plot_category An interger.More plot style.[under constrution]
#' @param geneSetID geneSet ID
#' @param seed An interger. Set the random seed to the plot.(default:123)
#' @param ... Additional parameters, see also \code{\link[enrichplot]{gseaplot2}}
#' @rdname EMP_gsea_plot
#' @importFrom withr with_seed
#'
#' @return EMPT object
#'
#' @examples
#' # add example
EMP_curveplot_enrich <- function(obj,plot_category = 1,seed =123,geneSetID,...) {
  call <- match.call()
  .get.plot_category.EMPT(obj) <- plot_category
  .get.history.EMPT(obj) <- call

  switch(.get.plot_category.EMPT(obj),
         "1" = {
           withr::with_seed(seed,EMP_curveplot_enrich_default(EMPT=obj,geneSetID,...))
         },
         "2" = {
           # withr::with_seed(seed,EMP_curveplot_enrich_2(EMPT,...))
         }

  )

}

.show_EMP_curveplot_enrich <- function(obj,plot) {
  switch(plot,
         "pic" = print(obj@plot_deposit$gsea_curve_plot$pic),
         "html" = print(obj@plot_deposit$gsea_curve_plot$html)
  )
}

