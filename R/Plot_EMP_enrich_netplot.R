#' @param obj EMP or EMPT object
#' @param plot_category An interger. More plot style.(under constrution)
#' @param seed An interger. Set the random seed to the plot.(default:123)
#' @param showCategory A number or a vector of terms. If it is a number, the first n terms will be displayed. If it is a vector of terms, the selected terms will be displayed.
#' @param ... Further parameters passed to enrichplot::cnetplot.
#' @rdname EMP_netplot
#' @return Enrichment netplot object
#' @export
#'
#' @examples
#' # add example
EMP_netplot_enrich <- function(obj,plot_category = 1,seed =123,showCategory=10,...) {
  call <- match.call()
  if (inherits(obj,c('EMP_multi_diff_enrich','EMP_multi_same_enrich'))) {  
    .get.plot_category.EMP(obj) <- plot_category
    .get.history.EMP(obj) <- call
    switch(.get.plot_category.EMP(obj),
         "1" = {
           withr::with_seed(seed,EMP_netplot_enrich_default(obj,showCategory=showCategory,...))
         },
         "2" = {
           # where is EMP_dotplot_enrich_2?
           # withr::with_seed(seed,EMP_dotplot_enrich_2(EMPT,...))
         }

    )    
  }else if(inherits(obj,'EMPT')) {
    .get.plot_category.EMPT(obj) <- plot_category
    .get.history.EMPT(obj) <- call
    switch(.get.plot_category.EMPT(obj),
         "1" = {
           withr::with_seed(seed,EMP_netplot_enrich_default(obj,showCategory=showCategory,...))
         },
         "2" = {
           # where is EMP_dotplot_enrich_2?
           # withr::with_seed(seed,EMP_dotplot_enrich_2(EMPT,...))
         }

    )    
  }else{
    stop("Please input the EMP or EMPT object!")
  } 

}


#' @importFrom enrichplot cnetplot
EMP_netplot_enrich_default <- function(obj,show='pic',showCategory=10,...) {
  enrich_plot <- list()

  if (inherits(obj,"EMP_multi_same_enrich")) {
    p <- obj@deposit[['multi_same_enrich']] %>% 
      cnetplot(showCategory=showCategory,...)

    enrich_plot[['pic']] <- p
    #enrich_plot[['html']] <- plotly::ggplotly(p)

    .get.plot_deposit.EMP(obj,info = 'multi_same_enrich_netplot') <- enrich_plot
    .get.plot_specific.EMP(obj) <- show
    .get.info.EMP(obj) <- 'EMP_multi_same_enrich_netplot'
    class(obj) <- 'EMP_multi_same_enrich_netplot'
  }else if (inherits(obj,"EMP_multi_diff_enrich")) {
    p <- obj@deposit[['multi_diff_enrich']] %>% 
      cnetplot(showCategory=showCategory,...)

    enrich_plot[['pic']] <- p
    #enrich_plot[['html']] <- plotly::ggplotly(p)

    .get.plot_deposit.EMP(obj,info = 'multi_diff_enrich_netplot') <- enrich_plot
    .get.plot_specific.EMP(obj) <- show
    .get.info.EMP(obj) <- 'EMP_multi_diff_enrich_netplot'
    class(obj) <- 'EMP_multi_diff_enrich_netplot'    
  }else if(inherits(obj,'EMPT')) {
    p <- EMP_result(obj,info='enrich_data') %>% 
      cnetplot(showCategory=showCategory,...)

    enrich_plot[['pic']] <- p
    #enrich_plot[['html']] <- plotly::ggplotly(p)

    .get.plot_deposit.EMPT(obj,info = 'enrich_analysis_netplot') <- enrich_plot
    .get.plot_specific.EMPT(obj) <- show
    .get.algorithm.EMPT(obj) <- 'enrich_analysis_netplot'
    .get.info.EMPT(obj) <- 'EMP_enrich_analysis_netplot'
    class(obj) <- 'EMP_enrich_analysis_netplot'
  }else{
    stop("Please input the EMP or EMPT object!")
  }

  return(obj)

}


.show_EMP_netplot_enrich <- function(obj,plot) {
  if (inherits(obj,"EMP_multi_same_enrich")) {
      result <- .get.plot_deposit.EMP(obj,info = 'multi_same_enrich_netplot')
  }else if (inherits(obj,"EMP_multi_diff_enrich")) {
      result <- .get.plot_deposit.EMP(obj,info = 'multi_diff_enrich_netplot')
  }else if (inherits(obj,"EMPT")) {
      result <- .get.plot_deposit.EMPT(obj,info = 'enrich_analysis_netplot')
  }
  switch(plot,
         "pic" = print(result$pic),
         "html" = print(result$html)
  )
}
