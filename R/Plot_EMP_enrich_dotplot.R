#' @param obj EMPT or EMP object
#' @param plot_category An interger.More plot style.(under constrution)
#' @param seed An interger. Set the random seed to the plot.(default:123)
#' @param color A character string. Variable that used to color enriched terms, e.g. 'pvalue', 'p.adjust' or 'qvalue'
#' @param showCategory A number or a list of terms. If it is a number, the first n terms will be displayed. If it is a list of terms, the selected terms will be displayed.
#' @param ... Further parameters passed to enrichplot::dotplot.
#' @rdname EMP_dotplot
#'
#' @examples
#' # add example
EMP_dotplot_enrich <- function(obj,plot_category = 1,seed =123,color='p.adjust',showCategory=10,...) {
  call <- match.call()
  if (is(obj,'EMP_multi_diff_enrich') | is(obj,'EMP_multi_same_enrich')) {  
    .get.plot_category.EMP(obj) <- plot_category
    .get.history.EMP(obj) <- call
    switch(.get.plot_category.EMP(obj),
         "1" = {
           withr::with_seed(seed,EMP_dotplot_enrich_default(obj,color=color,showCategory=showCategory,...))
         },
         "2" = {
           # where is EMP_dotplot_enrich_2?
           # withr::with_seed(seed,EMP_dotplot_enrich_2(EMPT,...))
         }

    )    
  }else if(is(obj,'EMPT')) {
    .get.plot_category.EMPT(obj) <- plot_category
    .get.history.EMPT(obj) <- call
    switch(.get.plot_category.EMPT(obj),
         "1" = {
           withr::with_seed(seed,EMP_dotplot_enrich_default(obj,color=color,showCategory=showCategory,...))
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


#' @importFrom enrichplot dotplot
EMP_dotplot_enrich_default <- function(obj,show='pic',color='p.adjust',showCategory=10,...) {
  enrich_plot <- list()

  if (is(obj,"EMP_multi_same_enrich")) {
    p <- obj@deposit[['multi_same_enrich']] %>% 
      dotplot(color=color,showCategory=showCategory,...)

    enrich_plot[['pic']] <- p
    #enrich_plot[['html']] <- plotly::ggplotly(p)

    .get.plot_deposit.EMP(obj,info = 'multi_same_enrich_dotplot') <- enrich_plot
    .get.plot_specific.EMP(obj) <- show
    .get.info.EMP(obj) <- 'EMP_multi_same_enrich_dotplot'
    class(obj) <- 'EMP_multi_same_enrich_dotplot'
  }else if (is(obj,"EMP_multi_diff_enrich")) {
    p <- obj@deposit[['multi_diff_enrich']] %>% 
      dotplot(color=color,showCategory=showCategory,...)

    enrich_plot[['pic']] <- p
    #enrich_plot[['html']] <- plotly::ggplotly(p)

    .get.plot_deposit.EMP(obj,info = 'multi_diff_enrich_dotplot') <- enrich_plot
    .get.plot_specific.EMP(obj) <- show
    .get.info.EMP(obj) <- 'EMP_multi_diff_enrich_dotplot'
    class(obj) <- 'EMP_multi_diff_enrich_dotplot'    
  }else if(is(obj,'EMPT')) {
    p <- EMP_result(obj,info='enrich_data') %>% 
      dotplot(color=color,showCategory=showCategory,...)

    enrich_plot[['pic']] <- p
    #enrich_plot[['html']] <- plotly::ggplotly(p)

    .get.plot_deposit.EMPT(obj,info = 'enrich_analysis_dotplot') <- enrich_plot
    .get.plot_specific.EMPT(obj) <- show
    .get.algorithm.EMPT(obj) <- 'enrich_analysis_dotplot'
    .get.info.EMPT(obj) <- 'EMP_enrich_analysis_dotplot'
    class(obj) <- 'EMP_enrich_analysis_dotplot'
  }else{
    stop("Please input the EMP or EMPT object!")
  }

  return(obj)

}

.show_EMP_dotplot_enrich <- function(obj,plot) {
  if (is(obj,"EMP_multi_same_enrich")) {
      result <- .get.plot_deposit.EMP(obj,info = 'multi_same_enrich_dotplot')
  }else if (is(obj,"EMP_multi_diff_enrich")) {
      result <- .get.plot_deposit.EMP(obj,info = 'multi_diff_enrich_dotplot')
  }else if (is(obj,"EMPT")) {
      result <- .get.plot_deposit.EMPT(obj,info = 'enrich_analysis_dotplot')
  }
  switch(plot,
         "pic" = print(result$pic),
         "html" = print(result$html)
  )
}

