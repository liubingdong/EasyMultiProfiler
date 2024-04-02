#' @param EMPT Object in EMPT or MultiAssayExperiment format.
#' @param plot_category An interger.More plot style.(under constrution)
#' @param seed An interger. Set the random seed to the plot.(default:123)
#' @param obj EMPT object
#' @param method A character string. The name of the statistical test that is applied to the values of the columns (e.g. t.test, wilcox.test etc.).
#' @param estimate_group A character string. Select the colname in the coldata to compare the data in the statistical test.
#' @param group_level A string vector. Set the group order in the plot.
#' @param ncol An interger. Set the col number in the facet plot.
#' @param select_metrics A series of character string. Select the alpha metrics to show. 
#' @param show A character string include pic (default), html.
#' @param palette A series of character string. Color palette.
#' @param html_width An interger. Set the html width.
#' @param html_height An interger. Set the html height.
#' @param mytheme Modify components of a theme according to the ggplot2::theme.
#' @rdname EMP_boxplot
#' @importFrom withr with_seed
#'
#' @export
#'
#' @examples
#' # add example
EMP_boxplot.EMP_alpha_analysis <- function(EMPT,plot_category = 1,seed =123,method = 'wilcox.test',
                               estimate_group = NULL,group_level = 'default',
                               ncol = NULL,select_metrics=NULL,show = 'pic',palette = NULL,
                               html_width=NULL,html_height=NULL,
                               mytheme = 'theme()') {
  #call <- match.call()
  .get.plot_category.EMPT(EMPT) <- plot_category
  #.get.history.EMPT(EMPT) <- call

  switch(.get.plot_category.EMPT(EMPT),
         "1" = {
           withr::with_seed(seed,EMP_boxplot_alpha_default(EMPT,method = method,
                               estimate_group = estimate_group,group_level = group_level,
                               ncol = ncol,select_metrics=select_metrics,show = show,palette = palette,
                               html_width=html_width,html_height=html_height,
                               mytheme = mytheme))
         },
         #"2" = {
         # withr::with_seed(seed,EMP_boxplot_alpha_2(EMPT,method = method,
         #                      estimate_group = estimate_group,group_level = group_level,
         #                      ncol = ncol,select_metrics=select_metrics,show = show,palette = palette,
         #                      html_width=html_width,html_height=html_height,
         #                      mytheme = mytheme))
         #}

  )

}


.show_EMP_alpha_boxplot <- function(obj,plot) {
  result <- .get.plot_deposit.EMPT(obj,info = 'alpha_analysis_plot')
  switch(plot,
         "pic" = print(result$pic),
         "html" = print(result$html)
  )
}



#' @import ggthemes
EMP_boxplot_alpha_default <- function (EMPT,method = 'wilcox.test',
                                       estimate_group = NULL,group_level = 'default',
                                       ncol = NULL,select_metrics = NULL,palette = NULL,
                                 show = 'pic',html_width=NULL,html_height=NULL,mytheme = 'theme()') {
  primary <- ID <- value <- NULL

  alpha_plot <- list()

  estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)

  if (is.null(palette)) {
    col_values <- .get.palette.EMPT(EMPT)
  }else {
    col_values = palette
  }

  EMPT %<>% .group_level_modified(estimate_group = estimate_group,group_level = group_level)

  mapping <- .get.mapping.EMPT(EMPT) %>% dplyr::select(primary,!!estimate_group)

  alpha_data <- .get.result.EMPT(EMPT) %>%
    dplyr::full_join(mapping,by= 'primary')

  group_combn <- combn(as.character(unique(mapping[[estimate_group]])),2)

  #compare <- plyr::alply(group_combn,2)
  compare <- list() 
  for (i in 1:ncol(group_combn)) {
    compare[[i]] <- group_combn[,i]
  }
  names(compare) <- 1:ncol(group_combn)  

  alpha_data %<>%  tidyr::pivot_longer(cols = c(-primary,-!!dplyr::sym(estimate_group)),
                                       names_to = 'ID',
                                       values_to = 'value')

  if (!is.null(select_metrics)) {
    if (all(select_metrics %in% unique(alpha_data[['ID']]))) {
      alpha_data %<>% dplyr::filter(ID %in% select_metrics)
    } else{
      stop('Please check the select_metrics parameter')
    }
  }


  alpha_plot[['pic']] <- ggplot(alpha_data, aes(x = !!dplyr::sym(estimate_group), y = value, fill = !!dplyr::sym(estimate_group))) +
    geom_boxplot(outlier.color=NA) +
    ggiraph::geom_jitter_interactive(aes(tooltip = paste0(primary,' : ',value)),shape=21,position = position_jitter(height = .00000001))+
    ggsignif::geom_signif(comparisons = compare,test = method,step_increase = 0.1) +
    facet_wrap(ID~., scales = 'free', strip.position = 'right',ncol = ncol) +
    xlab(NULL) +
    ylab("Alpha Metrics") +
    ggtitle('Alpha analysis Plot') +
    scale_fill_manual(values = col_values) +
    theme_bw() + eval(parse(text = paste0(mytheme)))


  alpha_plot[['html']] <- ggiraph::girafe(code = print(alpha_plot[['pic']]),width = html_width,height = html_height)

  .get.plot_deposit.EMPT(EMPT,info = 'alpha_analysis_plot') <- alpha_plot
  .get.plot_specific.EMPT(EMPT) <- show
  .get.estimate_group.EMPT(EMPT) <- estimate_group
  .get.algorithm.EMPT(EMPT) <- 'alpha_analysis_plot'
  .get.info.EMPT(EMPT) <- 'EMP_alpha_analysis_boxplot'
  class(EMPT) <- 'EMP_alpha_analysis_boxplot'
  EMPT
}


# EMP_boxplot_alpha_2 <- function(EMPT,method = 'wilcox.test',estimate_group = NULL,shape=21,
#                                 show = 'pic',html_width=NULL,html_height=NULL,mytheme = 'theme()') {
#   primary <- value <- NULL
#   call <- match.call()
#   alpha_plot <- list()
#   
#   estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)
#   
#   
#   mapping <- .get.mapping.EMPT(EMPT) %>% dplyr::select(primary,!!estimate_group)
#   
#   alpha_data <- .get.result.EMPT(EMPT) %>%
#     dplyr::full_join(mapping,by= 'primary')
#   
#   group_combn <- combn(as.character(unique(mapping[[estimate_group]])),2)
#   compare <- plyr::alply(group_combn,2)
#   
#   col_values <- .get.palette.EMPT(EMPT)
#   alpha_data %<>%  tidyr::pivot_longer(cols = c(-primary,-!!dplyr::sym(estimate_group)),
#                                        names_to = 'ID',
#                                        values_to = 'value')
#   
#   alpha_plot[['pic']] <- ggplot(alpha_data, aes(x = !!dplyr::sym(estimate_group), y = value, fill = !!dplyr::sym(estimate_group))) +
#     geom_boxplot(outlier.color=NA,shape = shape,alpha=0.5) +
#     ggiraph::geom_jitter_interactive(aes(tooltip = value),shape=21,position = position_jitter(height = .00000001))+
#     ggsignif::geom_signif(comparisons = compare,test = method,step_increase = 0.1) +
#     facet_wrap(ID~., scales = 'free', strip.position = 'right') +
#     xlab(NULL) +
#     ylab("Alpha Metrics") +
#     ggtitle('Alpha analysis Plot') +
#     scale_fill_manual(values = col_values) +
#     theme_bw() + eval(parse(text = paste0(mytheme))) +theme_few()
#   
#   
#   alpha_plot[['html']] <- ggiraph::girafe(code = print(alpha_plot[['pic']]),width = html_width,height = html_height)
#   
#   .get.plot_deposit.EMPT(EMPT,info = 'alpha_analysis_plot') <- alpha_plot
#   .get.plot_specific.EMPT(EMPT) <- show
#   .get.estimate_group.EMPT(EMPT) <- estimate_group
#   .get.algorithm.EMPT(EMPT) <- 'alpha_analysis_plot'
#   .get.info.EMPT(EMPT) <- 'EMP_alpha_analysis_boxplot'
#   .get.history.EMPT(EMPT) <- call
#   class(EMPT) <- 'EMP_alpha_analysis_boxplot'
#   EMPT
# }





