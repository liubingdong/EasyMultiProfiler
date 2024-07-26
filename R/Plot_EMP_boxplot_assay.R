#' @import ggthemes
#' @import ggplot2
#' @importFrom ggiraph geom_jitter_interactive
#' @importFrom ggsignif geom_signif
EMP_boxplot_assay_default <- function (EMPT,method = 'wilcox.test',
                               estimate_group = NULL,group_level = 'default',
                               ncol = NULL,show = 'pic',palette = NULL,
                               html_width=NULL,html_height=NULL,
                               mytheme = 'theme()') {

  primary <- value <- `.` <- NULL
  estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)

  if (is.null(palette)) {
    col_values <- .get.palette.EMPT(EMPT)
  }else {
    col_values = palette
  }

  EMPT %<>% .group_level_modified(estimate_group = estimate_group,
                                  group_level = group_level)

  mapping <- .get.mapping.EMPT(EMPT) %>% dplyr::select(primary,!!estimate_group)

  data <-.get.result.EMPT(EMPT,info = 'EMP_assay_data') %>% dplyr::left_join(mapping,by ='primary')

  ## clean the missing value in the group label
  if(any(is.na(data[[estimate_group]]))) {
    warning('Column ',estimate_group,' has beed deteced missing value, all related samples will be removed in the display!')
    data <- data %>% tidyr::drop_na(!!estimate_group)
  }

  group_combn <- combn(as.character(unique(mapping[[estimate_group]])),2)
  #compare <- plyr::alply(group_combn,2)
  compare <- list() 
  for (i in 1:ncol(group_combn)) {
    compare[[i]] <- group_combn[,i]
  }
  names(compare) <- 1:ncol(group_combn)

  data %<>%  tidyr::pivot_longer(cols = c(-primary,-!!dplyr::sym(estimate_group)),
                                       names_to = 'ID',
                                       values_to = 'value')

  data_plot <- list()
  data_plot[['pic']] <- data %>%
    ggplot(., aes(x = !!dplyr::sym(estimate_group), y = value, fill = !!dplyr::sym(estimate_group))) +
    geom_boxplot(outlier.color=NA) +
    ggiraph::geom_jitter_interactive(aes(tooltip = paste0(primary,' : ',value)),shape=21,position = position_jitter(height = .00000001))+
    ggsignif::geom_signif(comparisons = compare,test = method,step_increase = 0.1) +
    facet_wrap(ID~., scales = 'free', strip.position = 'top',ncol = ncol) +
    xlab(NULL) + 
    ggtitle('Feature Boxplot') +
    scale_fill_manual(values = col_values) +
    theme_bw() + 
    theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
    eval(parse(text = paste0(mytheme)))


  data_plot[['html']] <- ggiraph::girafe(code = print(data_plot[['pic']]),width = html_width,height = html_height)

  .get.plot_deposit.EMPT(EMPT,info = 'EMP_assay_boxplot') <- data_plot
  .get.plot_specific.EMPT(EMPT) <- show
  .get.estimate_group.EMPT(EMPT) <- estimate_group
  EMPT@algorithm <- 'EMP_assay_boxplot'
  .get.info.EMPT(EMPT) <- 'EMP_assay_boxplot'
  class(EMPT) <- 'EMP_assay_boxplot'
  EMPT
}

.show_EMP_assay_boxplot <- function(obj,plot) {
  result <- .get.plot_deposit.EMPT(obj,info = 'EMP_assay_boxplot')
  switch(plot,
         "pic" = print(result$pic),
         "html" = print(result$html)
  )
}


#' @param obj EMPT object
#' @param plot_category An interger.More plot style.(under constrution)
#' @param seed An interger. Set the random seed to the plot.
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

EMP_boxplot.EMP_assay_data <- function(obj,plot_category = 1,seed =123,method = 'wilcox.test',
                               estimate_group = NULL,group_level = 'default',
                               ncol = NULL,show = 'pic',palette = NULL,
                               html_width=NULL,html_height=NULL,
                               mytheme = 'theme()') {
  call <- match.call()
  .get.plot_category.EMPT(obj) <- plot_category
  .get.history.EMPT(obj) <- call
  switch(.get.plot_category.EMPT(obj),
         "1" = {
           withr::with_seed(seed,EMP_boxplot_assay_default(EMPT=obj,method = method,
                               estimate_group = estimate_group,group_level = group_level,
                               ncol = ncol,show = show,palette = palette,
                               html_width=html_width,html_height=html_height,
                               mytheme = mytheme))
         },
         "2" = {
           # where is EMP_boxplot_assay_2?
           # withr::with_seed(seed,EMP_boxplot_assay_2(EMPT,...))
         }

  )

}

EMP_boxplot.EMP_decostand <- function(obj,plot_category = 1,seed =123,method = 'wilcox.test',
                               estimate_group = NULL,group_level = 'default',
                               ncol = NULL,show = 'pic',palette = NULL,
                               html_width=NULL,html_height=NULL,
                               mytheme = 'theme()') {
  call <- match.call()
  .get.plot_category.EMPT(obj) <- plot_category
  .get.history.EMPT(obj) <- call
  switch(.get.plot_category.EMPT(obj),
         "1" = {
           withr::with_seed(seed,EMP_boxplot_assay_default(EMPT=obj,method = method,
                               estimate_group = estimate_group,group_level = group_level,
                               ncol = ncol,show = show,palette = palette,
                               html_width=html_width,html_height=html_height,
                               mytheme = mytheme))
         },
         "2" = {
           # where is EMP_boxplot_assay_2?
           # withr::with_seed(seed,EMP_boxplot_assay_2(EMPT,...))
         }

  )

}

