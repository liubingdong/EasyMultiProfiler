#' Title
#'
#' @param EMPT wait_for_add
#' @param method wait_for_add
#' @param estimate_group wait_for_add
#' @param group_level  wait_for_add
#' @param ncol wait_for_add
#' @param show wait_for_add
#' @param palette wait_for_add
#' @param html_width wait_for_add
#' @param html_height wait_for_add
#' @param mytheme wait_for_add
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_boxplot_assay_default <- function (EMPT,method = 'wilcox.test',
                               estimate_group = NULL,group_level = 'default',
                               ncol = NULL,show = 'pic',palette = NULL,
                               html_width=NULL,html_height=NULL,
                               mytheme = 'theme()') {


  estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)

  if (is.null(palette)) {
    col_values <- .get.palette.EMPT(EMPT)
  }else {
    col_values = palette
  }

  EMPT %<>% .group_level_modified(estimate_group = estimate_group,
                                  group_level = group_level)

  mapping <- .get.mapping.EMPT(EMPT) %>% dplyr::select(primary,!!estimate_group)

  data <-.get.result.EMPT(EMPT) %>% dplyr::left_join(mapping,by ='primary')

  group_combn <- combn(as.character(unique(mapping[[estimate_group]])),2)
  compare <- plyr::alply(group_combn,2)


  data %<>%  tidyr::pivot_longer(cols = c(-primary,-!!dplyr::sym(estimate_group)),
                                       names_to = 'ID',
                                       values_to = 'value')

  data_plot <- list()
  data_plot[['pic']] <- data %>%
    ggplot(., aes(x = !!dplyr::sym(estimate_group), y = value, fill = !!dplyr::sym(estimate_group))) +
    geom_boxplot(outlier.color=NA) +
    ggiraph::geom_jitter_interactive(aes(tooltip = paste0(primary,' : ',value)),shape=21,position = position_jitter(height = .00000001))+
    ggsignif::geom_signif(comparisons = compare,test = method,step_increase = 0.1) +
    facet_wrap(ID~., scales = 'free', strip.position = 'right',ncol = ncol) +
    xlab(NULL) +
    ggtitle('Feature Boxplot') +
    scale_fill_manual(values = col_values) +
    theme_bw() + eval(parse(text = paste0(mytheme)))


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
EMP_assay_boxplot <- function(EMPT,plot_category = 1,seed =123,...) {
  #call <- match.call()
  .get.plot_category.EMPT(EMPT) <- plot_category
  #.get.history.EMPT(EMPT) <- call
  switch(.get.plot_category.EMPT(EMPT),
         "1" = {
           withr::with_seed(seed,EMP_boxplot_assay_default(EMPT,...))
         },
         "2" = {
           withr::with_seed(seed,EMP_boxplot_assay_2(EMPT,...))
         }

  )

}



