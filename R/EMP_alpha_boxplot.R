#' Title
#'
#' @param EMPT wait_for_add
#' @param plot_category wait_for_add
#' @param seed wait_for_add
#' @param ... wait_for_add
#' @importFrom withr with_seed
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_boxplot_alpha <- function(EMPT,plot_category = 1,seed =123,...) {
  #call <- match.call()
  .get.plot_category.EMPT(EMPT) <- plot_category
  #.get.history.EMPT(EMPT) <- call

  switch(.get.plot_category.EMPT(EMPT),
         "1" = {
           withr::with_seed(seed,EMP_boxplot_alpha_default(EMPT,...))
         },
         "2" = {
           withr::with_seed(seed,EMP_boxplot_alpha_2(EMPT,...))
         }

  )

}


.show_EMP_alpha_boxplot <- function(obj,plot) {
  result <- .get.plot_deposit.EMPT(obj,info = 'alpha_analysis_plot')
  switch(plot,
         "pic" = print(result$pic),
         "html" = print(result$html)
  )
}



#' Title
#'
#' @param EMPT wait_for_add
#' @param method wait_for_add
#' @param estimate_group wait_for_add
#' @param group_level wait_for_add
#' @param ncol wait_for_add
#' @param select_metrics wait_for_add
#' @param palette wait_for_add
#' @param show wait_for_add
#' @param html_width wait_for_add
#' @param html_height wait_for_add
#' @param mytheme wait_for_add
#' @importFrom tidyr pivot_longer
#' @importFrom plyr alply
#' @importFrom dplyr sym
#' @importFrom ggiraph geom_jitter_interactive
#' @importFrom ggsignif geom_signif
#' @importFrom dplyr full_join
#' @importFrom plyr alply
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr filter
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_boxplot_alpha_default <- function (EMPT,method = 'wilcox.test',
                                       estimate_group = NULL,group_level = 'default',
                                       ncol = NULL,select_metrics = NULL,palette = NULL,
                                 show = 'pic',html_width=NULL,html_height=NULL,mytheme = 'theme()') {

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

  compare <- plyr::alply(group_combn,2)

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


  alpha_plot[['pic']] <- alpha_data %>%
    ggplot(., aes(x = !!dplyr::sym(estimate_group), y = value, fill = !!dplyr::sym(estimate_group))) +
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




##### only for test
#' Title
#'
#' @param EMPT wait_for_add
#' @param method wait_for_add
#' @param estimate_group wait_for_add
#' @param shape wait_for_add
#' @param show wait_for_add
#' @param html_width wait_for_add
#' @param html_height wait_for_add
#' @param mytheme wait_for_add
#' @noRd
EMP_boxplot_alpha_2 <- function(EMPT,method = 'wilcox.test',estimate_group = NULL,shape,
                                       show = 'pic',html_width=NULL,html_height=NULL,mytheme = 'theme()') {

  call <- match.call()
  alpha_plot <- list()

  estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)


  mapping <- .get.mapping.EMPT(EMPT) %>% dplyr::select(primary,!!estimate_group)

  alpha_data <- .get.result.EMPT(EMPT) %>%
    dplyr::full_join(mapping,by= 'primary')

  group_combn <- combn(as.character(unique(mapping[[estimate_group]])),2)
  compare <- plyr::alply(group_combn,2)

  col_values <- .get.palette.EMPT(EMPT)
  alpha_data %<>%  tidyr::pivot_longer(cols = c(-primary,-!!dplyr::sym(estimate_group)),
                                       names_to = 'ID',
                                       values_to = 'value')

  alpha_plot[['pic']] <- alpha_data %>%
    ggplot(., aes(x = !!dplyr::sym(estimate_group), y = value, fill = !!dplyr::sym(estimate_group))) +
    geom_boxplot(outlier.color=NA,shape = shape,alpha=0.5) +
    ggiraph::geom_jitter_interactive(aes(tooltip = value),shape=21,position = position_jitter(height = .00000001))+
    ggsignif::geom_signif(comparisons = compare,test = method,step_increase = 0.1) +
    facet_wrap(ID~., scales = 'free', strip.position = 'right') +
    xlab(NULL) +
    ylab("Alpha Metrics") +
    ggtitle('Alpha analysis Plot') +
    scale_fill_manual(values = col_values) +
    theme_bw() + eval(parse(text = paste0(mytheme))) +theme_few()


  alpha_plot[['html']] <- ggiraph::girafe(code = print(alpha_plot[['pic']]),width = html_width,height = html_height)

  .get.plot_deposit.EMPT(EMPT,info = 'alpha_analysis_plot') <- alpha_plot
  .get.plot_specific.EMPT(EMPT) <- show
  .get.estimate_group.EMPT(EMPT) <- estimate_group
  .get.algorithm.EMPT(EMPT) <- 'alpha_analysis_plot'
  .get.info.EMPT(EMPT) <- 'EMP_alpha_analysis_boxplot'
  .get.history.EMPT(EMPT) <- call
  class(EMPT) <- 'EMP_alpha_analysis_boxplot'
  EMPT
}





