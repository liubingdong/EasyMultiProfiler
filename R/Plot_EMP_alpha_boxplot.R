#' @param obj EMPT object
#' @param plot_category A character string including default and violin.
#' @param method A character string. The name of the statistical test that is applied to the values of the columns (e.g. t.test, wilcox.test etc.).
#' @param estimate_group A character string. Select the colname in the coldata to compare the data in the statistical test.
#' @param group_level A string vector. Set the group order in the plot.
#' @param step_increase A numeric vector with the increase in fraction of total height for every additional comparison to minimize overlap.
#' @param ref.group a character string specifying the reference group. If specified, for a given grouping variable, each of the group levels will be compared to the reference group (i.e. control group).
#' @param comparisons A list of length-2 vectors. The entries in the vector are either the names of 2 values on the x-axis or the 2 integers that correspond to the index of the columns of interest.(default:NULL)
#' @param ncol An interger. Set the col number in the facet plot.
#' @param select_metrics A series of character string. Select the alpha metrics to show. 
#' @param show A character string include pic (default), html.
#' @param palette A series of character string. Color palette.
#' @param html_width An interger. Set the html width.
#' @param html_height An interger. Set the html height.
#' @param mytheme Modify components of a theme according to the ggplot2::theme.
#' @param ... Additional parameters, see also \code{\link[ggsignif]{geom_signif}}
#' @rdname EMP_boxplot
#' @importFrom withr with_seed


EMP_boxplot.EMP_alpha_analysis <- function(obj,plot_category = 'default',method = 'wilcox.test',
                               estimate_group = NULL,group_level = 'default',
                               step_increase = 0.1,ref.group = NULL,comparisons = NULL,
                               ncol = NULL,select_metrics=NULL,show = 'pic',palette = NULL,
                               html_width=NULL,html_height=NULL,
                               mytheme = 'theme()',...) {
  call <- match.call()
  .get.plot_category.EMPT(obj) <- plot_category
  .get.history.EMPT(obj) <- call

  switch(.get.plot_category.EMPT(obj),
         "default" = {
           withr::with_seed(123,EMP_boxplot_alpha_default(EMPT=obj,method = method,
                               estimate_group = estimate_group,group_level = group_level,
                               step_increase = step_increase,ref.group = ref.group,comparisons = comparisons,
                               ncol = ncol,select_metrics=select_metrics,show = show,palette = palette,
                               html_width=html_width,html_height=html_height,
                               mytheme = mytheme,...))
         },
         "violin" = {
          withr::with_seed(123,EMP_boxplot_alpha_violin(EMPT=obj,method = method,
                               estimate_group = estimate_group,group_level = group_level,
                               ncol = ncol,select_metrics=select_metrics,show = show,palette = palette,
                               html_width=html_width,html_height=html_height,
                               mytheme = mytheme,...))
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



#' @import ggthemes
EMP_boxplot_alpha_default <- function (EMPT,method = 'wilcox.test',
                                       estimate_group = NULL,group_level = 'default',
                                       step_increase = 0.1,ref.group = NULL,comparisons = NULL,
                                       ncol = NULL,select_metrics = NULL,palette = NULL,
                                 show = 'pic',html_width=NULL,html_height=NULL,mytheme = 'theme()',...) {
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

  alpha_data <- .get.result.EMPT(EMPT,info='EMP_alpha_analysis') %>% suppressMessages() %>%
    dplyr::full_join(mapping,by= 'primary')

  ## clean the missing value in the group label
  if(any(is.na(alpha_data[[estimate_group]]))) {
    warning('Column ',estimate_group,' has beed deteced missing value, all related samples will be removed in the display!')
    alpha_data <- alpha_data %>% tidyr::drop_na(!!estimate_group)
  }

  # choose the compare group
  group_name <- unique(alpha_data[[estimate_group]])
  if (!is.null(ref.group)) {
    if (!any(ref.group %in% group_name)) {
      stop('Please check the parameter ref.group!')
    }
    group_combn <- combn(as.character(group_name),2) %>% as.data.frame() %>%
    dplyr::select(where(~ any(str_detect_multi(., ref.group))))
  }else{
    group_combn <- combn(as.character(group_name),2)
  }

  if (is.null(comparisons)) {
    comparisons <- list() 
    for (i in 1:ncol(group_combn)) {
      comparisons[[i]] <- group_combn[,i]
    }
    names(comparisons) <- 1:ncol(group_combn)    
  }else{
    comparisons <- comparisons
  } 

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
    ggsignif::geom_signif(comparisons = comparisons,test = method,step_increase = step_increase,...) +
    facet_wrap(ID~., scales = 'free', strip.position = 'top',ncol = ncol) +
    xlab(NULL) +
    ylab("Alpha Metrics") +
    ggtitle('Alpha analysis Plot') +
    scale_fill_manual(values = col_values) + 
    theme_bw() + 
    theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) + 
    eval(parse(text = paste0(mytheme)))


  alpha_plot[['html']] <- ggiraph::girafe(code = print(alpha_plot[['pic']]),width = html_width,height = html_height)

  .get.plot_deposit.EMPT(EMPT,info = 'alpha_analysis_plot') <- alpha_plot
  .get.plot_specific.EMPT(EMPT) <- show
  .get.estimate_group.EMPT(EMPT) <- estimate_group
  .get.algorithm.EMPT(EMPT) <- 'alpha_analysis_plot'
  .get.info.EMPT(EMPT) <- 'EMP_alpha_analysis_boxplot'
  class(EMPT) <- 'EMP_alpha_analysis_boxplot'
  EMPT
}


EMP_boxplot_alpha_violin <- function (EMPT,method = 'wilcox.test',
                                       estimate_group = NULL,group_level = 'default',
                                       step_increase = 0.1,ref.group = NULL,comparisons = NULL,
                                       ncol = NULL,select_metrics = NULL,palette = NULL,
                                 show = 'pic',html_width=NULL,html_height=NULL,mytheme = 'theme()',...) {
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

  alpha_data <- .get.result.EMPT(EMPT,info='EMP_alpha_analysis') %>% suppressMessages() %>%
    dplyr::full_join(mapping,by= 'primary')

  ## clean the missing value in the group label
  if(any(is.na(alpha_data[[estimate_group]]))) {
    warning('Column ',estimate_group,' has beed deteced missing value, all related samples will be removed in the display!')
    alpha_data <- alpha_data %>% tidyr::drop_na(!!estimate_group)
  }

  # choose the compare group
  group_name <- unique(alpha_data[[estimate_group]])
  if (!is.null(ref.group)) {
    if (!any(ref.group %in% group_name)) {
      stop('Please check the parameter ref.group!')
    }
    group_combn <- combn(as.character(group_name),2) %>% as.data.frame() %>%
    dplyr::select(where(~ any(str_detect_multi(., ref.group))))
  }else{
    group_combn <- combn(as.character(group_name),2)
  }

  if (is.null(comparisons)) {
    comparisons <- list() 
    for (i in 1:ncol(group_combn)) {
      comparisons[[i]] <- group_combn[,i]
    }
    names(comparisons) <- 1:ncol(group_combn)    
  }else{
    comparisons <- comparisons
  } 

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


  alpha_plot[['pic']] <-ggplot(alpha_data, aes(x = !!dplyr::sym(estimate_group), y = value, fill = !!dplyr::sym(estimate_group))) +
    geom_violin(position = position_dodge(width = 0.1), scale = 'width',alpha=0.8) +
    geom_boxplot(outlier.color=NA,fill="white", width=0.3) +
    ggiraph::geom_jitter_interactive(aes(tooltip = paste0(primary,' : ',value)),shape=21,position = position_jitter(height = .00000001))+
    ggsignif::geom_signif(comparisons = comparisons,test = method,step_increase = step_increase,...) +
    facet_wrap(ID~., scales = 'free', strip.position = 'top',ncol = ncol) +
    xlab(NULL) +
    ylab("Alpha Metrics") +
    ggtitle('Alpha analysis Plot') +
    scale_fill_manual(values = col_values) + 
    theme_bw() + 
    theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) + 
    eval(parse(text = paste0(mytheme)))


  alpha_plot[['html']] <- ggiraph::girafe(code = print(alpha_plot[['pic']]),width = html_width,height = html_height)

  .get.plot_deposit.EMPT(EMPT,info = 'alpha_analysis_plot') <- alpha_plot
  .get.plot_specific.EMPT(EMPT) <- show
  .get.estimate_group.EMPT(EMPT) <- estimate_group
  .get.algorithm.EMPT(EMPT) <- 'alpha_analysis_plot'
  .get.info.EMPT(EMPT) <- 'EMP_alpha_analysis_boxplot'
  class(EMPT) <- 'EMP_alpha_analysis_boxplot'
  EMPT
}





