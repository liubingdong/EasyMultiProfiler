#' Barplot for EMPT result
#'
#' @param obj EMPT object.
#' @param method statistics including t.test, wilcox.test, sign.test, emmeans.test, dunn.test, tukey.hsd.
#' @param estimate_group The primary grouping variable for estimation and comparison.Serves as the comparison group by default.
#' @param compare_group Optional subgroup analysis variable. When specified, the subgroup comparison performs.
#' @param group_level A string vector. Set the group order in the plot.
#' @param compare_group_level A string vector. Set the compare group order in the plot.
#' @param error_bar A character string including default and upper.
#' @param label A character string including p, p.adj, p.signif, p.adj.signif.
#' @param dot_size A numeric. Set the dot size.
#' @param dot_color A color character string.Dot color including white, black, ... or "group" for grouping-based colors.
#' @param bar_alpha A numeric. Set the bar alpha.
#' @param bar_color A color character to set the border of bar. Border color including black, ... or "group" for grouping-based colors.
#' @param step.increase A numeric vector with the increase in fraction of total height for every additional comparison to minimize overlap.
#' @param tip.length A numeric vector with the fraction of total height that the bar goes down to indicate the precise column
#' @param ref.group A character string specifying the reference group. If specified, for a given grouping variable, each of the group levels will be compared to the reference group (i.e. control group).
#' @param comparisons A list of length-2 vectors. The entries in the vector are either the names of 2 values on the x-axis or the 2 integers that correspond to the index of the columns of interest.(default:NULL)
#' @param hide_point A boolean. Whether show the dot or not.
#' @param paired_group A character string. Variable name corresponding to paired primary or sample.
#' @param paired_line A boolean. Whether show the paired line when paired test activated.(defalut:TRUE)
#' @param p.adjust A character string. Adjust P-values for Multiple Comparisons inluding fdr, holm, hochberg, hommel, bonferroni, BH, BY. (default:fdr)
#' @param ncol An interger. Set the col number in the facet plot.
#' @param show A character string include pic (default), html.
#' @param seed An interger. Set the random seed for the point layout.
#' @param palette A series of character string. Color palette.
#' @param test.arg Additional arguments for the test method.
#' @param html_width html_width
#' @param html_height html_height
#' @param mytheme Modify components of a theme according to the \code{\link[ggplot2]{theme}} and \code{\link[ggplot2]{ggtheme}}.
#' @param y_break Break point.
#' @param facet If TRUE (default), each variable is displayed in a separate panel. If FALSE, all variables are overlaid in a single panel.
#' @param ... additional parameters, see also \code{\link[ggpubr]{stat_pvalue_manual}}
#' @rdname EMP_barplot
#' @import ggplot2
#' @importFrom ggiraph geom_point_interactive
#' @importFrom ggpubr stat_pvalue_manual
#' @importFrom ggbreak scale_y_break
#' @importFrom ggiraph geom_jitter_interactive
#' @export
#' @examples
#' data(MAE)
#' # plot
#' MAE |>
#'   EMP_assay_extract('host_gene',pattern = 'A1BG') |> 
#'   EMP_barplot(estimate_group = 'Group')
#' 
#' # set error bar
#' MAE |>
#'   EMP_assay_extract('host_gene',pattern = 'A1BG') |> 
#'   EMP_barplot(estimate_group = 'Group',error_bar = 'upper')
#' 
#' # for subgroup
#' MAE |>
#'   EMP_assay_extract('host_gene',pattern = 'A1BG') |> 
#'   EMP_barplot(estimate_group = 'Group',compare_group='Sex')
#' 
#' # set signif
#' MAE |>
#'   EMP_assay_extract('host_gene',pattern = 'A1BG') |> 
#'   EMP_barplot(estimate_group = 'Group',compare_group='Sex',label = 'p.signif')
#' 
#' # combine the plot
#' MAE |>
#'   EMP_assay_extract('host_gene',pattern = 'A1BG') |> 
#'   EMP_barplot(estimate_group = 'Group',compare_group='Sex',facet = FALSE)
#' 
#' # set the comparision
#' MAE |>
#'   EMP_assay_extract('host_gene',pattern = 'A1BG') |> 
#'   EMP_barplot(estimate_group = 'Group',compare_group='sub_group',
#'               ref.group = 'A',
#'               facet = FALSE)
#' MAE |>
#'   EMP_assay_extract('host_gene',pattern = 'A1BG') |> 
#'   EMP_barplot(estimate_group = 'Group',compare_group='sub_group',
#'               comparisons = list(c('A','B')),
#'               facet = FALSE)
#' 
#' # for paired test
#' MAE |>
#'   EMP_assay_extract('host_gene',pattern = 'A1BG') |> 
#'   EMP_barplot(estimate_group = 'sub_group',
#'               paired_group = 'patient')
#' 
#' # plot as prism style
#' if (!requireNamespace("ggprism", quietly=TRUE)) install.packages("ggprism")
#' library(ggprism)
#' MAE |>
#'   EMP_assay_extract('host_gene',pattern = 'A1BG') |> 
#'   EMP_barplot(estimate_group = 'Group',compare_group='Status',
#'               palette = c('white','red','black'),error_bar = 'upper',
#'               method = 'tukey.hsd',facet = FALSE,
#'               mytheme ='theme_prism()')
#' MAE |>
#'   EMP_assay_extract('host_gene',pattern = 'A1BG') |> 
#'   EMP_barplot(estimate_group = 'Group',compare_group='Status',
#'               bar_alpha = 0,bar_color = 'group',
#'               palette = c('yellow','red','green'),error_bar = 'upper',
#'               method = 'tukey.hsd',dot_color = 'group',
#'               mytheme ='theme_prism()')
EMP_barplot <- function (obj,method = 'wilcox.test',
                           estimate_group = NULL,compare_group = NULL,
                           group_level = 'default',compare_group_level = 'default',
                           error_bar= 'default',label = 'p',
                           dot_size=3,dot_color='white',bar_alpha=1,bar_color='black',
                           step.increase=0.1,tip.length=0,ref.group=NULL,comparisons = NULL,
                           hide_point = FALSE,paired_group=NULL,paired_line=TRUE,p.adjust='fdr',
                           ncol = NULL,show = 'pic',seed=123,
                           palette = NULL,test.arg=list(),
                           html_width=NULL,html_height=NULL,
                           mytheme = 'theme()',y_break = NULL,
                           facet = TRUE,...) {

  rlang::check_installed(c('ggbreak'), reason = 'for EMP_barplot().', action = install.packages) 

  if (!is(obj,'EMPT')) {
    stop("Please check the input data in EMPT format!")
  }else{
    EMPT <- obj
  }

  label <- match.arg(label, c("p", "p.adj","p.signif","p.adj.signif"))

  error_bar <- match.arg(error_bar, c("default", "upper"))
  error_bar <- switch(
    error_bar,
    default = "errorbar",
    upper = "uperrorbar"
  )
  
  if (dot_color == 'group') {
    dot_color <- NULL
  }

  primary <- value <- `.` <- NULL
  estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)

  if (is.null(palette)) {
    col_values <- .get.palette.EMPT(EMPT)
  }else {
    col_values <- palette
  }

  if (length(bar_color) == 1 && bar_color == "group") {
    bar_color <- col_values
  } else if (length(bar_color) == 1) {
    bar_color <- rep(bar_color, length(col_values))
  }else{
    bar_color <-bar_color
  }

  EMPT %<>% .group_level_modified(estimate_group = estimate_group,
                                  group_level = group_level)
  
  if (!is.null(compare_group)) {
    EMPT %<>% .group_level_modified(estimate_group = compare_group,
                                  group_level = compare_group_level)
  }


  mapping <- .get.mapping.EMPT(EMPT) %>% dplyr::select(primary,dplyr::any_of(c(!!estimate_group,!!compare_group,!!paired_group)))
  
  data_info <- .get.info.EMPT(EMPT)
  data <-.get.result.EMPT(EMPT,info = data_info) %>% suppressMessages() %>% dplyr::left_join(mapping,by ='primary') 

  ## clean the missing value in the group label
  if(any(is.na(data[[estimate_group]]))) {
    warning('Column ',estimate_group,' has beed deteced missing value, all related samples will be removed in the display!')
    data <- data %>% tidyr::drop_na(!!estimate_group)
  }

  if (!is.null(compare_group)) {
    if(any(is.na(data[[compare_group]]))) {
      warning('Column ',estimate_group,' has beed deteced missing value, all related samples will be removed in the display!')
      data <- data %>% tidyr::drop_na(!!compare_group) 
    }   
  }


  if (!is.null(compare_group)) {
    compare_group <- compare_group
  }else{
    compare_group <- estimate_group
  }

  data %<>%  tidyr::pivot_longer(cols = c(-primary,-dplyr::any_of(c(!!estimate_group,!!compare_group,!!paired_group))),
                                       names_to = 'ID',
                                       values_to = 'value')
  n_ID <- data[['ID']] |> unique() |> length()

  if (facet == TRUE) {
    stat.test <- do.call(
      "stat_test",
      c(
        list(data = data,
             estimate_group = estimate_group,
             compare_group = compare_group,
             method = method,
             value = 'value',
             comparisons=comparisons,
             ref.group=ref.group,
             paired_group=paired_group,
             p.adjust=p.adjust,
             facet.by='ID',silent=TRUE),  # base parameter
        test.arg  # addtional parameter
      )
    )
  }else{
    if (n_ID > 1) {
      if (is.null(compare_group)) {
        compare_group <- estimate_group
      }
      estimate_group <- 'ID'
    }

    stat.test <- do.call(
      "stat_test",
      c(
        list(data = data,
             estimate_group = estimate_group,
             compare_group = compare_group,
             method = method,
             value = 'value',
             comparisons=comparisons,
             ref.group=ref.group,
             paired_group=paired_group,
             p.adjust=p.adjust,
             facet.by=NULL,silent=TRUE),  # base parameter
        test.arg  # addtional parameter
      )
    )
  }

  if (is.null(compare_group)) {
    group_expr <- sym(estimate_group)
  } else {
    group_expr <- expr(interaction(!!sym(estimate_group), !!sym(compare_group)))
  }

  data_plot <- list()
  if (!is.null(paired_group)) {
      data_plot[['pic']] <- (ggplot(data, aes(x = !!dplyr::sym(estimate_group), y = value, 
                  fill = !!dplyr::sym(compare_group),
                  color = !!dplyr::sym(compare_group))) +
        # 添加柱子
        stat_summary(fun = mean,geom = 'col',alpha=bar_alpha,
                     color = 'black',
                     position = position_dodge(width=0.8),
                     width = 0.8) +
        # 修改颜色
        scale_fill_manual(values = col_values) +
        scale_color_manual(values = bar_color) +
        # 添加误差线
        stat_summary(fun.data = 'mean_se',geom = error_bar,
                     #color = 'black',
                     position = position_dodge(width = 0.8),
                     width = 0.3,linewidth = 0.8) +
        geom_line(aes(group = !!dplyr::sym(paired_group)), color = 'gray', lwd = 0.5) +       
        # 添加抖动点
        geom_point_interactive(aes(color = factor(compare_group),
                   tooltip = paste0(primary,' : ',value),
                   group = !!group_expr),
                   show.legend = FALSE,# 不显示图例
                   position = position_jitterdodge(seed = seed,jitter.height = 0.000001,
                    jitter.width = 0, ### necessary
                    dodge.width = 0.8),
                   shape = 21,size = dot_size,fill = dot_color,color = 'black') +
        # 抖动点边框颜色
        #scale_color_manual(values = rep('black',3)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +    
        facet_wrap(ID~., scales = 'free', strip.position = 'top',ncol = ncol) +
        theme_bw() + 
        theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
        eval(parse(text = paste0(mytheme))) ) |> suppressWarnings()

      data_plot[['pic']] <- data_plot[['pic']] + stat_pvalue_manual(stat.test, inherit.aes = F,
                                  label = label,step.increase=step.increase,tip.length=tip.length,...)  
      if (facet == FALSE & n_ID > 1) {
        data_plot[['pic']] <- data_plot[['pic']] + facet_null() + xlab(label = NULL)
      }  

      if (!is.null(y_break)) {
        data_plot[['pic']] <- data_plot[['pic']] + ggbreak::scale_y_break(breaks = y_break)
      } 

      # Hide the line      
      if (paired_line == FALSE) {
        data_plot[['pic']][['layers']] <- data_plot[['pic']][['layers']][-3]
      }  

      # Hide the point      
      if (hide_point == TRUE) {
        data_plot[['pic']][['layers']] <- data_plot[['pic']][['layers']][-4]
      }  

   }else{

      data_plot[['pic']] <- (ggplot(data, aes(x = !!dplyr::sym(estimate_group), y = value, 
                  fill = !!dplyr::sym(compare_group),
                  color = !!dplyr::sym(compare_group))) +
        # 添加柱子
        stat_summary(fun = mean,geom = 'col',alpha=bar_alpha,
                     #color = 'black',
                     position = position_dodge(0.8),
                     width = 0.8) +
        # 修改颜色
        scale_fill_manual(values = col_values) +
        scale_color_manual(values = bar_color) +
        # 添加误差线
        stat_summary(fun.data = 'mean_se',geom = error_bar,
                     color = 'black',
                     position = position_dodge(width = 0.8),
                     width = 0.3,linewidth = 0.8) +
        # 添加抖动点
        geom_point_interactive(aes(color = factor(compare_group),
                   tooltip = paste0(primary,' : ',value),
                   group = !!group_expr),
                   show.legend = FALSE,# 不显示图例
                   position = position_jitterdodge(seed = seed,jitter.height = 0.000001,dodge.width = 0.8),
                   shape = 21,size = dot_size,fill = dot_color,color = 'black') +
        # 抖动点边框颜色
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +     
        facet_wrap(ID~., scales = 'free', strip.position = 'top',ncol = ncol) +
        theme_bw() + 
        theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
        eval(parse(text = paste0(mytheme))) ) |> suppressWarnings()

      if (facet == FALSE & n_ID > 1) {
        data_plot[['pic']] <- data_plot[['pic']] + facet_null() + xlab(label = NULL)
      }  

      # Hide the point      
      if (hide_point == TRUE) {
        data_plot[['pic']][['layers']] <- data_plot[['pic']][['layers']][-3]
      }  

      if (!is.null(y_break)) {
        data_plot[['pic']] <- data_plot[['pic']] + ggbreak::scale_y_break(breaks = y_break)
      }  

      data_plot[['pic']] <- data_plot[['pic']] + stat_pvalue_manual(stat.test, inherit.aes = F,
                                  label = label,step.increase=step.increase,tip.length=tip.length,...)  

   }

  data_plot[['html']] <- ggiraph::girafe(ggobj = print(data_plot[['pic']]),width = html_width,height = html_height)

  .get.plot_deposit.EMPT(EMPT,info = 'EMP_barplot') <- data_plot
  .get.plot_specific.EMPT(EMPT) <- show
  .get.estimate_group.EMPT(EMPT) <- estimate_group
  EMPT@algorithm <- 'EMP_barplot'
  .get.info.EMPT(EMPT) <- 'EMP_barplot'
  class(EMPT) <- 'EMP_barplot'
  return(EMPT)
}

.show_EMP_barplot <- function(obj,plot) {
  result <- .get.plot_deposit.EMPT(obj,info = 'EMP_barplot')
  switch(plot,
         "pic" = print(result$pic),
         "html" = print(result$html)
  )
}



