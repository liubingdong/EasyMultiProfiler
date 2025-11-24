#' @import ggthemes
#' @import ggprism
#' @import ggplot2
#' @importFrom ggiraph geom_point_interactive
#' @importFrom ggpubr stat_pvalue_manual
#' @importFrom ggbreak scale_y_break
#' @importFrom ggiraph geom_jitter_interactive
#' @export
EMP_barplot <- function (EMPT,method = 'wilcox.test',
                               estimate_group = NULL,compare_group = NULL,
                               group_level = 'default',compare_group_level = 'default',
                               error_bar= 'default',label = 'p',
                               dot_size=3,bar_alpha=1,
                               step_increase=0.1,tip.length=0,ref.group=NULL,comparisons = NULL,
                               hide_point = FALSE,paired_group=NULL,paired_line=TRUE,p.adjust='fdr',
                               ncol = NULL,show = 'pic',
                               palette = NULL,test.arg=list(),
                               html_width=NULL,html_height=NULL,
                               mytheme = 'theme()',y_break = NULL,
                               facet = TRUE,...) {

  rlang::check_installed(c('ggbreak'), reason = 'for EMP_barplot().', action = install.packages) 

  #if (facet==FALSE & is.null(compare_group)) {
  #  stop("When compare_group is null, the facet can not be FALSE!")
  #}

  label <- match.arg(label, c("p", "p.adj","p.signif","p.adj.signif"))

  error_bar <- match.arg(error_bar, c("default", "upper"))
  error_bar <- switch(
    error_bar,
    default = "errorbar",
    upper = "uperrorbar"
  )

  primary <- value <- `.` <- NULL
  estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)

  if (is.null(palette)) {
    col_values <- .get.palette.EMPT(EMPT)
  }else {
    col_values <- palette
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
      data_plot[['pic']] <- ggplot(data, aes(x = !!dplyr::sym(estimate_group), y = value, fill = !!dplyr::sym(compare_group))) +
        # 添加柱子
        stat_summary(fun = mean,geom = 'col',alpha=bar_alpha,
                     color = 'black',
                     position = position_dodge(0.8),
                     width = 0.8) +
        # 修改颜色
        scale_fill_manual(values = col_values) +
        # 添加误差线
        stat_summary(fun.data = 'mean_se',geom = error_bar,
                     color = 'black',
                     position = position_dodge(width = 0.8),
                     width = 0.3,linewidth = 0.8) +
        geom_line(aes(group = !!dplyr::sym(paired_group)), color = 'gray', lwd = 0.5) +       
        # 添加抖动点
        geom_point_interactive(aes(color = factor(compare_group),
                   tooltip = paste0(primary,' : ',value),
                   group = !!group_expr),
                   show.legend = FALSE,# 不显示图例
                   position = position_jitterdodge(seed = 123,jitter.height = 0.000001,
                    jitter.width = 0, ### necessary
                    dodge.width = 0.75),
                   shape = 21,size = dot_size,fill = 'white') +
        # 抖动点边框颜色
        scale_color_manual(values = rep('black',3)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +    
        facet_wrap(ID~., scales = 'free', strip.position = 'top',ncol = ncol) +
        theme_bw() + 
        theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
        eval(parse(text = paste0(mytheme)))

      data_plot[['pic']] <- data_plot[['pic']] + stat_pvalue_manual(stat.test, inherit.aes = F,
                                  label = label,step_increase=step_increase,tip.length=tip.length,...)  
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

      data_plot[['pic']] <- ggplot(data, aes(x = !!dplyr::sym(estimate_group), y = value, fill = !!dplyr::sym(compare_group))) +
        # 添加柱子
        stat_summary(fun = mean,geom = 'col',alpha=bar_alpha,
                     color = 'black',
                     position = position_dodge(0.8),
                     width = 0.8) +
        # 修改颜色
        scale_fill_manual(values = col_values) +
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
                   position = position_jitterdodge(seed = 123,jitter.height = 0.000001,dodge.width = 0.75),
                   shape = 21,size = dot_size,fill = 'white') +
        # 抖动点边框颜色
        scale_color_manual(values = rep('black',3)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +     
        facet_wrap(ID~., scales = 'free', strip.position = 'top',ncol = ncol) +
        theme_bw() + 
        theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
        eval(parse(text = paste0(mytheme)))

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
                                  label = label,step_increase=step_increase,tip.length=tip.length,...)  

   }

  data_plot[['html']] <- ggiraph::girafe(ggobj = print(data_plot[['pic']]),width = html_width,height = html_height)

  .get.plot_deposit.EMPT(EMPT,info = 'EMP_assay_barplot') <- data_plot
  .get.plot_specific.EMPT(EMPT) <- show
  .get.estimate_group.EMPT(EMPT) <- estimate_group
  EMPT@algorithm <- 'EMP_assay_barplot'
  .get.info.EMPT(EMPT) <- 'EMP_assay_barplot'
  class(EMPT) <- 'EMP_assay_barplot'
  return(EMPT)
}

.show_EMP_assay_barplot <- function(obj,plot) {
  result <- .get.plot_deposit.EMPT(obj,info = 'EMP_assay_barplot')
  switch(plot,
         "pic" = print(result$pic),
         "html" = print(result$html)
  )
}



