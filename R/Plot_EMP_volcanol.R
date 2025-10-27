#' @importFrom ggrepel geom_text_repel
EMP_volcanol_plot_default <- function(EMPT,y='pvalue',palette = NULL,show = 'pic',
                           html_width=NULL,html_height=NULL,key_feature=NULL,threshold_x = 0,dot_size = 1.75,dot_alpha= 1,
                           mytheme = 'theme()',...) {
  sign_group <- log2FC <- feature <- color <- key <- NULL
  if (is.null(palette)) {
    col_values <- .get.palette.EMPT(EMPT)
    col_values <- c(col_values[1:2],"#636363")
  }else {
    col_values = palette
    if (length(col_values) <= 1) {
      stop("Parameter palette needs more than two color at least!")
    }else if (length(col_values) == 2){
      col_values[3] <- "#636363"
    }else{
      col_values <- col_values
    }
  }

  #assay_name <- .get.assay_name.EMPT(EMPT)
  #if (!assay_name %in% c('counts','relative','integer')) {
  #  stop('Only counts, relative and integer data in two groups provide volcanol plot!')
  #}

  method <- .get.method.EMPT(EMPT)
  estimate_group_info <- .get.estimate_group_info.EMPT(EMPT)

  title_info <- paste0('Comparison ',estimate_group_info,' in ',method)


  data <- .get.result.EMPT(EMPT,info='EMP_diff_analysis') %>% suppressMessages()
  # check the group num
  group_num <-  data %>% dplyr::pull(sign_group) %>% na.omit() %>% unique() %>% length()
  if (group_num > 2) {
    stop('volcanol plot only support no more than 2 groups comparition!')
  }

  # check y
  y_support <- c("pvalue","holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                 "fdr", "none",'padj')
  if (!y %in% y_support) {
    stop('y axis only support pvalue,holm,hochberg,hommel,bonferroni,BH,BY,fdr,none,padj!')
  }

  data %<>%
    dplyr::mutate(color = ifelse(log2FC > threshold_x & -log10(!!dplyr::sym(y)) > 1.3,
                                 yes = "UP",
                                 no = ifelse(log2FC < -threshold_x & -log10(!!dplyr::sym(y)) > 1.3,
                                             yes = "DOWN",
                                             no = "none")))
  # confirm x y axis limits
  ylim_max <- -log10(data[[y]]) %>% max(na.rm = T)
  ylim_break <- c(0,1.3,5,10,20,30,40,50,60,100)
  idy <- (ylim_break < ylim_max) %>% sum() +1
  ylim_break <- ylim_break[1:idy]


  xlim_max <-data$log2FC %>% abs() %>% max(na.rm = T)
  xlim_break <-c(0,1,2,4,6,8,10,12,24,30,50,100,threshold_x) |> unique() |> sort()
  idx <- (xlim_break < xlim_max) %>% sum() +1
  xlim_break <- xlim_break[1:idx]
  xlim_break <- c(xlim_break,xlim_break*-1) %>% unique() %>% sort()

  # set the key feature
  data <- data %>% dplyr::mutate(key = ifelse(feature %in% key_feature,'key',NA))

  p <- ggplot(data, aes(x = log2FC, y = -log10(!!dplyr::sym(y)),label = feature)) +
      ggiraph::geom_jitter_interactive(aes(tooltip = paste0(feature,' : ',log2FC,' ',!!dplyr::sym(y)),
                                           color = color), size = dot_size, alpha = dot_alpha, na.rm = T)+
      # add gene points
      ggtitle(label = title_info) +  # add title
      xlab(expression(log[2]("fold change"))) + # x-axis label
      ylab(substitute(-log[10](x), list(x = as.name(y)))) + # y-axis label
      geom_hline(yintercept = 1.3, colour = "black",linetype="twodash",linewidth=0.5) + # p(0.05) = 1.3
      scale_x_continuous(breaks=xlim_break,limits = c(min(xlim_break),max(xlim_break))) +  # set x axis
      # 根据需要调节纵坐标
      scale_y_continuous(breaks = ylim_break,limits = c(min(ylim_break),max(ylim_break)),trans = "log1p") + # set y axis

      scale_color_manual(values = c("UP" = col_values[1],
                                    "DOWN" = col_values[2],
                                    "none" = col_values[3])) 
  if (!is.null(key_feature)) {
    p <- p + ggrepel::geom_text_repel(data=subset(data,key == 'key'),...)
  }

  if (threshold_x == 0) {
    p <- p + geom_vline(xintercept = 0, colour = "black") +
      theme_bw() + # clean up theme
      theme(legend.position = "none")+
      eval(parse(text = paste0(mytheme)))
  }else{
    p <- p + 
      geom_vline(aes(xintercept= threshold_x), colour="black", linetype="twodash",linewidth=0.5) +
      geom_vline(aes(xintercept= -threshold_x), colour="black", linetype="twodash",linewidth=0.5) +
      theme_bw() + # clean up theme
      theme(legend.position = "none")+
      eval(parse(text = paste0(mytheme)))
  }

  data_plot <- list()
  data_plot[['pic']] <- p
  data_plot[['html']] <- ggiraph::girafe(ggobj = print(p),width = html_width,height = html_height)

  .get.plot_deposit.EMPT(EMPT,info = 'EMP_diff_volcanol_plot') <- data_plot
  .get.plot_specific.EMPT(EMPT) <- show
  EMPT@algorithm <- 'EMP_diff_volcanol_plot'
  .get.info.EMPT(EMPT) <- 'EMP_diff_volcanol_plot'
  class(EMPT) <- 'EMP_diff_volcanol_plot'
  EMPT
}

#' EMP volcanol plot
#'
#' @param obj EMPT object
#' @param plot_category An interger.More plot style.(under constrution)
#' @param y A character string. Select the pvalue from the EMP_diff_analysis.
#' @param palette A series of character string. Color palette.
#' @param threshold_x Set the threshold for log2FC.
#' @param dot_size A number. Set the dot size.
#' @param dot_alpha A number. Set the dot alpha.
#' @param key_feature A series of character string. Label the some feature.
#' @param show A character string include pic (default), html.
#' @param html_width An interger. Set the html width.
#' @param html_height An interger. Set the html height.
#' @param mytheme Modify components of a theme according to the \code{\link[ggplot2]{theme}} and \code{\link[ggplot2]{ggtheme}}.
#' @param ... Further parameters passed to \code{\link[ggrepel]{geom_text_repel}}.
#' @export
#'
#' @examples
#' data(MAE)
#' MAE |>
#'   EMP_assay_extract(experiment = 'geno_ec') |>
#'   EMP_diff_analysis(method='DESeq2',.formula = ~Group)  |>
#'   EMP_volcanol_plot(show='html')
#'
#' # label feature
#' MAE |>
#'   EMP_decostand(experiment = 'geno_ec',method = 'integer') |>
#'   EMP_diff_analysis(method='DESeq2',.formula = ~Group)  |>
#'   EMP_volcanol_plot(key_feature = c('3.6.1.62','1.5.3.19')) 
#'
#' # Addtionl parameters will pass into ggrepel::geom_text_repel.
#' MAE |>
#'   EMP_decostand(experiment = 'geno_ec',method = 'integer') |>
#'   EMP_diff_analysis(method='DESeq2',.formula = ~Group)  |>
#'   EMP_volcanol_plot(key_feature = c('3.6.1.62','1.5.3.19'),color = "white",
#'                     bg.color = "grey30",bg.r = 0.15)
#' MAE |>
#'   EMP_decostand(experiment = 'geno_ec',method = 'integer') |>
#'   EMP_diff_analysis(method='DESeq2',.formula = ~Group)  |>
#'   EMP_volcanol_plot(key_feature = c('3.6.1.62','1.5.3.19'),
#'                     min.segment.length = 0, seed = 42, box.padding = 0.5) ## Add arrow
#'
#' # More parameter
#' MAE |>
#'   EMP_assay_extract(experiment = 'geno_ec') |>
#'   EMP_diff_analysis(method='DESeq2',.formula = ~Group)  |>
#'   EMP_volcanol_plot(show='html',key_feature = c('3.6.1.64','1.5.3.19'),
#'                     palette = c('#FA7F6F','#96C47D','#BEB8DC'),
#'                     dot_size = 3,threshold_x = 0.5,mytheme = "theme_light()",
#'                     min.segment.length = 0, seed = 42, box.padding = 0.5)

EMP_volcanol_plot <- function(obj,plot_category=1,y='pvalue',palette = NULL,show = 'pic',
                           html_width=NULL,html_height=NULL,key_feature=NULL,threshold_x = 0,dot_size = 1.75,dot_alpha= 1,
                           mytheme = 'theme()',...) {
  call <- match.call()
  .get.plot_category.EMPT(obj) <- plot_category
  .get.history.EMPT(obj) <- call

  switch(.get.plot_category.EMPT(obj),
         "1" = {
           EMP_volcanol_plot_default(EMPT=obj,y=y,palette = palette,show = show,
                           html_width=html_width,html_height=html_height,key_feature=key_feature,
                           threshold_x = threshold_x,dot_size = dot_size,dot_alpha= dot_alpha,
                           mytheme = mytheme,...)
         },
         "2" = {
          message("Under construction!")
         }

  )

}

.show_EMP_diff_volcanol_plot <- function(obj,plot) {
  result <- .get.plot_deposit.EMPT(obj,info = 'EMP_diff_volcanol_plot')
  switch(plot,
         "pic" = print(result$pic),
         "html" = print(result$html)
  )
}
