#' EMP volcanol plot
#'
#' @param EMPT EMPT object
#' @param y A character string. Select the pvalue from the EMP_diff_analysis.
#' @param palette A series of character string. Color palette.
#' @param key_feature A series of character string. Label the some feature.
#' @param show A character string include pic (default), html.
#' @param html_width An interger. Set the html width.
#' @param html_height An interger. Set the html height.
#' @param mytheme Modify components of a theme according to the ggplot2::theme.
#' @param ... Further parameters passed to ggrepel::geom_text_repel
#' @importFrom ggrepel geom_text_repel
#' @return EMPT object
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

EMP_volcanol_plot <- function(EMPT,y='pvalue',palette = NULL,show = 'pic',
                           html_width=NULL,html_height=NULL,key_feature=NULL,
                           mytheme = 'theme()',...) {
  sign_group <- log2FC <- feature <- color <- key <- NULL
  if (is.null(palette)) {
    col_values <- .get.palette.EMPT(EMPT)
  }else {
    col_values = palette
  }

  assay_name <- .get.assay_name.EMPT(EMPT)
  if (!assay_name %in% c('counts','relative','integer')) {
    stop('Only counts, relative and integer data in two groups provide volcanol plot!')
  }

  method <- .get.method.EMPT(EMPT)
  estimate_group_info <- .get.estimate_group_info.EMPT(EMPT)

  title_info <- paste0('Comparison ',estimate_group_info,' in ',method)


  data <- .get.result.EMPT(EMPT,info='EMP_diff_analysis') %>% suppressMessages()
  # check the group num
  group_num <-  data %>% dplyr::pull(sign_group) %>% na.omit() %>% unique() %>% length()
  if (group_num != 2) {
    stop('volcanol plot only support two groups comparition!')
  }

  # check y
  y_support <- c("pvalue","holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                 "fdr", "none",'padj')
  if (!y %in% y_support) {
    stop('y axis only support pvalue,holm,hochberg,hommel,bonferroni,BH,BY,fdr,none,padj!')
  }

  data %<>%
    dplyr::mutate(color = ifelse(log2FC > 0 & -log10(!!dplyr::sym(y)) > 1.3,
                                 yes = "UP",
                                 no = ifelse(log2FC < 0 & -log10(!!dplyr::sym(y)) > 1.3,
                                             yes = "DOWN",
                                             no = "none")))
  # confirm x y axis limits
  ylim_max <- -log10(data[[y]]) %>% max(na.rm = T)
  ylim_break <- c(0,1.3,5,10,20,30,40,50,60,100)
  idy <- (ylim_break < ylim_max) %>% sum() +1
  ylim_break <- ylim_break[1:idy]


  xlim_max <-data$log2FC %>% abs() %>% max(na.rm = T)
  xlim_break <-c(0,4,6,8,10,12,24,30,50,100)
  idx <- (xlim_break < xlim_max) %>% sum() +1
  xlim_break <- xlim_break[1:idx]
  xlim_break <- c(xlim_break,xlim_break*-1) %>% unique() %>% sort()

  # set the key feature
  data <- data %>% dplyr::mutate(key = ifelse(feature %in% key_feature,'key',NA))
  if (is.null(key_feature)) {
      p1 <- ggplot(data, aes(x = log2FC, y = -log10(!!dplyr::sym(y)))) +
          ggiraph::geom_jitter_interactive(aes(tooltip = paste0(feature,' : ',log2FC,' ',!!dplyr::sym(y)),
                                               color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T)+
          # add gene points
          ggtitle(label = title_info) +  # add title
          xlab(expression(log[2]("fold change"))) + # x-axis label
          ylab(substitute(-log[10](x), list(x = as.name(y)))) + # y-axis label
          geom_vline(xintercept = 0, colour = "black") + # add line at 0
          geom_hline(yintercept = 1.3, colour = "black") + # p(0.05) = 1.3
          scale_x_continuous(breaks=xlim_break,limits = c(min(xlim_break),max(xlim_break))) +  # set x axis
          # 根据需要调节纵坐标
          scale_y_continuous(breaks = ylim_break,limits = c(min(ylim_break),max(ylim_break)),trans = "log1p") + # set y axis

          scale_color_manual(values = c("UP" = col_values[1],
                                        "DOWN" = col_values[2],
                                        "none" = "#636363")) +
          theme_bw() + # clean up theme
          theme(legend.position = "none")
  }else{
      p1 <- ggplot(data, aes(x = log2FC, y = -log10(!!dplyr::sym(y)),label = feature)) +
            ggrepel::geom_text_repel(data=subset(data,key == 'key'),...) +
            ggiraph::geom_jitter_interactive(aes(tooltip = paste0(feature,' : ',log2FC,' ',!!dplyr::sym(y)),
                                                 color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T)+
            # add gene points
            ggtitle(label = title_info) +  # add title
            xlab(expression(log[2]("fold change"))) + # x-axis label
            ylab(substitute(-log[10](x), list(x = as.name(y)))) + # y-axis label
            geom_vline(xintercept = 0, colour = "black") + # add line at 0
            geom_hline(yintercept = 1.3, colour = "black") + # p(0.05) = 1.3
            scale_x_continuous(breaks=xlim_break,limits = c(min(xlim_break),max(xlim_break))) +  # set x axis
            # 根据需要调节纵坐标
            scale_y_continuous(breaks = ylim_break,limits = c(min(ylim_break),max(ylim_break)),trans = "log1p") + # set y axis

            scale_color_manual(values = c("UP" = col_values[1],
                                          "DOWN" = col_values[2],
                                          "none" = "#636363")) +
            theme_bw() + # clean up theme
            theme(legend.position = "none")
  }
 

  data_plot <- list()
  data_plot[['pic']] <- p1
  data_plot[['html']] <- ggiraph::girafe(code = print(p1),width = html_width,height = html_height)

  .get.plot_deposit.EMPT(EMPT,info = 'EMP_diff_volcanol_plot') <- data_plot
  .get.plot_specific.EMPT(EMPT) <- show
  EMPT@algorithm <- 'EMP_diff_volcanol_plot'
  .get.info.EMPT(EMPT) <- 'EMP_diff_volcanol_plot'
  class(EMPT) <- 'EMP_diff_volcanol_plot'
  EMPT
}


.show_EMP_diff_volcanol_plot <- function(obj,plot) {
  result <- .get.plot_deposit.EMPT(obj,info = 'EMP_diff_volcanol_plot')
  switch(plot,
         "pic" = print(result$pic),
         "html" = print(result$html)
  )
}
