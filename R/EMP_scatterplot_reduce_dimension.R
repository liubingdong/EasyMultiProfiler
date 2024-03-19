.show_EMP_dimension_analysis_scatterplot <- function(obj,plot) {
  result <- .get.plot_deposit.EMPT(obj,info = 'dimension_analysis_scatterplot')
  switch(plot,
         "p12" = print(result$pic$p12),
         "p23" = print(result$pic$p23),
         "p13" = print(result$pic$p13),
         "p12html" = print(result$html$p12_html),
         "p23html" = print(result$html$p23_html),
         "p13html" = print(result$html$p13_html)
  )
}



#' Title
#'
#' @param EMPT wait_for_add
#' @param seed wait_for_add
#' @param group_level wait_for_add
#' @param show wait_for_add
#' @param distance_for_adonis wait_for_add
#' @param estimate_group wait_for_add
#' @param palette wait_for_add
#' @param method wait_for_add
#' @param key_samples wait_for_add
#' @param ellipse wait_for_add
#' @param width wait_for_add
#' @param height wait_for_add
#' @param force_adonis force_adonis
#' @importFrom ggpubr stat_compare_means
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 geom_boxplot
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_scatterplot_reduce_dimension  <- function(EMPT,seed=123,group_level='default',
                                           show='p12',distance_for_adonis=NULL,force_adonis=FALSE,
                                           estimate_group=NULL,palette=NULL,
                                           method='t.test',key_samples = NULL,
                                           ellipse = NULL,width=15,height=15){
  primary <- Group <- NULL
  call <- match.call()
  deposit <- list()

  palette <- .get.palette.EMPT(EMPT)

  estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)

  if (is.null(palette)) {
    col_values <- .get.palette.EMPT(EMPT)
  }else {
    col_values = palette
  }

  EMPT %<>% .group_level_modified(estimate_group = estimate_group,group_level = group_level)
  mapping <- .get.mapping.EMPT(EMPT) %>% dplyr::select(primary,!!estimate_group)

  colnames(mapping) <- c("primary","Group")
  length=length(unique(as.character(mapping$primary)))
  times1=length%/%8
  res1=length%%8
  times2=length%/%5
  res2=length%%5
  col1=rep(1:8,times1)
  col=c(col1,1:res1)
  pich1=rep(c(21:24),times2)
  pich=c(pich1,15:(15+res2))
  ###########bray
  ###PCoA分析--轴12
  data <- .get.result.EMPT(EMPT,info='EMP_dimension_analysis')[['dimension_coordinate']] %>% suppressMessages()
  plotdata <- dplyr::left_join(data,mapping,by='primary')
  pc1 <-.get.result.EMPT(EMPT,info='EMP_dimension_analysis')[['dimension_axis']][1] %>% suppressMessages()
  pc2 <-.get.result.EMPT(EMPT,info='EMP_dimension_analysis')[['dimension_axis']][2] %>% suppressMessages()
  try(pc3 <-.get.result.EMPT(EMPT,info='EMP_dimension_analysis')[['dimension_axis']][3] %>% suppressMessages(),silent = T)
  #plotdata$Group <- factor(plotdata$Group,levels = name_group)



  #PC1和PC2的显著性检验(PC1,PC2,PC3进行组间差异检验)
  if(!method %in% c('t.test','wilcox.test')){
    stop("method should be t.test or wilcox.test! ")
  }
  group_combn=combn(as.character(unique(plotdata$Group)),2)
  compare=plyr::alply(group_combn,2)


  axis_name <- colnames(data)[-1]
  axis_num <- length(unique(axis_name))


  #相须图绘制
  p1 <- ggplot(plotdata,aes(Group,!!dplyr::sym(axis_name[1]))) +
    geom_boxplot(aes(fill = Group),outlier.colour = NA) +scale_fill_manual(values=col_values)+
    ggpubr::stat_compare_means(comparisons = compare,method= method)+
    coord_flip() +ggiraph::geom_point_interactive(aes(tooltip = paste0(primary,' : ',round(!!dplyr::sym(axis_name[1]),2))),position = "jitter")+
    theme_bw()+
    theme(axis.ticks.length = unit(0.4,"lines"),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_text(colour='black',size=20,face = "bold"),
          axis.text.x=element_blank(),
          legend.position = "none")

  p2 <- ggplot(plotdata,aes(Group,!!dplyr::sym(axis_name[2]))) +
    geom_boxplot(aes(fill = Group),outlier.colour = NA) +scale_fill_manual(values=col_values)+
    ggpubr::stat_compare_means(comparisons = compare,method=method)+ggiraph::geom_point_interactive(aes(tooltip = paste0(primary,' : ',round(!!dplyr::sym(axis_name[2]),2))),position = "jitter")+
    theme_bw()+
    theme(axis.ticks.length = unit(0.4,"lines"),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_text(colour='black',size=20,angle = 45,
                                   vjust = 1,hjust = 1,face = "bold"),
          axis.text.y=element_blank(),
          legend.position = "none")

  p2_r <- ggplot(plotdata,aes(Group,!!dplyr::sym(axis_name[2]))) +
    geom_boxplot(aes(fill = Group),outlier.colour = NA) +scale_fill_manual(values=col_values)+
    ggpubr::stat_compare_means(comparisons = compare,method=method)+coord_flip() +ggiraph::geom_point_interactive(aes(tooltip = paste0(primary,' : ',round(!!dplyr::sym(axis_name[2]),2))),position = "jitter")+
    theme_bw()+
    theme(axis.ticks.length = unit(0.4,"lines"),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_text(colour='black',size=20,face = "bold"),
          axis.text.x=element_blank(),
          legend.position = "none")


  if (axis_num == 3) {
    p3 <- ggplot(plotdata,aes(Group,!!dplyr::sym(axis_name[3]))) + scale_fill_manual(values=col_values) +
      geom_boxplot(aes(fill = Group),outlier.colour = NA) +
      ggpubr::stat_compare_means(comparisons = compare,method=method) +ggiraph::geom_point_interactive(aes(tooltip = paste0(primary,' : ',round(!!dplyr::sym(axis_name[3]),2))),position = "jitter")+
      theme_bw()+
      theme(axis.ticks.length = unit(0.4,"lines"),
            axis.ticks = element_line(color='black'),
            axis.line = element_line(colour = "black"),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_text(colour='black',size=20,angle = 45,
                                     vjust = 1,hjust = 1,face = "bold"),
            axis.text.y=element_blank(),
            legend.position = "none")
  }



  if (!is.null(key_samples)) {
    group_info <- plotdata$Group %>% unique()
    plotdata  %<>% dplyr::mutate(Group = dplyr::if_else(primary %in% key_samples, "Key", Group))
    plotdata$Group <- factor(plotdata$Group,levels = c(group_info,'Key'))
  }

  #PCoA结果图绘制
  p12<-ggplot(plotdata, aes(!!dplyr::sym(axis_name[1]), !!dplyr::sym(axis_name[2]))) +
    ggiraph::geom_point_interactive(aes(fill=Group,tooltip = paste0(primary,'\n','x: ',round(!!dplyr::sym(axis_name[1]),2),'\n','y: ',round(!!dplyr::sym(axis_name[2]),2))),size=8,pch = 21)+
    scale_fill_manual(values=col_values,name = "Group") +
    xlab(paste(axis_name[1], "( ",pc1,"%"," )",sep="")) +
    ylab(paste(axis_name[2], "( ",pc2,"%"," )",sep="")) +
    xlim(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range) +
    ylim(ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range) +
    theme(text=element_text(size=30))+
    #geom_vline(aes(xintercept = 0),linetype="dotted")+
    #geom_hline(aes(yintercept = 0),linetype="dotted")+
    theme(panel.background = element_rect(fill='white', colour='black'),
          panel.grid=element_blank(),
          axis.title = element_text(color='black',size=34),
          axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_text(colour='black', size=34),
          axis.title.y=element_text(colour='black', size=34),
          axis.text=element_text(colour='black',size=28),
          legend.title=element_text(size = 24,face = "bold"),
          legend.text=element_text(size=20),
          legend.key=element_blank(),legend.position = c('left'),
          legend.background = element_rect(colour = "black"),
          legend.key.height=unit(1,"cm")) +
    guides(fill = guide_legend(ncol = 1))



  if (axis_num >= 3) {
    p3 <- ggplot(plotdata,aes(Group,!!dplyr::sym(axis_name[3]))) + scale_fill_manual(values=col_values) +
      geom_boxplot(aes(fill = Group),outlier.colour = NA) +
      ggpubr::stat_compare_means(comparisons = compare,method=method) +ggiraph::geom_point_interactive(aes(tooltip = paste0(primary,' : ',round(!!dplyr::sym(axis_name[3]),2))),position = "jitter")+
      theme_bw()+
      theme(axis.ticks.length = unit(0.4,"lines"),
            axis.ticks = element_line(color='black'),
            axis.line = element_line(colour = "black"),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_text(colour='black',size=20,angle = 45,
                                     vjust = 1,hjust = 1,face = "bold"),
            axis.text.y=element_blank(),
            legend.position = "none")

    p13<-ggplot(plotdata, aes(!!dplyr::sym(axis_name[1]), !!dplyr::sym(axis_name[3]))) +
      ggiraph::geom_point_interactive(aes(fill=Group,tooltip = paste0(primary,'\n','x: ',round(!!dplyr::sym(axis_name[1]),2),'\n','y: ',round(!!dplyr::sym(axis_name[3]),2))),size=8,pch = 21)+
      scale_fill_manual(values=col_values,name = "Group")+
      xlab(paste(axis_name[1], "( ",pc1,"%"," )",sep="")) +
      ylab(paste(axis_name[3], "( ",pc3,"%"," )",sep="")) +
      xlim(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range) +
      ylim(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range) +
      theme(text=element_text(size=30))+
      theme(panel.background = element_rect(fill='white', colour='black'),
            panel.grid=element_blank(),
            axis.title = element_text(color='black',size=34),
            axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
            axis.line = element_line(colour = "black"),
            axis.title.x=element_text(colour='black', size=34),
            axis.title.y=element_text(colour='black', size=34),
            axis.text=element_text(colour='black',size=28),
            legend.title=element_text(size = 24,face = "bold"),
            legend.text=element_text(size=20),
            legend.key=element_blank(),legend.position = c('left'),
            legend.background = element_rect(colour = "black"),
            legend.key.height=unit(1,"cm")) +
      guides(fill = guide_legend(ncol = 1))

    p23<-ggplot(plotdata, aes(!!dplyr::sym(axis_name[2]), !!dplyr::sym(axis_name[3]))) +
      ggiraph::geom_point_interactive(aes(fill=Group,tooltip = paste0(primary,'\n','x: ',round(!!dplyr::sym(axis_name[2]),2),'\n','y: ',round(!!dplyr::sym(axis_name[3]),2))),size=8,pch = 21)+
      scale_fill_manual(values=col_values,name = "Group")+
      xlab(paste(axis_name[2], "( ",pc2,"%"," )",sep="")) +
      ylab(paste(axis_name[3], "( ",pc3,"%"," )",sep="")) +
      xlim(ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range) +
      ylim(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range) +
      theme(text=element_text(size=30))+
      theme(panel.background = element_rect(fill='white', colour='black'),
            panel.grid=element_blank(),
            axis.title = element_text(color='black',size=34),
            axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
            axis.line = element_line(colour = "black"),
            axis.title.x=element_text(colour='black', size=34),
            axis.title.y=element_text(colour='black', size=34),
            axis.text=element_text(colour='black',size=28),
            legend.title=element_text(size = 24,face = "bold"),
            legend.text=element_text(size=20),
            legend.key=element_blank(),legend.position = c('left'),
            legend.background = element_rect(colour = "black"),
            legend.key.height=unit(1,"cm")) +
      guides(fill = guide_legend(ncol = 1))

  }


  if (!is.null(ellipse)) {
    p12 <-p12+ggplot2::stat_ellipse(aes(color=Group,group=Group),level=ellipse)+scale_colour_manual(values=palette)+guides(colour = "none")
    if (axis_num >= 3) {
      p13 <-p13+ggplot2::stat_ellipse(aes(color=Group,group=Group),level=ellipse)+scale_colour_manual(values=palette)+guides(colour = "none")
      p23 <-p23+ggplot2::stat_ellipse(aes(color=Group,group=Group),level=ellipse)+scale_colour_manual(values=palette)+guides(colour = "none")
    }
  }


 # set the warning
 assay_data  <- EMPT %>%
    .get.assay.EMPT() %>% tibble::column_to_rownames('primary')

 check_dim <- dim(assay_data)[1] * dim(assay_data)[2]
 if (force_adonis == F) {
   if (check_dim > 8.1e+07) {
   message_wrap("Inputting large-scale data may lead to extended computation time for adonis. If necessary, please enable the force_adonis = TRUE.")
   }
 }else {
   check_dim <- 1
 } 

 
 if (check_dim > 8.1e+07) {
  if (!is.null(.get.algorithm.EMPT(EMPT))) {
      distance <- .get.algorithm.EMPT(EMPT)
  } else {
      distance <- 'NULL'
  }
  
  dimension_method <- .get.method.EMPT(EMPT)
  p5 <- ggplot() +
  geom_text(aes(x = -0.5,y = 0.6,
                label = paste("\nReduce Dimension:\nmethod = ",
                              dimension_method,
                              "\ndistance = ",distance)),
            size = 5.5) +
  theme_bw() +
  xlab("") + ylab("") +
  theme(panel.grid=element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

 }else {
  #PERMANOVA分析
  if (is.null(distance_for_adonis)) {
    if (!is.null(.get.algorithm.EMPT(EMPT))) {
      distance_for_adonis <- .get.algorithm.EMPT(EMPT)
    } else {
      distance_for_adonis <- 'bray'
    }
  }else{
    distance_for_adonis <- distance_for_adonis
  }

  set.seed(seed)
  adonis_data <- .get.assay.EMPT(EMPT) %>% tibble::column_to_rownames('primary')
  adonis_result <- vegan::adonis2(adonis_data~Group,data = mapping,method = distance_for_adonis)
  p5 <- ggplot() +
    geom_text(aes(x = -0.5,y = 0.6,
                  label = paste(distance_for_adonis,
                                "\nPERMANOVA:\ndf = ",
                                adonis_result$Df[1],
                                "\nR2 = ",round(adonis_result$R2[1],4),
                                "\np-value = ",adonis_result$`Pr(>F)`[1],sep = "")),
              size = 7) +
    theme_bw() +
    xlab("") + ylab("") +
    theme(panel.grid=element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())

 }



  #图像拼接-使用patchwork包将4幅图拼在一起
  p12 <- p1 + p5 + p12 + p2 +
    patchwork::plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)
  set.seed(seed)
  p12_html=ggiraph::girafe(code = print(p12),width_svg = width,height_svg = height)
  if (axis_num >= 3) {
  p23 <- p2_r + p5 + p23 + p3 +
    patchwork::plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)
  p13 <- p1 + p5 + p13 + p3 +
    patchwork::plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)
  set.seed(seed)
  p13_html=ggiraph::girafe(code = print(p13),width_svg = width,height_svg = height)
  set.seed(seed)
  p23_html=ggiraph::girafe(code = print(p23),width_svg = width,height_svg = height)
  }


  #存储数据
  deposit$pic$p12=p12
  deposit$html$p12_html=p12_html
  deposit$pc_data=plotdata
  if (axis_num >= 3) {
    deposit$pic$p23=p23
    deposit$pic$p13=p13
    deposit$html$p13_html=p13_html
    deposit$html$p23_html=p23_html
  }



  .get.plot_deposit.EMPT(EMPT,info = 'dimension_analysis_scatterplot') <- deposit
  .get.algorithm.EMPT(EMPT) <- 'dimension_analysis_scatterplot'
  .get.info.EMPT(EMPT) <- 'EMP_dimension_analysis_scatterplot'
  .get.plot_specific.EMPT(EMPT) <- show
  .get.history.EMPT(EMPT) <- call
  class(EMPT) <- 'EMP_dimension_analysis_scatterplot'
  return(EMPT)
}
