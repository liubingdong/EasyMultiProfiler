.show_EMP_dimension_analysis_scatterplot <- function(obj,plot) {
  result <- .get.plot_deposit.EMPT(obj,info = 'EMP_dimension_analysis_scatterplot')
  switch(plot,
         "p12" = print(result$pic$p12),
         "p23" = print(result$pic$p23),
         "p13" = print(result$pic$p13),
         "p12html" = print(result$html$p12_html),
         "p23html" = print(result$html$p23_html),
         "p13html" = print(result$html$p13_html)
  )
}



#' @param obj EMPT object
#' @param seed An interger. Set the random seed for the adonis permutations.
#' @param group_level A string vector. Set the group order in the plot.
#' @param show A character string include pic (default), html. This could display graphical results on 3 axes. eg. p12,p12html,p23,p23html
#' @param distance_for_adonis A character string.Set the distance for adonis. Detailed in \code{\link[vegan]{adonis}}.
#' @param estimate_group A character string. Select the colname in the coldata to compare the data in the statistical test.
#' @param paired_group  A character string. Variable name corresponding to paired primary or sample.
#' @param paired_box_line A boolean. Whether show the paired line in the boxplot when paired test activated.(defalut:TRUE)
#' @param paired_dot_line A boolean. Whether show the paired line in the scatterplot when paired test activated.(defalut:TRUE)
#' @param palette A series of character string. Color palette.
#' @param method A character string. The name of the statistical test that is applied to barplot columns (default:wilcox.test).
#' @param key_samples A series of character string. To highlight your interested samples.
#' @param dot_size A numeric. Set the dot size.
#' @param box_width A numeric. Set the box width.
#' @param box_alpha A numeric. Set the box alpha.
#' @param step_increase A numeric vector with the increase in fraction of total height for every additional comparison to minimize overlap.
#' @param ref.group a character string specifying the reference group. If specified, for a given grouping variable, each of the group levels will be compared to the reference group (i.e. control group).
#' @param comparisons A list of length-2 vectors. The entries in the vector are either the names of 2 values on the x-axis or the 2 integers that correspond to the index of the columns of interest.(default:NULL)
#' @param ellipse A number from 0 to 1. Set the ellipse in the plot.
#' @param html_width An interger. Set the html width.
#' @param html_height An interger. Set the html height.
#' @param force_adonis Force the function run adnois analysis always.(default:FALSE)
#' @param adonis_permutations Permutations for the adonis2.(default:999)
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
#' @param ... Additional parameters for adjust the boxplot, see also \code{\link[ggsignif]{geom_signif}}.
#' @rdname EMP_scatterplot
#' @return EMPT object
EMP_scatterplot.EMP_dimension_analysis <- function(obj,seed=123,group_level='default',
                                           show='p12',distance_for_adonis=NULL,force_adonis=FALSE,adonis_permutations=999,
                                           estimate_group=NULL,palette=NULL,
                                           paired_group=NULL,paired_box_line=TRUE,paired_dot_line=TRUE,
                                           method='wilcox.test',key_samples = NULL,
                                           dot_size=8,box_width=NULL,box_alpha=1,
                                           step_increase=0.1,ref.group=NULL,comparisons=NULL,
                                           ellipse = NULL,html_width=15,html_height=15,use_cached=TRUE,...) {
  call <- match.call()

 
  if (use_cached == FALSE) {
    memoise::forget(.EMP_decostand_m) %>% invisible()
  }
  
  EMPT <- .EMP_scatterplot.EMP_dimension_analysis_m(obj=obj,seed=seed,group_level=group_level,
                                           show=show,distance_for_adonis=distance_for_adonis,force_adonis=force_adonis,adonis_permutations=adonis_permutations,
                                           estimate_group=estimate_group,palette=palette,
                                           paired_group=paired_group,paired_box_line=paired_box_line,paired_dot_line=paired_dot_line,
                                           method=method,key_samples = key_samples,
                                           dot_size=dot_size,box_width=box_width,box_alpha=box_alpha,
                                           step_increase=step_increase,ref.group=ref.group,comparisons=comparisons,
                                           ellipse = ellipse,html_width=html_width,html_height=html_height,...)
  .get.history.EMPT(EMPT) <- call
  return(EMPT)
}



#' @importFrom ggsignif geom_signif
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 geom_boxplot
#' @importFrom patchwork plot_layout
#' @importFrom vegan adonis2
.EMP_scatterplot.EMP_dimension_analysis  <- function(obj,seed=123,group_level='default',
                                           show='p12',distance_for_adonis=NULL,force_adonis=FALSE,adonis_permutations=999,
                                           estimate_group=NULL,palette=NULL,
                                           paired_group=NULL,paired_box_line=TRUE,paired_dot_line=TRUE,
                                           method='wilcox.test',key_samples = NULL,
                                           dot_size=8,box_width=NULL,box_alpha=1,
                                           step_increase=0.1,ref.group=NULL,comparisons=NULL,
                                           ellipse = NULL,html_width=15,html_height=15,...){
  primary <- Group <- NULL
  #call <- match.call()
  
  # just change to unified naming.
  if (is(obj,"EMPT")) {
    EMPT <- obj
  }else{
    stop("please check the input data for EMP_scatterplot.EMP_dimension_analysis!")
  }

  deposit <- list()

  estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)

  if (is.null(palette)) {
    col_values <- .get.palette.EMPT(EMPT)
  }else {
    col_values <- palette
  }

  EMPT %<>% .group_level_modified(estimate_group = estimate_group,group_level = group_level)
  mapping <- .get.mapping.EMPT(EMPT) %>% dplyr::select(primary,dplyr::any_of(c(!!estimate_group,!!paired_group)))

  length=length(unique(as.character(mapping$primary)))
  times1=length%/%8
  res1=length%%8
  times2=length%/%5
  res2=length%%5
  col1=rep(1:8,times1)
  col=c(col1,1:res1)
  pich1=rep(c(21:24),times2)
  pich=c(pich1,15:(15+res2))

  data <- .get.result.EMPT(EMPT,info='EMP_dimension_analysis')[['dimension_coordinate']] %>% suppressMessages()

  if (!is.null(paired_group)) {
    colnames(mapping) <- c("primary","Group","paired_group")
    plotdata <- dplyr::left_join(data,mapping,by='primary')
    ## check the missing value in the group label
    if(any(is.na(plotdata[['paired_group']]))) {
      stop('Column ',paired_group,' has beed deteced missing value, please check and filter them!')
    }
    ## check the paired sample and size
    check_paired_result <- .check_group_consistency(data=plotdata,group_col='Group',patient_col='paired_group')
    if(check_paired_result$all_patients_equal == FALSE) {
      stop(check_paired_result$message)
    }    
    plotdata <- plotdata %>% dplyr::arrange(Group,paired_group) 
  }else{
    colnames(mapping) <- c("primary","Group")
    plotdata <- dplyr::left_join(data,mapping,by='primary')
  }


  ## Due to the purpose of the module, missing value in the group label is not allowed!
  if(any(is.na(plotdata$Group))) {
    stop('Column ',estimate_group,' has beed deteced missing value, this sample should be filtered by EMP_filter firstly!')
  }

  axis_name <- colnames(data)[-1]
  axis_num <- length(unique(axis_name))

  dimension_axis_result <- .get.result.EMPT(EMPT,info='EMP_dimension_analysis')[['dimension_axis']] %>% suppressMessages()

  if (!is.null(dimension_axis_result)) {
    pc1 <-dimension_axis_result[1]
    pc2 <-dimension_axis_result[2]
    try(pc3 <-dimension_axis_result[3] %>% suppressMessages(),silent = T)
    pc1_text <- paste(axis_name[1], "( ",pc1,"%"," )",sep="")
    pc2_text <- paste(axis_name[2], "( ",pc2,"%"," )",sep="")
    try(pc3_text <- paste(axis_name[3], "( ",pc3,"%"," )",sep=""),silent = T)
  }else{
    pc1_text <- axis_name[1]
    pc2_text <- axis_name[2]
    try(pc3_text <- axis_name[3],silent = T)
  }


  #PC1和PC2的显著性检验(PC1,PC2,PC3进行组间差异检验)
  if(!method %in% c('t.test','wilcox.test')){
    stop("method should be t.test or wilcox.test! ")
  }

  # choose the compare group
  group_name <- unique(plotdata[['Group']])
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

  #相须图绘制
  if (!is.null(paired_group)) {

    p1 <- (ggplot(plotdata,aes(Group,!!dplyr::sym(axis_name[1]))) +
      geom_boxplot(aes(fill = Group),outlier.colour = NA,width=box_width,alpha=box_alpha) + scale_fill_manual(values=col_values) +
      ggsignif::geom_signif(comparisons = comparisons,test = method,test.args=list(paired=TRUE),step_increase = step_increase,...) +
      coord_flip() + 
      geom_line(aes(group = paired_group), color = 'gray', lwd = 0.5) +
      ggiraph::geom_point_interactive(aes(tooltip = paste0(primary,' : ',round(!!dplyr::sym(axis_name[1]),2))))+
      theme_bw() +
      theme(axis.ticks.length = unit(0.4,"lines"),
            axis.ticks = element_line(color='black'),
            axis.line = element_line(colour = "black"),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_text(colour='black',size=20,face = "bold"),
            axis.text.x=element_blank(),
            legend.position = "none") ) |> suppressWarnings()
    
    p2 <- (ggplot(plotdata,aes(Group,!!dplyr::sym(axis_name[2]))) +
      geom_boxplot(aes(fill = Group),outlier.colour = NA,width=box_width,alpha=box_alpha) + scale_fill_manual(values=col_values) +
      ggsignif::geom_signif(comparisons = comparisons,test = method,test.args=list(paired=TRUE),step_increase = step_increase,...) + 
      geom_line(aes(group = paired_group), color = 'gray', lwd = 0.5) +
      ggiraph::geom_point_interactive(aes(tooltip = paste0(primary,' : ',round(!!dplyr::sym(axis_name[2]),2)))) +
      theme_bw()+
      theme(axis.ticks.length = unit(0.4,"lines"),
            axis.ticks = element_line(color='black'),
            axis.line = element_line(colour = "black"),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_text(colour='black',size=20,angle = 45,
                                     vjust = 1,hjust = 1,face = "bold"),
            axis.text.y=element_blank(),
            legend.position = "none") ) |> suppressWarnings()

    p2_r <- (ggplot(plotdata,aes(Group,!!dplyr::sym(axis_name[2]))) +
      geom_boxplot(aes(fill = Group),outlier.colour = NA,width=box_width,alpha=box_alpha) + scale_fill_manual(values=col_values)+
      ggsignif::geom_signif(comparisons = comparisons,test = method,test.args=list(paired=TRUE),step_increase = step_increase,...) + coord_flip() + 
      geom_line(aes(group = paired_group), color = 'gray', lwd = 0.5) +
      ggiraph::geom_point_interactive(aes(tooltip = paste0(primary,' : ',round(!!dplyr::sym(axis_name[2]),2)))) +
      theme_bw()+
      theme(axis.ticks.length = unit(0.4,"lines"),
            axis.ticks = element_line(color='black'),
            axis.line = element_line(colour = "black"),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_text(colour='black',size=20,face = "bold"),
            axis.text.x=element_blank(),
            legend.position = "none") ) |> suppressWarnings()

    if (axis_num == 3) {
      p3 <- (ggplot(plotdata,aes(Group,!!dplyr::sym(axis_name[3]))) + scale_fill_manual(values=col_values) +
        geom_boxplot(aes(fill = Group),outlier.colour = NA,width=box_width,alpha=box_alpha) +
        ggsignif::geom_signif(comparisons = comparisons,test = method,test.args=list(paired=TRUE),step_increase = step_increase,...) +
        geom_line(aes(group = paired_group), color = 'gray', lwd = 0.5) + 
        ggiraph::geom_point_interactive(aes(tooltip = paste0(primary,' : ',round(!!dplyr::sym(axis_name[3]),2)))) +
        theme_bw()+
        theme(axis.ticks.length = unit(0.4,"lines"),
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.x=element_text(colour='black',size=20,angle = 45,
                                       vjust = 1,hjust = 1,face = "bold"),
              axis.text.y=element_blank(),
              legend.position = "none") ) |> suppressWarnings()
    }

    if (paired_box_line == FALSE) {
      p1[['layers']] <- p1[['layers']][-3]
      p2[['layers']] <- p2[['layers']][-3]
      p2_r[['layers']] <- p2_r[['layers']][-3]
      try(p3[['layers']] <- p3[['layers']][-3],silent=TRUE)
    }
 
  }else{
    p1 <- (ggplot(plotdata,aes(Group,!!dplyr::sym(axis_name[1]))) +
      geom_boxplot(aes(fill = Group),outlier.colour = NA,width=box_width,alpha=box_alpha) +scale_fill_manual(values=col_values)+
      ggsignif::geom_signif(comparisons = comparisons,test = method,step_increase = step_increase,...)+
      coord_flip() +ggiraph::geom_point_interactive(aes(tooltip = paste0(primary,' : ',round(!!dplyr::sym(axis_name[1]),2))),position = "jitter")+
      theme_bw()+
      theme(axis.ticks.length = unit(0.4,"lines"),
            axis.ticks = element_line(color='black'),
            axis.line = element_line(colour = "black"),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_text(colour='black',size=20,face = "bold"),
            axis.text.x=element_blank(),
            legend.position = "none") ) |> suppressWarnings()

    p2 <- (ggplot(plotdata,aes(Group,!!dplyr::sym(axis_name[2]))) +
      geom_boxplot(aes(fill = Group),outlier.colour = NA,width=box_width,alpha=box_alpha) +scale_fill_manual(values=col_values)+
      ggsignif::geom_signif(comparisons = comparisons,test = method,step_increase = step_increase,...)+ggiraph::geom_point_interactive(aes(tooltip = paste0(primary,' : ',round(!!dplyr::sym(axis_name[2]),2))),position = "jitter")+
      theme_bw()+
      theme(axis.ticks.length = unit(0.4,"lines"),
            axis.ticks = element_line(color='black'),
            axis.line = element_line(colour = "black"),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_text(colour='black',size=20,angle = 45,
                                     vjust = 1,hjust = 1,face = "bold"),
            axis.text.y=element_blank(),
            legend.position = "none") ) |> suppressWarnings()

    p2_r <- (ggplot(plotdata,aes(Group,!!dplyr::sym(axis_name[2]))) +
      geom_boxplot(aes(fill = Group),outlier.colour = NA,width=box_width,alpha=box_alpha) +scale_fill_manual(values=col_values)+
      ggsignif::geom_signif(comparisons = comparisons,test = method,step_increase = step_increase,...)+coord_flip() +ggiraph::geom_point_interactive(aes(tooltip = paste0(primary,' : ',round(!!dplyr::sym(axis_name[2]),2))),position = "jitter")+
      theme_bw()+
      theme(axis.ticks.length = unit(0.4,"lines"),
            axis.ticks = element_line(color='black'),
            axis.line = element_line(colour = "black"),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_text(colour='black',size=20,face = "bold"),
            axis.text.x=element_blank(),
            legend.position = "none") ) |> suppressWarnings()

    if (axis_num == 3) {
      p3 <- (ggplot(plotdata,aes(Group,!!dplyr::sym(axis_name[3]))) + scale_fill_manual(values=col_values) +
        geom_boxplot(aes(fill = Group),outlier.colour = NA,width=box_width,alpha=box_alpha) +
        ggsignif::geom_signif(comparisons = comparisons,test = method,step_increase = step_increase,...) +ggiraph::geom_point_interactive(aes(tooltip = paste0(primary,' : ',round(!!dplyr::sym(axis_name[3]),2))),position = "jitter")+
        theme_bw()+
        theme(axis.ticks.length = unit(0.4,"lines"),
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.x=element_text(colour='black',size=20,angle = 45,
                                       vjust = 1,hjust = 1,face = "bold"),
              axis.text.y=element_blank(),
              legend.position = "none") ) |> suppressWarnings()
    }
  }


  if (!is.null(key_samples)) {
    group_info <- plotdata$Group %>% unique()
    plotdata  %<>% dplyr::mutate(Group = dplyr::if_else(primary %in% key_samples, "Key", Group))
    plotdata$Group <- factor(plotdata$Group,levels = c(group_info,'Key'))
  }

  #PCoA结果图绘制

  if (!is.null(paired_group)) {
    p12<-ggplot(plotdata, aes(!!dplyr::sym(axis_name[1]), !!dplyr::sym(axis_name[2]))) +
      geom_line(aes(group = paired_group), color = 'gray', lwd = 0.5) +
      ggiraph::geom_point_interactive(aes(fill=Group,tooltip = paste0(primary,'\n','x: ',round(!!dplyr::sym(axis_name[1]),2),'\n','y: ',round(!!dplyr::sym(axis_name[2]),2))),size=dot_size,pch = 21)+
      scale_fill_manual(values=col_values,name = "Group") +
      xlab(pc1_text) +
      ylab(pc2_text) +
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
      p13<-ggplot(plotdata, aes(!!dplyr::sym(axis_name[1]), !!dplyr::sym(axis_name[3]))) +
        geom_line(aes(group = paired_group), color = 'gray', lwd = 0.5) +
        ggiraph::geom_point_interactive(aes(fill=Group,tooltip = paste0(primary,'\n','x: ',round(!!dplyr::sym(axis_name[1]),2),'\n','y: ',round(!!dplyr::sym(axis_name[3]),2))),size=dot_size,pch = 21)+
        scale_fill_manual(values=col_values,name = "Group")+
        xlab(pc1_text) +
        ylab(pc3_text) +
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
        geom_line(aes(group = paired_group), color = 'gray', lwd = 0.5) +
        ggiraph::geom_point_interactive(aes(fill=Group,tooltip = paste0(primary,'\n','x: ',round(!!dplyr::sym(axis_name[2]),2),'\n','y: ',round(!!dplyr::sym(axis_name[3]),2))),size=dot_size,pch = 21)+
        scale_fill_manual(values=col_values,name = "Group")+
        xlab(pc2_text) +
        ylab(pc3_text) +
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

   if (paired_dot_line == FALSE) {
      p12[['layers']] <- p12[['layers']][-1]
      try(p13[['layers']] <- p13[['layers']][-1],silent=TRUE)
      try(p23[['layers']] <- p23[['layers']][-1],silent=TRUE)
    }    

  }else{
    p12<-ggplot(plotdata, aes(!!dplyr::sym(axis_name[1]), !!dplyr::sym(axis_name[2]))) +
      ggiraph::geom_point_interactive(aes(fill=Group,tooltip = paste0(primary,'\n','x: ',round(!!dplyr::sym(axis_name[1]),2),'\n','y: ',round(!!dplyr::sym(axis_name[2]),2))),size=dot_size,pch = 21)+
      scale_fill_manual(values=col_values,name = "Group") +
      xlab(pc1_text) +
      ylab(pc2_text) +
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
      p13<-ggplot(plotdata, aes(!!dplyr::sym(axis_name[1]), !!dplyr::sym(axis_name[3]))) +
        ggiraph::geom_point_interactive(aes(fill=Group,tooltip = paste0(primary,'\n','x: ',round(!!dplyr::sym(axis_name[1]),2),'\n','y: ',round(!!dplyr::sym(axis_name[3]),2))),size=dot_size,pch = 21)+
        scale_fill_manual(values=col_values,name = "Group")+
        xlab(pc1_text) +
        ylab(pc3_text) +
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
        ggiraph::geom_point_interactive(aes(fill=Group,tooltip = paste0(primary,'\n','x: ',round(!!dplyr::sym(axis_name[2]),2),'\n','y: ',round(!!dplyr::sym(axis_name[3]),2))),size=dot_size,pch = 21)+
        scale_fill_manual(values=col_values,name = "Group")+
        xlab(pc2_text) +
        ylab(pc3_text) +
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
  }
  

  if (!is.null(ellipse)) {
    p12 <-p12+ggplot2::stat_ellipse(aes(color=Group,group=Group),level=ellipse)+scale_colour_manual(values=col_values)+guides(colour = "none")
    if (axis_num >= 3) {
      p13 <-p13+ggplot2::stat_ellipse(aes(color=Group,group=Group),level=ellipse)+scale_colour_manual(values=col_values)+guides(colour = "none")
      p23 <-p23+ggplot2::stat_ellipse(aes(color=Group,group=Group),level=ellipse)+scale_colour_manual(values=col_values)+guides(colour = "none")
    }
  }


 # set the warning
 assay_data  <- assay(EMPT) %>% t()

 check_dim <- dim(assay_data)[1] * dim(assay_data)[2]
 check_sample_num <- dim(assay_data)[1]

 if (force_adonis == FALSE) {
   if (check_dim > 8.1e+07 | check_sample_num > 500) {
    EMP_message("Large-scale data require more time for adonis.\nIf need this, please enable the force_adonis = TRUE.",color = 31,order = 1,show='message')
   }
  }else if (force_adonis == TRUE) {
   check_dim <- 1
   check_sample_num <- 1
  }else {
  stop('force_adonis must be TRUE or FALSE')
 }

 if (check_dim > 8.1e+07 | check_sample_num > 500) {
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
  adonis_data <- assay(EMPT) %>% t()
  adonis_result <- NULL
  try(spsUtil::quiet(adonis_result <- adonis2(formula=adonis_data~Group,data = mapping,
    method = distance_for_adonis,permutations = adonis_permutations),print_cat = FALSE, message = TRUE, warning = TRUE))

  if (is.null(adonis_result)) {
    stop("Function adonis failed, try distance_for_adonis = 'euclidean' again.")
  }

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
  p12_html=ggiraph::girafe(ggobj = print(p12),width_svg = html_width,height_svg = html_height)
  if (axis_num >= 3) {
  p23 <- p2_r + p5 + p23 + p3 +
    patchwork::plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)
  p13 <- p1 + p5 + p13 + p3 +
    patchwork::plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)
  set.seed(seed)
  p13_html=ggiraph::girafe(ggobj = print(p13),width_svg = html_width,height_svg = html_height)
  set.seed(seed)
  p23_html=ggiraph::girafe(ggobj = print(p23),width_svg = html_width,height_svg = html_height)
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



  .get.plot_deposit.EMPT(EMPT,info = 'EMP_dimension_analysis_scatterplot') <- deposit
  .get.algorithm.EMPT(EMPT) <- 'EMP_dimension_analysis_scatterplot'
  .get.info.EMPT(EMPT) <- 'EMP_dimension_analysis_scatterplot'
  .get.plot_specific.EMPT(EMPT) <- show
  #.get.history.EMPT(EMPT) <- call
  class(EMPT) <- 'EMP_dimension_analysis_scatterplot'
  return(EMPT)
}

#' @importFrom memoise memoise
.EMP_scatterplot.EMP_dimension_analysis_m <- memoise::memoise(.EMP_scatterplot.EMP_dimension_analysis,cache = cachem::cache_mem(max_size = 4096 * 1024^2))


