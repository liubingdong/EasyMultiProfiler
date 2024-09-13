
#' @param obj EMP object
#' @param palette 1-3 character string. Color palette. (default: steelblue, white, darkred).
#' @param show A character string. Show inluding all, sig and pvalue.(Only supported for EMP_cor_analysis and EMP_WGCNA_cor_analysis)
#' @param mytheme Modify components of a theme according to the ggplot2::theme.
#' @param clust_row A boolean. Whether the function clust the row or not. (default:FALSE) 
#' @param clust_col A boolean. Whether the function clust the row or not. (default:FALSE) 
#' @param dist_method A character string. More see stats::dist. (default: euclidean) 
#' @param clust_method A character string. More see fastcluster::hclust (default: complete) 
#' @param label_size A number. Set the label size. (default:4) 
#' @param tree_size A number between 0 and 1. Set the clust tree size. (default:0.1) 
#' @section Palettes:
#' The following palettes are available for use with these scales:
#' \describe{
#'   BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral,
#'   Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2, Set3
#'   Blues, BuGn, BuPu, GnBu, Greens, Greys, Oranges, OrRd, PuBu, 
#'   PuBuGn, PuRd, Purples, RdPu, Reds, YlGn, YlGnBu, YlOrBr, YlOrRd
#' }
#'
#' Modify the palette through the `palette` argument.
#' @rdname EMP_heatmap_plot
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 theme_minimal
#'
#'
EMP_heatmap.EMP_cor_analysis <- function(obj,palette=c("steelblue","white","darkred"),
                            clust_row=FALSE,clust_col=FALSE,dist_method='euclidean',clust_method='complete',tree_size=0.1,
                            show='all',label_size=4,mytheme = 'theme()'){
  var1 <- var2 <- coefficient <- EMP <- NULL
  call <- match.call()
  if (is(obj,"EMP")) {
    EMP <- obj
  }else{
    stop('Please check the input data for EMP_heatmap.EMP_cor_analysis!')
  }  

  result <- .get.result.EMP(EMP,info = 'EMP_cor_analysis')

  experiment_name <- result[["cor_info"]]
  if (length(experiment_name) == 1) {
    experiment_name <- c(experiment_name,experiment_name)
  }

  df <- result[['cor_p']]
  df_r <- result[['correlation']]

  var1_raw_order <- rownames(df_r)
  var2_raw_order <- colnames(df_r) %>% rev()

  if (any(is.na(df))) {
    stop("The NA value has been detected, please check the EMP_cor_analysis result!")
  }

  if (clust_col == TRUE) {
    var1_clust <- fastcluster::hclust(dist(df_r,method=dist_method),method = clust_method)
    var1_order <- df_r[var1_clust$order,] %>% rownames()  
    df$var1 <- factor(df$var1,levels=var1_order)
  }else{
    df$var1 <- factor(df$var1,levels=var1_raw_order)
  }
  
  if (clust_row == TRUE) {
    var2_clust <- fastcluster::hclust(dist(t(df_r),method=dist_method),method = clust_method)  
    var2_order <- df_r[,var2_clust$order] %>% colnames()     
    df$var2 <- factor(df$var2,levels=var2_order)
  }else{
    df$var2 <- factor(df$var2,levels=var2_raw_order)
  } 

  ra<-abs(df$pvalue)
  NN<-nrow(df)

  prefix1<-rep('',NN)
  for(i in 1:NN){
    prefix1[i] <- paste0('(',ra[i],')')
  }

  prefix2<-rep('',NN)
  for(i in 1:NN){
    if(ra[i]<=0.1)    prefix2[i] <- '.'
    if(ra[i]<=0.05)   prefix2[i] <- '*'
    if(ra[i]<=0.01)   prefix2[i] <- '**'
    if(ra[i]<=0.001)  prefix2[i] <- '***'
  }
  if (show == 'all') {
    ra<-paste(df$coefficient,prefix1,sep='\n')
    ra<-paste(ra,prefix2,sep=' ')
  }else if (show == 'sig') {
    ra<-paste(df$coefficient,prefix2,sep='\n')
  }else if (show == 'pvalue') {
    ra<-paste(df$coefficient,prefix1,sep='\n')
  }

  df$rp<-ra
  df.label<-df$rp


  if (length(palette) >= 3) {
    p1 <- ggplot(df,aes(x=var1,y=var2,fill=coefficient,label=df.label))+
      geom_tile(color = "white") +
      geom_text(size=label_size) + scale_fill_steps2(low = palette[1], mid=palette[2],high = palette[3],show.limits = T) + labs(x=experiment_name[1],y=experiment_name[2],fill='coefficient') +
      theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
      #guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
      eval(parse(text = paste0(mytheme)))
  }else if(length(palette) == 2){
    p1 <- ggplot(df,aes(x=var1,y=var2,fill=coefficient,label=df.label))+
      geom_tile(color = "white") +
      geom_text(size=label_size) + scale_fill_steps(low = palette[1],high = palette[2],show.limits = T) + labs(x=experiment_name[1],y=experiment_name[2],fill='coefficient') +
      theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
      #guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
      eval(parse(text = paste0(mytheme)))
  }else if(length(palette) == 1){
    check_palette <- palette %in% c('BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral',
                                    'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3',
                                    'Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd', 'PuBu',
                                    'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd')
    if (check_palette) {
      p1 <- ggplot(df,aes(x=var1,y=var2,fill=coefficient,label=df.label))+
        geom_tile(color = "white") +
        geom_text(size=label_size) + scale_fill_distiller(palette=palette) +  labs(x=experiment_name[1],y=experiment_name[2],fill='coefficient') +
        theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
        #guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
        eval(parse(text = paste0(mytheme)))

    }else{
      p1 <- ggplot(df,aes(x=var1,y=var2,fill=coefficient,label=df.label))+
        geom_tile(color = "white") +
        geom_text(size=label_size) + scale_fill_steps(low = 'white',high = palette) + labs(x=experiment_name[1],y=experiment_name[2],fill='coefficient') +
        theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
        #guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
        eval(parse(text = paste0(mytheme)))
    }
  }

  if (clust_row == TRUE) {
    row_tree <- ggtree::ggtree(var2_clust,layout = "rectangular",branch.length = "none")
    p1 <- p1 %>% aplot::insert_left(row_tree,width = tree_size)
  }
  
  if (clust_col == TRUE) {
    col_tree <- ggtree::ggtree(var1_clust,branch.length = "none") + ggtree::layout_dendrogram()
    p1 <- p1 %>% aplot::insert_top(col_tree,height = tree_size)
  } 

  .get.plot_deposit.EMP(EMP,info='EMP_cor_heatmap') <- p1
  .get.info.EMP(EMP) <- 'EMP_cor_heatmap'
  .get.history.EMP(EMP) <- call
  class(EMP) <- 'EMP_cor_heatmap'
  return(EMP)
}



#' @param obj EMPT or EMP object
#' @param palette 1-3 character string. Color palette. (default: steelblue, white, darkred)
#' @param show A character string. Show inluding all, sig and pvalue.
#' @param mytheme Modify components of a theme according to the ggplot2::theme.
#' @param clust_row A boolean. Whether the function clust the row or not. (default:FALSE) 
#' @param clust_col A boolean. Whether the function clust the row or not. (default:FALSE) 
#' @param dist_method A character string. More see stats::dist. (default: euclidean) 
#' @param clust_method A character string. More see fastcluster::hclust (default: complete) 
#' @param label_size A number. Set the label size. (default:4) 
#' @param tree_size A number between 0 and 1. Set the clust tree size. (default:0.1) 
#' @rdname EMP_heatmap_plot
#' @importFrom dplyr mutate
#'
#'
EMP_heatmap.WGCNA <- function(obj,palette=c("steelblue","white","darkred"),
                                clust_row=FALSE,clust_col=FALSE,dist_method='euclidean',clust_method='complete',tree_size=0.1,
                                show='all',label_size=4,mytheme = 'theme()'){
  WGCNA_color <- WGCNA_module_elements <- `.` <- var2 <- var1 <- coefficient <- NULL
  call <- match.call()
  if (is(obj,"EMP")) {
    WGCNA_cluster_result <- .get.result.EMPT(obj@ExperimentList[[1]],info = 'EMP_WGCNA_cluster_analysis')
    result <- obj@deposit[['WGCNA_cor_analysis_result']]
  }else if(is(obj,'EMPT')) {
    result <- .get.result.EMPT(obj,info = 'EMP_WGCNA_cor_analysis')
    WGCNA_cluster_result <- .get.result.EMPT(obj,info = 'EMP_WGCNA_cluster_analysis')
  }else {
    stop('Please check the input data for EMP_heatmap.WGCNA!')
  }

  experiment_name <- result[["cor_info"]]
  if (length(experiment_name) == 1) {
    experiment_name <- c(experiment_name,experiment_name)
  }

  df <- result[['cor_p']]
  df_r <- result[['correlation']]

  var1_raw_order <- rownames(df_r)
  var2_raw_order <- colnames(df_r)

  if (any(is.na(df))) {
    stop("The NA value has been detected, please check the EMP_WGCNA_cor_analysis result!")
  }

  if (clust_col == TRUE) {
    var1_clust <- fastcluster::hclust(dist(df_r,method=dist_method),method = clust_method)
    var1_order <- df_r[var1_clust$order,] %>% rownames()  
    df$var1 <- factor(df$var1,levels=var1_order)
  }else{
    df$var1 <- factor(df$var1,levels=var1_raw_order)
  }
  
  ra<-abs(df$pvalue)
  NN<-nrow(df)

  prefix1<-rep('',NN)
  for(i in 1:NN){
    prefix1[i] <- paste0('(',ra[i],')')
  }

  prefix2<-rep('',NN)
  for(i in 1:NN){
    if(ra[i]<=0.1)    prefix2[i] <- '.'
    if(ra[i]<=0.05)   prefix2[i] <- '*'
    if(ra[i]<=0.01)   prefix2[i] <- '**'
    if(ra[i]<=0.001)  prefix2[i] <- '***'
  }
  if (show == 'all') {
    ra<-paste(df$coefficient,prefix1,sep='\n')
    ra<-paste(ra,prefix2,sep=' ')
  }else if (show == 'sig') {
    ra<-paste(df$coefficient,prefix2,sep='\n')
  }else if (show == 'pvalue') {
    ra<-paste(df$coefficient,prefix1,sep='\n')
  }

  df$rp<-ra
  df.label<-df$rp

 ## add module elements nums to y label
 net <- WGCNA_cluster_result[['WGCNA_cluster_result']]
 feature_modules <-WGCNA_cluster_result[['WGCNA_cluster_df']] %>%
    dplyr::select(WGCNA_color,WGCNA_module_elements) %>%
    dplyr::mutate(var2=paste0('ME',WGCNA_color)) %>%
    dplyr::select(-WGCNA_color) %>% dplyr::distinct()
 df %<>% dplyr::left_join(.,feature_modules,by='var2') %>%
  dplyr::mutate(var2 = paste(var2,WGCNA_module_elements,sep='\n'))


  if (clust_row == TRUE) {
    df_r2 <- df %>%  dplyr::select(var1,var2,coefficient) %>%
      tidyr::pivot_wider(names_from = 'var2',values_from = 'coefficient') %>%
      tibble::column_to_rownames('var1')
    
    var2_clust <- fastcluster::hclust(dist(t(df_r2),method=dist_method),method = clust_method)  
    var2_order <- df_r2[,var2_clust$order] %>% colnames()     
    df$var2 <- factor(df$var2,levels=var2_order)
  }


  if (length(palette) >= 3) {
    p1 <- ggplot(df,aes(x=var1,y=var2,fill=coefficient,label=df.label))+
      geom_tile(color = "white") +
      geom_text(size=label_size) + scale_fill_steps2(low = palette[1], mid=palette[2],high = palette[3],show.limits = T) + labs(x=experiment_name[2],y=experiment_name[1],fill='coefficient') +
      theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
      #guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
      eval(parse(text = paste0(mytheme)))
  }else if(length(palette) == 2){
    p1 <- ggplot(df,aes(x=var1,y=var2,fill=coefficient,label=df.label))+
      geom_tile(color = "white") +
      geom_text(size=label_size) + scale_fill_steps(low = palette[1],high = palette[2],show.limits = T) + labs(x=experiment_name[2],y=experiment_name[1],fill='coefficient') +
      theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
      #guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
      eval(parse(text = paste0(mytheme)))
  }else if(length(palette) == 1){
    check_palette <- palette %in% c('BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral',
                                    'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3',
                                    'Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd', 'PuBu',
                                    'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd')
    if (check_palette) {
      p1 <- ggplot(df,aes(x=var1,y=var2,fill=coefficient,label=df.label))+
        geom_tile(color = "white") +
        geom_text(size=label_size) + scale_fill_distiller(palette=palette) +  labs(x=experiment_name[2],y=experiment_name[1],fill='coefficient') +
        theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
        #guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
        eval(parse(text = paste0(mytheme)))

    }else{
      p1 <- ggplot(df,aes(x=var1,y=var2,fill=coefficient,label=df.label))+
        geom_tile(color = "white") +
        geom_text(size=label_size) + scale_fill_steps(low = 'white',high = palette) + labs(x=experiment_name[2],y=experiment_name[1],fill='coefficient') +
        theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
        #guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
        eval(parse(text = paste0(mytheme)))
    }
  }

  if (clust_row == TRUE) {
    row_tree <- ggtree::ggtree(var2_clust,layout = "rectangular",branch.length = "none")
    p1 <- p1 %>% aplot::insert_left(row_tree,width = tree_size)
  }
  
  if (clust_col == TRUE) {
    col_tree <- ggtree::ggtree(var1_clust,branch.length = "none") + ggtree::layout_dendrogram()
    p1 <- p1 %>% aplot::insert_top(col_tree,height = tree_size)
  } 

  if (is(obj,"EMP")) {
    EMP <- obj
    .get.plot_deposit.EMP(EMP,info='EMP_WGCNA_cor_heatmap') <- p1
    .get.info.EMP(EMP) <- 'EMP_WGCNA_cor_heatmap2'
    .get.history.EMP(EMP) <- call
    class(EMP) <- 'EMP_WGCNA_cor_heatmap2'
    return(EMP)
  }else if(is(obj,'EMPT')) {
    EMPT <- obj
    .get.plot_deposit.EMPT(EMPT,info='EMP_WGCNA_cor_heatmap') <- p1
    .get.info.EMPT(EMPT) <- 'EMP_WGCNA_cor_heatmap'
    .get.history.EMPT(EMPT) <- call
    class(EMPT) <- 'EMP_WGCNA_cor_heatmap'
    return(EMPT)
  }
}


#' @param obj EMPT or EMP object
#' @param palette 1-3 character string. Color palette. (default: steelblue, white, darkred)
#' @param rotate A boolean. Whether rotate the heatmap or not. (Only activated for EMP_assay_data)
#' @param palette 1-3 character string. Color palette. (default: steelblue, white, darkred)
#' @param clust_row A boolean. Whether the function clust the row or not. (default:FALSE) 
#' @param clust_col A boolean. Whether the function clust the row or not. (default:FALSE) 
#' @param dist_method A character string. More see stats::dist. (default: euclidean) 
#' @param clust_method A character string. More see fastcluster::hclust (default: complete) 
#' @param tree_size A number between 0 and 1. Set the clust tree size. (default:0.1) 
#' @param label_size A number. Set the label size. (default:4) 
#' @param mytheme Modify components of a theme according to the ggplot2::theme.
#' @rdname EMP_heatmap_plot
#' @importFrom forcats fct_relevel
#' @importFrom aplot insert_top
#' @importFrom aplot insert_left
#' @importFrom ggtree ggtree
#' @importFrom stats dist
#' @importFrom fastcluster hclust

EMP_heatmap.EMP_assay_data <- function(obj,palette=c("steelblue","white","darkred"),rotate=FALSE,
                                         clust_row=FALSE,clust_col=FALSE,dist_method='euclidean',clust_method='complete',tree_size=0.1,label_size=4,
                                         mytheme = 'theme()'){
  call <- match.call()

  primary <- value <- NULL
  
  if (is(obj,"EMPT")) {
    EMPT <- obj
  }else{
    stop('Please check the input data for EMP_heatmap.EMP_assay_data!')
  } 
  
  result <- .get.result.EMPT(EMPT,info = 'EMP_assay_data') %>% suppressMessages()
  
  if (any(is.na(result))) {
    stop("The NA value has been detected, please check the assay data!")
  }

  if (clust_row == TRUE | clust_col == TRUE) {
    result_clust <- result %>%
       tibble::column_to_rownames('primary')

    primary_clust <- fastcluster::hclust(dist(result_clust,method=dist_method),method = clust_method)
    feature_clust <- fastcluster::hclust(dist(t(result_clust),method=dist_method),method = clust_method)  
  
    primary_order <- result_clust[primary_clust$order,] %>% rownames() 
    feature_order <- result_clust[,feature_clust$order] %>% colnames()     
  }

  df <- result %>% tidyr::pivot_longer(
    cols =  -primary,
    names_to = 'feature',
    values_to = "value"
  )

  if (rotate == TRUE) {
    xy_name <- c('feature','primary')
    if (clust_row == TRUE) {
      df$primary <- factor(df$primary,levels = primary_order)
    }
    if (clust_col == TRUE) {
      df$feature <- factor(df$feature,levels = feature_order)      
    }
  }else{
    xy_name <- c('primary','feature')  
    if (clust_col == TRUE) {
      df$primary <- factor(df$primary,levels = primary_order)
    }
    if (clust_row == TRUE) {
      df$feature <- factor(df$feature,levels = feature_order)
    }
  }

  if (length(palette) >= 3) {
    p1 <- ggplot(df,aes(x=!!sym(xy_name[1]),y=!!sym(xy_name[2]),fill=value)) +
      geom_tile(color = "white") +
      scale_fill_steps2(low = palette[1], mid=palette[2],high = palette[3],midpoint = median(df$value),show.limits = T) + 
      xlab(NULL) + ylab(NULL) +
      theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
      #guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
      eval(parse(text = paste0(mytheme)))
  }else if(length(palette) == 2){
    p1 <- ggplot(df,aes(x=!!sym(xy_name[1]),y=!!sym(xy_name[2]),fill=value)) +
      geom_tile(color = "white") +
      scale_fill_steps(low = palette[1],high = palette[2],midpoint = median(df$value),show.limits = T) + 
      xlab(NULL) + ylab(NULL) +
      theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
      #guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
      eval(parse(text = paste0(mytheme)))
  }else if(length(palette) == 1){
    check_palette <- palette %in% c('BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral',
                                    'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3',
                                    'Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd', 'PuBu',
                                    'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd')
    if (check_palette) {
      p1 <- ggplot(df,aes(x=!!sym(xy_name[1]),y=!!sym(xy_name[2]),fill=value)) +
        geom_tile(color = "white") +
        scale_fill_distiller(palette=palette) + 
        xlab(NULL) + ylab(NULL) + 
        theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
        #guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
        eval(parse(text = paste0(mytheme)))
      
    }else{
      p1 <- ggplot(df,aes(x=!!sym(xy_name[1]),y=!!sym(xy_name[2]),fill=value)) +
        geom_tile(color = "white") +
        scale_fill_steps(low = 'white',high = palette) + 
        xlab(NULL) + ylab(NULL) +
        theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
        #guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
        eval(parse(text = paste0(mytheme)))
    }
  }
  
  if (rotate == TRUE) {
    #xy_name <- c('feature','primary')
    if (clust_row == TRUE) {
      row_tree <- ggtree::ggtree(primary_clust,layout = "rectangular",branch.length = "none")
      p1 <- p1 %>% aplot::insert_left(row_tree,width = tree_size)
    }
    if (clust_col == TRUE) {
      col_tree <- ggtree::ggtree(feature_clust,branch.length = "none") + ggtree::layout_dendrogram()
      p1 <- p1 %>% aplot::insert_top(col_tree,height = tree_size)
    }   
  }else{
    #xy_name <- c('primary','feature')  
    if (clust_row == TRUE) {
      row_tree <- ggtree::ggtree(feature_clust,layout = "rectangular",branch.length = "none")
      p1 <- p1 %>% aplot::insert_left(row_tree,width = tree_size)
    }
    if (clust_col == TRUE) {
      col_tree <- ggtree::ggtree(primary_clust,branch.length = "none") + ggtree::layout_dendrogram()
      p1 <- p1 %>% aplot::insert_top(col_tree,height = tree_size)
    }     
  }

  .get.plot_deposit.EMPT(EMPT,info='EMP_assay_heatmap') <- p1
  .get.info.EMPT(EMPT) <- 'EMP_assay_heatmap'
  .get.history.EMPT(EMPT) <- call
  class(EMPT) <- 'EMP_assay_heatmap'
  return(EMPT)
}  


