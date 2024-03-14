#' Title
#'
#' @param data wait_for_add
#' @param meta_data wait_for_add
#' @param method wait_for_add
#' @param width wait_for_add
#' @param height wait_for_add
#' @param cor_output wait_for_add
#' @param cluster_rows wait_for_add
#' @param cluster_cols wait_for_add
#' @param file_name wait_for_add
#'
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom graphics abline
#' @importFrom methods new
#' @importFrom psych corr.test
#' @importFrom dplyr if_any
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_COR_HEAT_ <-function(data,meta_data,method = 'spearman',width=10,height=10,
                         cor_output=F,cluster_rows = F,cluster_cols = F,
                         file_name='cor_plot'){
  Group <- getSig <- NULL                      
  deposit <- list()
  real_sample <- intersect(meta_data$primary,data$primary)
  data <- data[data$primary%in%real_sample, ]
  meta_data <- meta_data[meta_data$primary%in%real_sample, ]
  data=data[order(data$primary),]
  meta_data=meta_data[order(meta_data$primary),]
  try(data <- subset(data,select = -c(Group)), silent = T)
  data.corr <- psych::corr.test(data[,-1], meta_data[,-1],method = method,adjust='none')
  data.r <- data.corr$r
  data.p <- data.corr$p

  data.p %>% as.data.frame() %>% filter(if_any(everything(),~. <0.05)) %>%
    select_if(~ min(.) < 0.05) -> data.p_filter
  data.r[rownames(data.p_filter),colnames(data.p_filter)] -> data.r_filter


  sig.mat <- matrix(sapply(as.matrix(data.p_filter), getSig), nrow=nrow(data.p_filter))
  pheatmap(data.r_filter, clustering_method="average", cluster_rows=cluster_rows,cluster_cols = cluster_cols, display_numbers=sig.mat)

  if (cor_output==T) {
    pheatmap(data.r_filter, clustering_method="average", cluster_rows=cluster_rows,cluster_cols = cluster_cols, display_numbers=sig.mat,filename=paste0(file_name,'.pdf'),width=width,height=height)
    dev.off()
  }
  deposit$r <- data.r
  deposit$p <- data.p
  deposit$p_filter <-data.p_filter
  deposit$r_filter <-data.r_filter
  return(deposit)
}

#' Title
#'
#' @param EMP wait_for_add
#' @param palette wait_for_add
#' @param show wait_for_add
#' @param mytheme wait_for_add
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 theme_minimal
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_heatmap_cor <- function(EMP,palette=c("steelblue","white","darkred"),
                            show='all',mytheme = 'theme()'){

  result <- .get.result.EMP(EMP,info = 'EMP_cor_analysis')

  experiment_name <- result[["cor_info"]]
  if (length(experiment_name) == 1) {
    experiment_name <- c(experiment_name,experiment_name)
  }

  df <- result[['cor_p']]

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
      geom_text() + scale_fill_gradient2(low = palette[1], mid=palette[2],high = palette[3]) + labs(x=experiment_name[1],y=experiment_name[2],fill='coefficient') +
      theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
      guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
      eval(parse(text = paste0(mytheme)))
  }else if(length(palette) == 2){
    p1 <- ggplot(df,aes(x=var1,y=var2,fill=coefficient,label=df.label))+
      geom_tile(color = "white") +
      geom_text() + scale_fill_gradient(low = palette[1],high = palette[2]) + labs(x=experiment_name[1],y=experiment_name[2],fill='coefficient') +
      theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
      guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
      eval(parse(text = paste0(mytheme)))
  }else if(length(palette) == 1){
    check_palette <- palette %in% c('BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral',
                                    'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3',
                                    'Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd', 'PuBu',
                                    'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd')
    if (check_palette) {
      p1 <- ggplot(df,aes(x=var1,y=var2,fill=coefficient,label=df.label))+
        geom_tile(color = "white") +
        geom_text() + scale_fill_distiller(palette=palette) +  labs(x=experiment_name[1],y=experiment_name[2],fill='coefficient') +
        theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
        guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
        eval(parse(text = paste0(mytheme)))

    }else{
      p1 <- ggplot(df,aes(x=var1,y=var2,fill=coefficient,label=df.label))+
        geom_tile(color = "white") +
        geom_text() + scale_fill_gradient(low = 'white',high = palette) + labs(x=experiment_name[1],y=experiment_name[2],fill='coefficient') +
        theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
        guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
        eval(parse(text = paste0(mytheme)))
    }
  }

  return(p1)


}


#' Title
#'
#' @param x wait_for_add
#' @param palette wait_for_add
#' @param show wait_for_add
#' @param mytheme wait_for_add
#' @importFrom dplyr mutate
#' @importFrom ggplot2 guides
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_heatmap_WGCNA<- function(x,palette=c("steelblue","white","darkred"),show='all',mytheme = 'theme()'){
  call <- match.call()
  if (inherits(x,"EMP")) {
    WGCNA_cluster_result <- .get.result.EMPT(x@ExperimentList[[1]],info = 'EMP_WGCNA_cluster_analysis')
    result <- x@deposit[['WGCNA_cor_analysis_result']]
  }else if(inherits(x,'EMPT')) {
    result <- .get.result.EMPT(x,info = 'EMP_WGCNA_cor_analysis')
    WGCNA_cluster_result <- .get.result.EMPT(x,info = 'EMP_WGCNA_cluster_analysis')
  }else {
    stop('Please check the input data!')
  }

  experiment_name <- result[["cor_info"]]
  if (length(experiment_name) == 1) {
    experiment_name <- c(experiment_name,experiment_name)
  }

  df <- result[['cor_p']]

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


  if (length(palette) >= 3) {
    p1 <- ggplot(df,aes(x=var1,y=var2,fill=coefficient,label=df.label))+
      geom_tile(color = "white") +
      geom_text() + scale_fill_gradient2(low = palette[1], mid=palette[2],high = palette[3]) + labs(x=experiment_name[1],y=experiment_name[2],fill='coefficient') +
      theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
      guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
      eval(parse(text = paste0(mytheme)))
  }else if(length(palette) == 2){
    p1 <- ggplot(df,aes(x=var1,y=var2,fill=coefficient,label=df.label))+
      geom_tile(color = "white") +
      geom_text() + scale_fill_gradient(low = palette[1],high = palette[2]) + labs(x=experiment_name[1],y=experiment_name[2],fill='coefficient') +
      theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
      guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
      eval(parse(text = paste0(mytheme)))
  }else if(length(palette) == 1){
    check_palette <- palette %in% c('BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral',
                                    'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3',
                                    'Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd', 'PuBu',
                                    'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd')
    if (check_palette) {
      p1 <- ggplot(df,aes(x=var1,y=var2,fill=coefficient,label=df.label))+
        geom_tile(color = "white") +
        geom_text() + scale_fill_distiller(palette=palette) +  labs(x=experiment_name[1],y=experiment_name[2],fill='coefficient') +
        theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
        guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
        eval(parse(text = paste0(mytheme)))

    }else{
      p1 <- ggplot(df,aes(x=var1,y=var2,fill=coefficient,label=df.label))+
        geom_tile(color = "white") +
        geom_text() + scale_fill_gradient(low = 'white',high = palette) + labs(x=experiment_name[1],y=experiment_name[2],fill='coefficient') +
        theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
        guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
        eval(parse(text = paste0(mytheme)))
    }
  }


  if (inherits(x,"EMP")) {
    EMP <- x
    .get.plot_deposit.EMP(EMP,info='EMP_WGCNA_cor_heatmap') <- p1
    .get.info.EMP(EMP) <- 'EMP_WGCNA_cor_heatmap2'
    .get.history.EMP(EMP) <- call
    class(EMP) <- 'EMP_WGCNA_cor_heatmap2'
    return(EMP)
  }else if(inherits(x,'EMPT')) {
    EMPT <- x
    .get.plot_deposit.EMPT(EMPT,info='EMP_WGCNA_cor_heatmap') <- p1
    .get.info.EMPT(EMPT) <- 'EMP_WGCNA_cor_heatmap'
    .get.history.EMPT(EMPT) <- call
    class(EMPT) <- 'EMP_WGCNA_cor_heatmap'
    return(EMPT)
  }
}
