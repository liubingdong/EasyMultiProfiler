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
  var1 <- var2 <- coefficient <- NULL
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
EMP_heatmap_WGCNA <- function(x,palette=c("steelblue","white","darkred"),show='all',mytheme = 'theme()'){
  WGCNA_color <- WGCNA_module_elements <- `.` <- var2 <- var1 <- coefficient <- NULL
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
