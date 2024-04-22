
#' @param obj EMP object
#' @param palette 1-3 character string. Color palette. (default: steelblue, white, darkred)
#' @param show A character string. Show inluding all, sig and pvalue.
#' @param mytheme Modify components of a theme according to the ggplot2::theme.
#' @rdname EMP_heatmap_plot
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 theme_minimal
#'
#'
EMP_heatmap.EMP_cor_analysis <- function(obj,palette=c("steelblue","white","darkred"),
                            show='all',mytheme = 'theme()'){
  var1 <- var2 <- coefficient <- EMP <- NULL
  
  if (inherits(obj,"EMP")) {
    EMP <- obj
  }else{
    stop('Please check the input data!')
  }  

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
      geom_text() + scale_fill_steps2(low = palette[1], mid=palette[2],high = palette[3],show.limits = T) + labs(x=experiment_name[1],y=experiment_name[2],fill='coefficient') +
      theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
      #guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
      eval(parse(text = paste0(mytheme)))
  }else if(length(palette) == 2){
    p1 <- ggplot(df,aes(x=var1,y=var2,fill=coefficient,label=df.label))+
      geom_tile(color = "white") +
      geom_text() + scale_fill_steps(low = palette[1],high = palette[2],show.limits = T) + labs(x=experiment_name[1],y=experiment_name[2],fill='coefficient') +
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
        geom_text() + scale_fill_distiller(palette=palette) +  labs(x=experiment_name[1],y=experiment_name[2],fill='coefficient') +
        theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
        #guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
        eval(parse(text = paste0(mytheme)))

    }else{
      p1 <- ggplot(df,aes(x=var1,y=var2,fill=coefficient,label=df.label))+
        geom_tile(color = "white") +
        geom_text() + scale_fill_steps(low = 'white',high = palette) + labs(x=experiment_name[1],y=experiment_name[2],fill='coefficient') +
        theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
        #guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
        eval(parse(text = paste0(mytheme)))
    }
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
#' @rdname EMP_heatmap_plot
#' @importFrom dplyr mutate
#'
#'
EMP_heatmap.WGCNA <- function(obj,palette=c("steelblue","white","darkred"),show='all',mytheme = 'theme()'){
  WGCNA_color <- WGCNA_module_elements <- `.` <- var2 <- var1 <- coefficient <- NULL
  call <- match.call()
  if (inherits(obj,"EMP")) {
    WGCNA_cluster_result <- .get.result.EMPT(obj@ExperimentList[[1]],info = 'EMP_WGCNA_cluster_analysis')
    result <- obj@deposit[['WGCNA_cor_analysis_result']]
  }else if(inherits(obj,'EMPT')) {
    result <- .get.result.EMPT(obj,info = 'EMP_WGCNA_cor_analysis')
    WGCNA_cluster_result <- .get.result.EMPT(obj,info = 'EMP_WGCNA_cluster_analysis')
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
      geom_text() + scale_fill_steps2(low = palette[1], mid=palette[2],high = palette[3],show.limits = T) + labs(x=experiment_name[2],y=experiment_name[1],fill='coefficient') +
      theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
      #guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
      eval(parse(text = paste0(mytheme)))
  }else if(length(palette) == 2){
    p1 <- ggplot(df,aes(x=var1,y=var2,fill=coefficient,label=df.label))+
      geom_tile(color = "white") +
      geom_text() + scale_fill_steps(low = palette[1],high = palette[2],show.limits = T) + labs(x=experiment_name[2],y=experiment_name[1],fill='coefficient') +
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
        geom_text() + scale_fill_distiller(palette=palette) +  labs(x=experiment_name[2],y=experiment_name[1],fill='coefficient') +
        theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
        #guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
        eval(parse(text = paste0(mytheme)))

    }else{
      p1 <- ggplot(df,aes(x=var1,y=var2,fill=coefficient,label=df.label))+
        geom_tile(color = "white") +
        geom_text() + scale_fill_steps(low = 'white',high = palette) + labs(x=experiment_name[2],y=experiment_name[1],fill='coefficient') +
        theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
        #guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
        eval(parse(text = paste0(mytheme)))
    }
  }


  if (inherits(obj,"EMP")) {
    EMP <- obj
    .get.plot_deposit.EMP(EMP,info='EMP_WGCNA_cor_heatmap') <- p1
    .get.info.EMP(EMP) <- 'EMP_WGCNA_cor_heatmap2'
    .get.history.EMP(EMP) <- call
    class(EMP) <- 'EMP_WGCNA_cor_heatmap2'
    return(EMP)
  }else if(inherits(obj,'EMPT')) {
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
#' @param mytheme Modify components of a theme according to the ggplot2::theme.
#' @rdname EMP_heatmap_plot

EMP_heatmap.EMP_assay_data <- function(obj,palette=c("steelblue","white","darkred"),rotate=FALSE,
                                         mytheme = 'theme()'){
  
  primary <- value <- NULL
  
  if (inherits(obj,"EMPT")) {
    EMPT <- obj
  }else{
    stop('Please check the input data!')
  } 
  
  result <- .get.result.EMPT(EMPT,info = 'EMP_assay_data')
  
  df <- result %>% tidyr::pivot_longer(
    cols =  -primary,
    names_to = 'feature',
    values_to = "value"
  )
  
  
  if (rotate == FALSE) {
    xy_name <- c('feature','primary')
  }else{
    xy_name <- c('primary','feature')
  }
  
  if (length(palette) >= 3) {
    p1 <- ggplot(df,aes(x=!!sym(xy_name[1]),y=!!sym(xy_name[2]),fill=value)) +
      geom_tile(color = "white") +
      scale_fill_steps2(low = palette[1], mid=palette[2],high = palette[3],show.limits = T) + 
      theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
      #guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
      eval(parse(text = paste0(mytheme)))
  }else if(length(palette) == 2){
    p1 <- ggplot(df,aes(x=!!sym(xy_name[1]),y=!!sym(xy_name[2]),fill=value)) +
      geom_tile(color = "white") +
      scale_fill_steps(low = palette[1],high = palette[2],show.limits = T) + 
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
        theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
        #guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
        eval(parse(text = paste0(mytheme)))
      
    }else{
      p1 <- ggplot(df,aes(x=!!sym(xy_name[1]),y=!!sym(xy_name[2]),fill=value)) +
        geom_tile(color = "white") +
        scale_fill_steps(low = 'white',high = palette) + 
        theme_minimal() +theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
        #guides(fill = guide_colorsteps(title.position = "top",show.limits = TRUE), color="none") +
        eval(parse(text = paste0(mytheme)))
    }
  }
  
  .get.plot_deposit.EMPT(EMPT,info='EMP_assay_heatmap') <- p1
  .get.info.EMPT(EMPT) <- 'EMP_assay_heatmap'
  .get.history.EMPT(EMPT) <- call
  class(EMPT) <- 'EMP_assay_heatmap'
  return(EMPT)
}  


