
#' @import ggplot2
#' @importFrom ggplot2 geom_smooth
.EMP_fitline_plot_default.EMPT <- function(EMPT,var_select,formula = y~poly(x,1,raw = FALSE),se=FALSE,eq_size=3,show='pic',estimate_group=NULL,palette=NULL,mytheme='theme()',html_width = 5, html_height = 5,...){
  
  rlang::check_installed(c('BiocManager'), reason = 'for EMP_fitline_plot().', action = install.packages) 
  rlang::check_installed(c('ggpmisc'), reason = 'for EMP_fitline_plot().', action = BiocManager::install)     
  
  library('ggplot2')

  primary <- eq.label <- rr.label <- adj.rr.label <-p.value.label <- NULL

  assay_data <- .get.assay.EMPT(EMPT)
  coldata <- .get.mapping.EMPT(EMPT)
  
  data <- dplyr::inner_join(assay_data,coldata,by='primary') 
    
  if (is.null(palette)) {
    palette <- .get.palette.EMPT(EMPT)
  }else {
    palette = palette
  }

  # 确保方程不会被点覆盖
  ymax <- max(data[,var_select[2]])*1.1
  ymin <- min(data[,var_select[2]])*0.9
  
  data_plot <- list()
  
  if (is.null(estimate_group)) {

    data <- data %>% dplyr::filter(dplyr::if_all(dplyr::all_of(var_select), ~ !is.na(.)))

    data_plot[['pic']] <- ggplot(data=data,aes(x=!!sym(var_select[1]),y=!!sym(var_select[2]))) +
        ggpmisc::stat_poly_line(formula =formula,colour='red',se=se)+
        ggiraph::geom_jitter_interactive(aes(tooltip = paste0(primary,'\n','x: ',!!sym(var_select[1]),'\n','y: ',!!sym(var_select[2]))),
                                         position = position_jitter(height = .00000001))+
        theme_bw() +  
        xlab(paste0(var_select[1])) + ylab(paste0(var_select[2])) +
        ggpmisc::stat_poly_eq(formula = formula, 
                              aes(label = paste(after_stat(eq.label), "*\", \"*", 
                                                after_stat(rr.label), "*\", \"*", 
                                                after_stat(adj.rr.label), "*\", and \"*",
                                                after_stat(p.value.label), "*\".\"",
                                                sep = "")
                              ),size=eq_size,
                              parse = TRUE) + coord_cartesian(ylim=c(ymin,ymax)) + eval(parse(text = paste0(mytheme))) 
  }else{
    
    if(!estimate_group %in% colnames(data)){
      stop('Please check the parameter estimate_group, ',estimate_group,' is not in the data!')
    }
    
    ## check the missing value in the group label
    if(any(is.na(data[[estimate_group]]))) {
      stop('Column ',estimate_group,' has beed deteced missing value, please check and filter them!')
    }

    ymax <- ymax*1.3
    
    data <- data %>% dplyr::filter(dplyr::if_all(dplyr::all_of(c(var_select,estimate_group)), ~ !is.na(.)))

    data_plot[['pic']]  <- ggplot(data=data,aes(x=!!sym(var_select[1]),y=!!sym(var_select[2]),colour = !!sym(estimate_group)))+
      ggpmisc::stat_poly_line(formula = formula,se=se)+  
      ggiraph::geom_jitter_interactive(aes(tooltip = paste0(primary,'\n','x: ',!!sym(var_select[1]),'\n','y: ',!!sym(var_select[2]))),
                                       position = position_jitter(height = .00000001),size=3)+
      theme_bw() +  
      xlab(paste0(var_select[1])) + ylab(paste0(var_select[2])) +
      ggpmisc::stat_poly_eq(formula = formula, vstep = NULL,
                            aes(label = paste(after_stat(eq.label), "*\", \"*", 
                                              after_stat(rr.label), "*\", \"*", 
                                              after_stat(adj.rr.label), "*\", and \"*",
                                              after_stat(p.value.label), "*\".\"",
                                              sep = "")
                            ),size=eq_size,
                            parse = TRUE) + coord_cartesian(ylim=c(ymin,ymax)) + 
      scale_colour_manual(values=palette,name = estimate_group) + eval(parse(text = paste0(mytheme))) 
    
  }
  data_plot[['html']]  <- ggiraph::girafe(code = print(data_plot[['pic']] ),width = html_width,height = html_height)

  .get.plot_deposit.EMPT(EMPT,info = 'EMP_fitline_plot') <- data_plot
  .get.plot_specific.EMPT(EMPT) <- show
  .get.estimate_group.EMPT(EMPT) <- estimate_group
  EMPT@algorithm <- 'EMP_fitline_plot'
  .get.info.EMPT(EMPT) <- 'EMP_fitline_plot'
  class(EMPT) <- 'EMP_fitline_plot'
  EMPT
}


.EMP_fitline_plot_default.EMP <- function(EMP,var_select,select=NULL,formula = y~poly(x,1,raw = FALSE),se=FALSE,eq_size=3,show='pic',estimate_group=NULL,palette=NULL,mytheme='theme()',html_width = 5, html_height = 5,...){
  
  rlang::check_installed(c('BiocManager'), reason = 'for EMP_fitline_plot().', action = install.packages) 
  rlang::check_installed(c('ggpmisc'), reason = 'for EMP_fitline_plot().', action = BiocManager::install)  

  library('ggplot2')

  primary <- eq.label <- rr.label <- adj.rr.label <-p.value.label <- NULL

  experiment_name <- names(EMP@ExperimentList)
  experiment_num <- length(experiment_name)

  if (length(var_select)==1) {
    var_select <- c(var_select,var_select)
  }
  
  if (experiment_num >=2) {
    if (is.null(select)) {
      assay_data1 <- EMP@ExperimentList[[1]] %>% .get.assay.EMPT()
      assay_data2 <- EMP@ExperimentList[[2]] %>% .get.assay.EMPT()
      
      coldata1 <- EMP@ExperimentList[[1]] %>% .get.mapping.EMPT()
      
      x_data <- assay_data1 %>% dplyr::select(primary,dplyr::all_of(var_select[1]))
      y_data <- assay_data2 %>% dplyr::select(primary,dplyr::all_of(var_select[2]))

      data <- dplyr::inner_join(x_data,y_data,by='primary') %>%
        dplyr::inner_join(coldata1,by='primary')      
      
    }else{
      if (!all(select %in% experiment_name)) {
        stop("Pararmeter select in not in the ExperimentList,please check!")
      }
      if (length(select)==1) {
        select <- c(select,select)
      }   
         
      assay_data1 <- EMP@ExperimentList[[select[1]]] %>% .get.assay.EMPT()
      assay_data2 <- EMP@ExperimentList[[select[2]]] %>% .get.assay.EMPT()
      
      coldata1 <- EMP@ExperimentList[[1]] %>% .get.mapping.EMPT()
      
      x_data <- assay_data1 %>% dplyr::select(primary,dplyr::all_of(var_select[1]))
      y_data <- assay_data2 %>% dplyr::select(primary,dplyr::all_of(var_select[2]))
      data <- dplyr::inner_join(x_data,y_data,by='primary') %>%
        dplyr::inner_join(coldata1,by='primary')      
    }
  }else{
    assay_data <- .get.assay.EMPT(EMP@ExperimentList[[1]]) 
    coldata <- .get.mapping.EMPT(EMP@ExperimentList[[1]])
    data <- dplyr::inner_join(assay_data,coldata,by='primary') 
    experiment_name <- c(experiment_name,experiment_name)
    if (var_select[[1]] == var_select[[2]]) {
      stop('For one experiment, paramter var_select must be two different variable!')
    }
  }

  if (var_select[[1]] == var_select[[2]]) {
    new_x_name <- paste0(var_select[[1]],' from ',experiment_name[[1]])
    new_y_name <- paste0(var_select[[2]],' from ',experiment_name[[2]])
    data <- data %>% dplyr::rename({{new_x_name}} := paste0(var_select[[1]],'.x')) %>%
      dplyr::rename({{new_y_name}} := paste0(var_select[[2]],'.y'))
    var_select <- c(new_x_name,new_y_name)
  }

  if (is.null(palette)) {
    palette <- .get.palette.EMP(EMP)
  }else {
    palette = palette
  }
  
  # 确保方程不会被点覆盖
  ymax <- max(data[,var_select[2]])*1.1
  ymin <- min(data[,var_select[2]])*0.9
  
  data_plot <- list()
  
  if (is.null(estimate_group)) {
    
    data <- data %>% dplyr::filter(dplyr::if_all(dplyr::all_of(var_select), ~ !is.na(.)))
    
    data_plot[['pic']] <- ggplot(data=data,aes(x=!!sym(var_select[1]),y=!!sym(var_select[2]))) +
      ggpmisc::stat_poly_line(formula =formula,colour='red',se=se)+
      ggiraph::geom_jitter_interactive(aes(tooltip = paste0(primary,'\n','x: ',!!sym(var_select[1]),'\n','y: ',!!sym(var_select[2]))),
                                       position = position_jitter(height = .00000001))+
      theme_bw() +  
      xlab(paste0(var_select[1])) + ylab(paste0(var_select[2])) +
      ggpmisc::stat_poly_eq(formula = formula, 
                            aes(label = paste(after_stat(eq.label), "*\", \"*", 
                                              after_stat(rr.label), "*\", \"*", 
                                              after_stat(adj.rr.label), "*\", and \"*",
                                              after_stat(p.value.label), "*\".\"",
                                              sep = "")
                            ),size=eq_size,
                            parse = TRUE) + coord_cartesian(ylim=c(ymin,ymax)) + eval(parse(text = paste0(mytheme))) 
  }else{
    
    if(!estimate_group %in% colnames(data)){
      stop('Please check the parameter estimate_group, ',estimate_group,' is not in the data!')
    }
    
    ymax <- ymax*1.3
    
    data <- data %>% dplyr::filter(dplyr::if_all(dplyr::all_of(c(var_select,estimate_group)), ~ !is.na(.)))
    
    data_plot[['pic']]  <- ggplot(data=data,aes(x=!!sym(var_select[1]),y=!!sym(var_select[2]),colour = !!sym(estimate_group)))+
      ggpmisc::stat_poly_line(formula = formula,se=se)+
      ggiraph::geom_jitter_interactive(aes(tooltip = paste0(primary,'\n','x: ',!!sym(var_select[1]),'\n','y: ',!!sym(var_select[2]))),
                                       position = position_jitter(height = .00000001),size=3)+      theme_bw() +  
      xlab(paste0(var_select[1])) + ylab(paste0(var_select[2])) +
      ggpmisc::stat_poly_eq(formula = formula, vstep = NULL,
                            aes(label = paste(after_stat(eq.label), "*\", \"*", 
                                              after_stat(rr.label), "*\", \"*", 
                                              after_stat(adj.rr.label), "*\", and \"*",
                                              after_stat(p.value.label), "*\".\"",
                                              sep = "")
                            ),size=eq_size,
                            parse = TRUE) + coord_cartesian(ylim=c(ymin,ymax)) + 
      scale_colour_manual(values=palette,name = estimate_group) + eval(parse(text = paste0(mytheme)))
    
  }
  data_plot[['html']]  <- ggiraph::girafe(code = print(data_plot[['pic']]),width = html_width,height = html_height)
  
  .get.plot_deposit.EMP(EMP,info = 'EMP_fitline_plot') <- data_plot
  .get.plot_specific.EMP(EMP) <- show
  .get.info.EMP(EMP) <- 'EMP_fitline_plot'
  class(EMP) <- 'EMP_fitline_plot2'
  EMP
}


.show_EMP_fitplot<- function(obj,plot) {
  if(is(obj,'EMPT')){
    result <- .get.plot_deposit.EMPT(obj,info = 'EMP_fitline_plot')
  }else if (is(obj,'EMP')) {
    result <- .get.plot_deposit.EMP(obj,info = 'EMP_fitline_plot')
  }
  switch(plot,
         "pic" = print(result$pic),
         "html" = print(result$html)
  )
}


#' @param obj EMPT or EMP object
#' @param plot_category An interger.More plot style.(under constrution)
#' @param seed An interger. Set the random seed to the plot.
#' @param var_select An series of string. Select two feature to regression.
#' @param formula a formula object. Using aesthetic names x and y instead of original variable names. [y~poly(x,1,raw = FALSE)]. Detailed in ggpmisc::stat_poly_line
#' @param se An boolean. Display confidence interval around smooth. [Default:False]
#' @param eq_size An number. Equation label size. [Default:3]
#' @param estimate_group A character string. Select the colname in the coldata to compare the data in the statistical test.
#' @param show A character string include pic (default), html.
#' @param palette A series of character string. Color palette.
#' @param html_width An interger. Set the html width.
#' @param html_height An interger. Set the html height.
#' @param mytheme Modify components of a theme according to the theme.
#' @rdname EMP_fitline_plot

EMP_fitline_plot.EMPT <- function(obj,plot_category = 1,seed =123,var_select,
                               estimate_group = NULL,formula=y~poly(x,1,raw = FALSE),se=FALSE,
                               show = 'pic',palette = NULL,eq_size=3,mytheme='theme()',
                               html_width=NULL,html_height=NULL) {
  call <- match.call()
  .get.plot_category.EMPT(obj) <- plot_category
  .get.history.EMPT(obj) <- call
  switch(.get.plot_category.EMPT(obj),
         "1" = {
           withr::with_seed(seed,.EMP_fitline_plot_default.EMPT(EMPT=obj,var_select=var_select,formula=formula,se=se,
            eq_size=eq_size,show=show,estimate_group=estimate_group,palette=palette,mytheme=mytheme,html_width=html_width,html_height=html_height))
         },
         "2" = {
           # where is EMP_boxplot_assay_2?
           # withr::with_seed(seed,EMP_boxplot_assay_2(EMPT,...))
         }

  )

}


#' @param obj EMPT or EMP object
#' @param plot_category An interger.More plot style.(under constrution)
#' @param seed An interger. Set the random seed to the plot.
#' @param var_select An seris of string. Select two variable to regression.
#' @param select An series of string. The first feature of var_select will be chosen from the first experiment, and the second feature will be chosen from the second experiment.
#' @param estimate_group A character string. Select the colname in the coldata to compare the data in the statistical test.
#' @param show A character string include pic (default), html.
#' @param palette A series of character string. Color palette.
#' @param html_width An interger. Set the html width.
#' @param html_height An interger. Set the html height.
#' @param mytheme Modify components of a theme according to the theme.
#' @rdname EMP_fitline_plot

EMP_fitline_plot.EMP <- function(obj,plot_category = 1,seed =123,var_select,select=NULL,
                               estimate_group = NULL,formula=y~poly(x,1,raw = FALSE),se=FALSE,
                               show = 'pic',palette = NULL,eq_size=3,mytheme='theme()',
                               html_width=NULL,html_height=NULL) {
  #call <- match.call()
  .get.plot_category.EMP(obj) <- plot_category
  #.get.history.EMPT(EMPT) <- call
  switch(.get.plot_category.EMP(obj),
         "1" = {
           withr::with_seed(seed,.EMP_fitline_plot_default.EMP(EMP=obj,var_select=var_select,formula=formula,se=se,select=select,
            eq_size=eq_size,show=show,estimate_group=estimate_group,palette=palette,mytheme=mytheme,html_width=html_width,html_height=html_height))
         },
         "2" = {
           # where is EMP_boxplot_assay_2?
           # withr::with_seed(seed,EMP_boxplot_assay_2(EMPT,...))
         }

  )

}


