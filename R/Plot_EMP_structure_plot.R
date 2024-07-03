
#' @importFrom forcats fct_relevel
EMP_structure_plot_default <-function(EMPT,method = 'mean',top_num = 10,estimate_group = NULL,show = 'pic',palette = NULL,
                                      ncol = NULL,mytheme = 'theme()'){

  primary <- feature <- value <- assay_data <- cols <- mapping <- str_data <- NULL
  
  assay_data <- .get.assay.EMPT(EMPT)

  if (is.null(palette)) {
    col_values <- .get.palette.EMPT(EMPT)
  }else {
    col_values <- palette
  }

  if (!is.null(estimate_group)) {
    mapping <- .get.mapping.EMPT(EMPT) %>% dplyr::select(primary,!!estimate_group)
    str_data <- .get.assay.EMPT(EMPT) %>% dplyr::left_join(mapping,by ='primary')
  }else{
    str_data <- .get.assay.EMPT(EMPT)
  }

  top_str_data <- top_exper_caculate(str_data,top_num = top_num,estimate_group = estimate_group) %>%
    tidyr::pivot_longer(
      cols = !dplyr::any_of(c('primary',!!estimate_group)),
      names_to = 'feature',
      values_to = "value"
    ) 

  ## In case that data contains negative value of feature
  if (any(top_str_data$value <0)) {
    stop("EMP_structure_plot dont allow negative value, please check the assay data!")
  }

  ## In case that data only contain one feature
  if (dplyr::n_distinct(top_str_data$feature) >=2 & 'Others' %in% top_str_data$feature) {
    top_str_data <- top_str_data %>% dplyr::mutate(feature = forcats::fct_relevel(factor(feature), "Others", after = Inf))
  }
  
  data_plot <- list()
  if (!is.null(estimate_group)) {
    data_plot[['pic']]  <- ggplot(top_str_data, aes(x=primary,y=value, fill = feature))+  geom_col(position = 'stack', width = 0.8) +
      scale_fill_manual(values =col_values) +
      facet_wrap({{estimate_group}}, scales = 'free_x', ncol = ncol,strip.position = 'top') +
      theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
      theme(axis.text = element_text(size = 5), axis.title = element_text(size = 13), legend.title = element_blank(), legend.text = element_text(size = 11)) +
      xlab(NULL) + theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
      eval(parse(text = paste0(mytheme)))
  }else{
    data_plot[['pic']]  <- ggplot(top_str_data, aes(x=primary,y=value, fill = feature))+  geom_col(position = 'stack', width = 0.8) +
      scale_fill_manual(values =col_values) +
      theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
      theme(axis.text = element_text(size = 5), axis.title = element_text(size = 13), legend.title = element_blank(), legend.text = element_text(size = 11)) +
      xlab(NULL) + theme(axis.text.x =element_text(angle = 45, hjust = 1,size = 10)) +
      eval(parse(text = paste0(mytheme)))
  }
  
  .get.plot_deposit.EMPT(EMPT,info = 'EMP_structure_plot') <- data_plot
  .get.plot_specific.EMPT(EMPT) <- show
  .get.estimate_group.EMPT(EMPT) <- estimate_group
  EMPT@algorithm <- 'EMP_structure_plot'
  .get.info.EMPT(EMPT) <- 'EMP_structure_plot'
  class(EMPT) <- 'EMP_structure_plot'
  EMPT
  
}

.show_EMP_structure_plot <- function(obj,plot) {
  result <- .get.plot_deposit.EMPT(obj,info = 'EMP_structure_plot')
  switch(plot,
         "pic" = print(result$pic),
         "html" = print(result$html)
  )
}


top_exper_caculate <- function(assay_data,estimate_group=NULL,structure_method = 'mean',top_num = 10){
  
  primary <- NULL

  data_abundace <- assay_data %>% dplyr::select(-primary,-dplyr::any_of(!!estimate_group))
  
  
  feature_num <- ncol(data_abundace)
  if(top_num >10){
    stop("Top number could not exceed 10.")
  }
  
  if (top_num < feature_num & 0 < top_num) {
    if (structure_method=='median') {
      estimate_index <- apply(data_abundace,2,median)
    }else if (structure_method=='mean') {
      estimate_index <- apply(data_abundace,2,mean)
    }else if (structure_method=='max'){
      estimate_index <- apply(data_abundace,2,max)
    }else if (structure_method=='min'){
      estimate_index <- apply(data_abundace,2,min)
    }else{
      warning('structure_method is wrong, please check !')
    }
    id <- names(sort(estimate_index,decreasing = T)[1:top_num])
    data_select <- assay_data %>%
      dplyr::select(primary,dplyr::all_of(!!estimate_group),id)
    data_others <- data_abundace %>% 
      dplyr::select(-{{id}}) %>% 
      rowSums()
    data_select_combie <- data_select %>% dplyr::mutate(Others=data_others)
  }else{
    data_select_combie <- assay_data
  }
  return(data_select_combie)
} 


#' @param obj EMPT object
#' @param plot_category An interger.More plot style.(under constrution)
#' @param method A character string including mean,median,max and min. The name of the statistical test that is applied to select the top feature.
#' @param top_num An interger. Max and default number is 10.
#' @param estimate_group A character string. Select the colname in the coldata to divide the data in the plot.
#' @param ncol An interger. Set the col number in the facet plot.
#' @param show A character string include pic (default), html(under constrution).
#' @param palette A series of character string. Color palette.
#' @param mytheme Modify components of a theme according to the ggplot2::theme.
#' @rdname EMP_structure_plot


EMP_structure_plot.EMP_assay_data <- function(obj,plot_category = 1,method = 'mean',top_num=10,
                               estimate_group = NULL,ncol = NULL,show = 'pic',palette = NULL,
                               mytheme = 'theme()') {
  call <- match.call()
  .get.plot_category.EMPT(obj) <- plot_category
  .get.history.EMPT(obj) <- call
  switch(.get.plot_category.EMPT(obj),
         "1" = {
           EMP_structure_plot_default(EMPT=obj,method = method,top_num = top_num,
            estimate_group=estimate_group,show=show,palette=palette,ncol=ncol,mytheme=mytheme)
         },
         "2" = {
           # where is EMP_boxplot_assay_2?
           # withr::with_seed(seed,EMP_boxplot_assay_2(EMPT,...))
         }

  )

}
#' @rdname EMP_structure_plot
EMP_structure_plot.EMP_decostand <- function(obj,plot_category = 1,method = 'mean',top_num = 10,
                               estimate_group = NULL,ncol = NULL,show = 'pic',palette = NULL,
                               mytheme = 'theme()') {
  call <- match.call()
  .get.plot_category.EMPT(obj) <- plot_category
  .get.history.EMPT(obj) <- call
  switch(.get.plot_category.EMPT(obj),
         "1" = {
           EMP_structure_plot_default(EMPT=obj,method = method,top_num = top_num,
            estimate_group=estimate_group,show=show,palette=palette,ncol=ncol,mytheme=mytheme)
         },
         "2" = {
           # where is EMP_boxplot_assay_2?
           # withr::with_seed(seed,EMP_boxplot_assay_2(EMPT,...))
         }

  )

}
