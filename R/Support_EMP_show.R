#' @importFrom methods show    
setMethod("show", "MultiAssayExperiment", function(object) {
  o_class <- class(object)
  o_len <- length(object)
  o_names <- names(object)
  if (!length(o_names)) {
    o_names <- "none"
  }
  c_elist <- class(MultiAssayExperiment::experiments(object))
  c_mp <- class(MultiAssayExperiment::colData(object))
  c_sm <- class(MultiAssayExperiment::sampleMap(object))
  cat(sprintf("A %s", o_class),
      "object of", o_len, "listed\n",
      ifelse(o_len == 1L, "experiment", "experiments"),
      "with",
      ifelse(identical(o_names, "none"), "no user-defined names",
             ifelse(length(o_names) == 1L, "a user-defined name",
                    "user-defined names")),
      ifelse(length(o_len) == 0L, "or", "and"),
      ifelse(length(o_len) == 0L, "classes.",
             ifelse(o_len == 1L,
                    "respective class.\n", "respective classes.\n")),
      "Containing an ")
  show(experiments(object))
  cat("Functionality:\n ","For more information, please type help(EasyMultiProfiler)"
)
})


#' @importFrom methods show 
setMethod("show", "EMPT",
          function(object) {
            info <- .get.info.EMPT(object)
            switch(info,
                   "EMP_assay_data" = {
                     for (str in object@message_info) {
                      EMP_message(str,color = 32,order = 1,show='message')
                     }
                     object %>% enhance_print(n=10)
                   },   
                   "EMP_mutate_assay" = {
                     print(.get.assay.EMPT(object))
                   }, 
                   "EMP_mutate_row" = {
                     print(.get.row_info.EMPT(object))
                   },
                   "EMP_mutate_col" = {
                     print(.get.mapping.EMPT(object))
                   },                              
                   "EMP_decostand" = {
                     for (str in object@message_info) {
                      EMP_message(str,color = 32,order = 1,show='message')
                     }
                     object %>% enhance_print(n=10)
                   },
                   "EMP_rrarefy" = {
                     for (str in object@message_info) {
                      EMP_message(str,color = 32,order = 1,show='message')
                     }
                     object %>% enhance_print(n=10)
                   },
                   "EMP_diff_analysis" = {
                     for (str in object@message_info) {
                      EMP_message(str,color = 32,order = 1,show='message')
                     }                    
                     object %>% enhance_print(n=10)
                   },
                   "EMP_enrich_analysis" = {
                    try(object@deposit[["enrich_data"]]@compareClusterResult %>% tibble::as_tibble() %>% print(n=Inf),silent=TRUE)  
                    try(object@deposit[["enrich_data"]]@result %>% tibble::as_tibble() %>% print(n=Inf),silent=TRUE)                   
                   },
                   "EMP_cluster_analysis" = {
                     cluster_re <- .get.result.EMPT(object)
                     for (i in names(cluster_re)) {
                       cluster_re[[i]]['cluster'] %>% table() %>% print()
                       cluster_re[i] %>% print()
                     }
                   },
                   "EMP_WGCNA_cluster_analysis" = {
                     deposit <- .get.result.EMPT(object)
                     net <- deposit[['WGCNA_cluster_result']]
                     Module <- WGCNA::labels2colors(net$colors)
                     table(Module) %>% print()
                     deposit[['WGCNA_cluster_df']] %>% print()
                     mergedColors = WGCNA::labels2colors(net$colors)
                     WGCNA::plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                                        dendroLabels = FALSE, hang = 0.03,
                                        addGuide = TRUE, guideHang = 0.05,
                                        abHeight=net$mergeCutHeight)

                   },
                   "EMP_WGCNA_cor_analysis" = {
                      object <- .get.result.EMPT(object,info = 'EMP_WGCNA_cor_analysis')
                      str <- paste0('EMP_WGCNA_cor_analysis:','\n',
                                     'Cor-relationship matrix:',dim(object$correlation)[1],'x',dim(object$correlation)[2],'\n',
                                     object$cor_info[1],'observation:',object$n.obs[1],'\n',
                                     object$cor_info[2],'observation:',object$n.obs[2],'\n',
                                     'Intersect observation:',object$n.obs[3])
                      EMP_message(str,color = 32,order = 1,show='message')
                   },
                    "EMP_WGCNA_cor_heatmap" = {
                      .get.result.EMPT(object) %>% print()
                   },
                   "EMP_alpha_analysis" = {
                     for (str in object@message_info) {
                      EMP_message(str,color = 32,order = 1,show='message')
                     }
                     object %>% enhance_print(n=10)
                   },
                   "EMP_marker_analysis" = {
                     .get.result.EMPT(object) %>% print()                  
                   },
                   "EMP_dimension_analysis" = {
                     .get.result.EMPT(object) %>% print()
                   },
                   "EMP_network_analysis" = {
                    .network_print(object)
                    deposit <- .get.result.EMPT(object) |> suppressMessages()
                    deposit[['net_centrality']] %>% print()
                   },
                    "EMP_dimension_analysis_scatterplot" = {
                      set.seed(123)
                     .show_EMP_dimension_analysis_scatterplot(object,.get.plot_specific.EMPT(object))
                   },
                   "EMP_alpha_analysis_boxplot" = {
                     set.seed(123)
                     .show_EMP_alpha_boxplot(object,.get.plot_specific.EMPT(object))
                   },
                   "EMP_assay_boxplot" = {
                     set.seed(123)
                     .show_EMP_assay_boxplot(object,.get.plot_specific.EMPT(object))
                   },
                   "EMP_structure_plot" = {
                     .show_EMP_structure_plot(object,.get.plot_specific.EMPT(object))
                   },
                   "EMP_enrich_analysis_dotplot" = {
                     .show_EMP_dotplot_enrich(object,.get.plot_specific.EMPT(object))
                   },
                   "EMP_enrich_analysis_netplot" = {
                     .show_EMP_netplot_enrich(object,.get.plot_specific.EMPT(object))
                   },
                   "EMP_enrich_analysis_curveplot" = {
                     .show_EMP_curveplot_enrich(object,.get.plot_specific.EMPT(object))
                   },
                   "EMP_diff_volcanol_plot" = {
                     .show_EMP_diff_volcanol_plot(object,.get.plot_specific.EMPT(object))
                   },
                   "EMP_assay_heatmap" = {
                     .get.result.EMPT(object) %>% print()
                   },
                   "EMP_fitline_plot" = {
                     .show_EMP_fitplot(object,.get.plot_specific.EMPT(object)) %>% suppressWarnings()
                   },
                   "EMP_network_plot" = {
                    show <-.get.plot_specific.EMPT(object)
                    return(.show_EMP_network_plot(object,plot=show))                   
                  }
                   )
          }
)

.get.result.EMPT <- function(object,info=NULL) {
  if (is.null(info)) {
    info <- .get.info.EMPT(object)
  }else{
    info <- info
  }

  switch(info,
                   "EMP_assay_data" = {
                     for (str in object@message_info) {
                      EMP_message(str,color = 32,order = 1,show='message')
                     }
                     return(.get.assay.EMPT(object))
                   },
                   "EMP_mutate_assay" = {
                     return(.get.assay.EMPT(object))
                   }, 
                   "EMP_mutate_row" = {
                     return(.get.row_info.EMPT(object))
                   },
                   "EMP_mutate_col" = {
                     return(.get.mapping.EMPT(object))
                   },                                                            
                   "EMP_decostand" = {
                     for (str in object@message_info) {
                      EMP_message(str,color = 32,order = 1,show='message')
                     }
                     return(.get.assay.EMPT(object))
                   },
                   "EMP_rrarefy" = {
                     for (str in object@message_info) {
                      EMP_message(str,color = 32,order = 1,show='message')
                     }
                     return(.get.assay.EMPT(object))
                   },
                   "EMP_diff_analysis" = {
                     for (str in object@message_info) {
                      EMP_message(str,color = 32,order = 1,show='message')
                     }
                     return(object@deposit[["diff_analysis_result"]])
                   },
                   "EMP_dimension_analysis" = {
                     for (str in object@message_info) {
                      EMP_message(str,color = 32,order = 1,show='message')
                     }
                     deposit <- list()
                     deposit[['dimension_coordinate']] <- object@deposit$dimension_coordinate
                     deposit[['dimension_VIP']] <- object@deposit$dimension_VIP
                     deposit[['dimension_axis']] <- object@deposit$dimension_axis
                     return(deposit)
                   },
                    "EMP_network_analysis" = {
                     for (str in object@message_info) {
                      EMP_message(str,color = 32,order = 1,show='message')
                     }
                     deposit <- list()
                     deposit[['net']] <- object@deposit[['net']] 
                     deposit[['net_feature_info']] <- object@deposit[['net_feature_info']] 
                     deposit[['net_centrality']] <- object@deposit[['net_centrality']] 
                     return(deposit)
                   },   
                   "EMP_cluster_analysis" = {
                     check_sample_cluster <- is.null(object@deposit[['sample_cluster_result']])
                     check_feature_cluster <- is.null(object@deposit[['feature_cluster_result']])
                     deposit_list <- list()
                     if(!check_sample_cluster){
                      deposit_list[['sample_cluster_result']] <- object@deposit[['sample_cluster_result']]
                     }
                     if(!check_feature_cluster){
                      deposit_list[['feature_cluster_result']] <- object@deposit[['feature_cluster_result']]
                     }
                     return(deposit_list)
                   },
                   "EMP_WGCNA_cluster_analysis" = {
                     deposit <- list()
                     deposit[['WGCNA_cluster_result']] <- object@deposit_append[['feature_WGCNA_cluster_result']]
                     deposit[['WGCNA_cluster_df']] <- object@deposit[['feature_WGCNA_cluster_result']]
                     return(deposit)
                   },
                   "EMP_enrich_analysis" = {
                     return(object@deposit[["enrich_data"]])
                   },
                   "EMP_alpha_analysis" = {
                    return(object@deposit[["diversity_result"]])
                   },
                   "EMP_WGCNA_cor_analysis" = {
                    return(object@deposit_append[["WGCNA_cor_result"]])
                   },
                   "EMP_WGCNA_cor_heatmap" = {
                    return(.get.plot_deposit.EMPT(object,info='EMP_WGCNA_cor_heatmap'))
                   },
                    "EMP_marker_analysis" = {
                            method <- .get.method.EMPT(object)
                            deposit <- list()
                            if(method == 'boruta') {
                              deposit <- list(Boruta_model=object@deposit[['Boruta_model']],
                                              Boruta_feature_importance=object@deposit[['Boruta_feature_importance']])
                              return(deposit)
                              }else if(method == 'randomforest'){
                                deposit <- list(rf_model=object@deposit[['rf_model']],
                                              rf_feature_importance=object@deposit[['rf_feature_importance']])
                                return(deposit)                                
                              }else if(method == 'xgboost'){
                                deposit <- list(xgb_model=object@deposit[['xgb_model']],
                                              xgb_feature_importance=object@deposit[['xgb_feature_importance']])
                                return(deposit) 
                              }else if(method == 'lasso'){
                                deposit <- list(lasso_model=object@deposit[['lasso_model']],
                                              lasso_feature_importance=object@deposit[['lasso_feature_importance']])
                                return(deposit)                                 
                              }else{
                                return(NULL)
                              }
                   },
                   "EMP_assay_heatmap" = {
                      return(.get.plot_deposit.EMPT(object,info='EMP_assay_heatmap'))
                   },                
                  {
                   print('No info is matched!')
                  }
  )
}


setMethod("show", "EMP",
          function(object) {
            info <- .get.info.EMP(object)
            switch(info,
                   "EMP_list_data" = {
                     experiment_names <- names(object@ExperimentList)
                     experiment_num <- length(object@ExperimentList)
                     str1 <- paste0('EMP object contains',experiment_num, 'experiment list:','\n')
                     str2 <- paste(names(object@ExperimentList),collapse ='\n')
                     str3 <- paste0(str1,str2)
                     EMP_message(str3,color = 32,order = 1,show='message')
                   },
                   "EMP_cor_analysis" = {
                      cor_method <- .get.method.EMP(object)
                      object <- .get.result.EMP(object,info = 'EMP_cor_analysis')
                      str <- paste0('EMP_cor_analysis:','\n',
                             'Cor-relationship observation:',paste0(object$n.obs,collapse = ' x '),'\n',
                             'Cor-relationship method: ',cor_method)
                      EMP_message(str,color = 32,order = 1,show='message')
                   },
                   "EMP_network_analysis2" = {
                    .network_print(object)
                    deposit <- .get.result.EMP(object) |> suppressMessages()
                    deposit[['net_centrality']] %>% print()
                   },
                   "EMP_network_plot2" = {
                    show <-.get.plot_specific.EMP(object)
                    return(.show_EMP_network_plot(object,plot=show))
                   },
                   "EMP_WGCNA_cor_analysis2" = {
                     object <- object@deposit$WGCNA_cor_analysis_result
                     str <- paste0('EMP_WGCNA_cor_analysis:','\n',
                                    'Cor-relationship matrix:',dim(object$correlation)[1],'x',dim(object$correlation)[2],'\n',
                                    object$cor_info[1],'observation:',object$n.obs[1],'\n',
                                    object$cor_info[2],'observation:',object$n.obs[2],'\n',
                                    'Intersect observation:',object$n.obs[3])
                     EMP_message(str,color = 32,order = 1,show='message')
                   },
                    "EMP_WGCNA_cor_heatmap2" = {
                      .get.result.EMP(object) %>% print()
                   },
                    "EMP_cor_heatmap" = {
                      .get.result.EMP(object) %>% print()
                   },
                    "EMP_cor_sankey" = {
                      .get.result.EMP(object) %>% print()
                   },
                   "EMP_multi_same_df" = {
                      .get.result.EMP(object) %>% print()
                   },
                   "EMP_multi_diff_enrich" = {
                      .get.result.EMP(object) %>% tibble::as_tibble() %>% print(n=Inf)
                   },
                   "EMP_multi_same_enrich" = {
                      .get.result.EMP(object) %>% tibble::as_tibble() %>% print(n=Inf)
                   },
                   "EMP_multi_diff_enrich_dotplot" = {
                      .get.result.EMP(object) %>% print()
                   },
                   "EMP_multi_same_enrich_dotplot" = {
                      .get.result.EMP(object) %>% print()
                   },
                   "EMP_multi_diff_enrich_netplot" = {
                      .get.result.EMP(object) %>% print()
                   },
                   "EMP_multi_same_enrich_netplot" = {
                      .get.result.EMP(object) %>% print()
                   },
                   "EMP_fitline_plot" = {
                     .show_EMP_fitplot(object,.get.plot_specific.EMP(object)) %>% suppressWarnings()
                   },
                   {
                     print('No info is matched!')
                   }
            )
          }
)


.get.result.EMP <- function(object,info=NULL) {
  if (is.null(info)) {
    info <- .get.info.EMP(object)
  }else{
    info <- info
  }

  switch(info,
                  "EMP_cor_analysis" = {
                     for (str in object@message_info) {
                      EMP_message(str,color = 32,order = 1,show='message')
                     }
                     return(object@deposit[["cor_analysis_result"]])
                   },
                  "EMP_cor_heatmap" = {
                     for (str in object@message_info) {
                      EMP_message(str,color = 32,order = 1,show='message')
                     }
                     return(.get.plot_deposit.EMP(object,info='EMP_cor_heatmap'))
                   },
                  "EMP_cor_sankey" = {
                     for (str in object@message_info) {
                      EMP_message(str,color = 32,order = 1,show='message')
                     }
                     return(.get.plot_deposit.EMP(object,info='EMP_cor_sankey'))
                   },
                  "EMP_network_analysis2" = {
                     for (str in object@message_info) {
                      EMP_message(str,color = 32,order = 1,show='message')
                     }
                     return(object@deposit[["EMP_network_analysis"]])
                   },
                  "EMP_WGCNA_cor_analysis2" = {
                     for (str in object@message_info) {
                      EMP_message(str,color = 32,order = 1,show='message')
                     }
                     return(object@deposit[["WGCNA_cor_analysis_result"]])
                   },
                  "EMP_WGCNA_cor_heatmap2" = {
                     return(.get.plot_deposit.EMP(object,info='EMP_WGCNA_cor_heatmap')) ## here not EMP_WGCNA_cor_heatmap2
                   },
                  "EMP_multi_same_df" = {
                     return(object@deposit[['multi_same_df']]) 
                   },
                   "EMP_multi_diff_enrich" = {
                     try(return(object@deposit[['multi_diff_enrich']]@compareClusterResult),silent=TRUE)
                     try(return(object@deposit[['multi_diff_enrich']]@result),silent=TRUE)
                   },
                   "EMP_multi_same_enrich" = {
                     try(return(object@deposit[['multi_same_enrich']]@compareClusterResult),silent=TRUE)
                     try(return(object@deposit[['multi_same_enrich']]@result),silent=TRUE)                   
                   },
                  "EMP_multi_diff_enrich_dotplot" = {
                     .show_EMP_dotplot_enrich(object,.get.plot_specific.EMP(object))
                   },
                  "EMP_multi_same_enrich_dotplot" = {
                     .show_EMP_dotplot_enrich(object,.get.plot_specific.EMP(object))
                   },
                  "EMP_multi_diff_enrich_netplot" = {
                     .show_EMP_netplot_enrich(object,.get.plot_specific.EMP(object))
                   },
                  "EMP_multi_same_enrich_netplot" = {
                     .show_EMP_netplot_enrich(object,.get.plot_specific.EMP(object))
                   }  

  )
}



