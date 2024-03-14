
#' Title
#'
#' @param object wait_for_add
#'
#' @return xx object
##' @importFrom methods show
##' @exportMethod show
#' @examples
#' # add example
setMethod("show", "EMPT",
          function(object) {
            info <- .get.info.EMPT(object)
            switch(info,
                   "EMP_decostand" = {
                     .get.result.EMPT(object) %>% print()
                   },
                   "EMP_rrarefy" = {
                     .get.result.EMPT(object) %>% print()
                   },
                   "EMP_diff_analysis" = {
                     .get.result.EMPT(object) %>% print()
                   },
                   "EMP_enrich_analysis" = {
                     object@deposit[["enrich_data"]] %>% tibble::as_tibble() %>% print(n=Inf)
                   },
                   "EMP_assay_data" = {
                     .get.result.EMPT(object) %>% print()
                   },
                   "EMP_cluster_analysis" = {
                     .get.result.EMPT(object) %>% print()
                   },
                   "EMP_WGCNA_cluster_analysis" = {
                     deposit <- .get.result.EMPT(object)
                     net <- deposit[['WGCNA_cluster_result']]
                     Module <- WGCNA::labels2colors(net$colors)
                     table(Module) %>% print()
                     deposit[['WGCNA_cluster_df']] %>% print()
                   },
                   "EMP_WGCNA_cor_analysis" = {
                      object <- .get.result.EMPT(object)
                      cat('EMP_WGCNA_cor_analysis:','\n')
                      cat('Cor-relationship matrix:',dim(object$correlation)[1],'x',dim(object$correlation)[2],'\n')
                      cat(object$cor_info[1],'observation:',object$n.obs[1],'\n')
                      cat(object$cor_info[2],'observation:',object$n.obs[2],'\n')
                      cat('Intersect observation:',object$n.obs[3])
                   },
                    "EMP_WGCNA_cor_heatmap" = {
                      .get.result.EMPT(object) %>% print()
                   },
                   "EMP_alpha_analysis" = {
                     object@deposit$diversity_result %>% print()
                   },
                   "EMP_beta_analysis" = {
                     object@deposit$distance_result %>% print()
                   },
                   "EMP_dimension_analysis" = {
                     .get.result.EMPT(object) %>% print()
                   },
                    "EMP_dimension_analysis_scatterplot" = {
                     .show_EMP_dimension_analysis_scatterplot(object,.get.plot_specific.EMPT(object))
                   },
                   "EMP_beta_analysis_boxplot" = {
                     .show_EMP_beta_boxplot(object,.get.plot_specific.EMPT(object))
                   },
                   "EMP_alpha_analysis_boxplot" = {
                     .show_EMP_alpha_boxplot(object,.get.plot_specific.EMPT(object))
                   },
                   "EMP_assay_boxplot" = {
                     .show_EMP_assay_boxplot(object,.get.plot_specific.EMPT(object))
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
                   "EMP_decostand" = {
                     for (str in object@message_info) {
                      message_wrap(str)
                     }
                     return(.get.assay.EMPT(object))
                   },
                   "EMP_rrarefy" = {
                     for (str in object@message_info) {
                      message_wrap(str)
                     }
                     return(.get.assay.EMPT(object))
                   },
                   "EMP_diff_analysis" = {
                     for (str in object@message_info) {
                      message_wrap(str)
                     }
                     return(object@deposit[["diff_analysis_result"]])
                   },
                   "EMP_dimension_analysis" = {
                     for (str in object@message_info) {
                      message_wrap(str)
                     }
                     deposit <- list()
                     deposit[['dimension_coordinate']] <- object@deposit$dimension_coordinate
                     deposit[['dimension_VIP']] <- object@deposit$dimension_VIP
                     deposit[['dimension_axis']] <- object@deposit$dimension_axis
                     return(deposit)
                   },
                   "EMP_assay_data" = {
                     for (str in object@message_info) {
                      message_wrap(str)
                     }
                     return(.get.assay.EMPT(object))
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
                  {
                   print('No info is matched!')
                  }
  )
}


#' Title
#'
#' @param object wait_for_add
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
setMethod("show", "EMP",
          function(object) {
            info <- .get.info.EMP(object)
            switch(info,
                   "EMP_list_data" = {
                     experiment_names <- names(object@ExperimentList)
                     experiment_num <- length(object@ExperimentList)
                     cat('EMP object contains',experiment_num, 'experiment list:','\n')
                     for (i in names(object@ExperimentList)) {
                       cat(i,'\n')
                     }
                   },
                   "EMP_cor_analysis" = {
                     object <- object@deposit$cor_analysis_result
                     cat('EMP_cor_analysis:','\n')
                     cat('Cor-relationship matrix:',dim(object$correlation)[1],'x',dim(object$correlation)[2],'\n')
                     cat(object$cor_info[1],'observation:',object$n.obs[1],'\n')
                     cat(object$cor_info[2],'observation:',object$n.obs[2],'\n')
                     cat('Intersect observation:',object$n.obs[3])
                   },
                   "EMP_WGCNA_cor_analysis2" = {
                     object <- object@deposit$WGCNA_cor_analysis_result
                     cat('EMP_WGCNA_cor_analysis:','\n')
                     cat('Cor-relationship matrix:',dim(object$correlation)[1],'x',dim(object$correlation)[2],'\n')
                     cat(object$cor_info[1],'observation:',object$n.obs[1],'\n')
                     cat(object$cor_info[2],'observation:',object$n.obs[2],'\n')
                     cat('Intersect observation:',object$n.obs[3])
                   },
                    "EMP_WGCNA_cor_heatmap2" = {
                      .get.result.EMP(object) %>% print()
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
                      message_wrap(str)
                     }
                     return(object@deposit[["cor_analysis_result"]])
                   },
                  "EMP_WGCNA_cor_analysis2" = {
                     for (str in object@message_info) {
                      message_wrap(str)
                     }
                     return(object@deposit[["WGCNA_cor_analysis_result"]])
                   },
                  "EMP_WGCNA_cor_heatmap2" = {
                     return(.get.plot_deposit.EMP(object,info='EMP_WGCNA_cor_heatmap')) ## here not EMP_WGCNA_cor_heatmap2
                   }

  )
}




