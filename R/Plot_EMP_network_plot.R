#' Network plot
#'
#' @param obj EMP object.
#' @param layout Placement of the nodes including spring and circle.
#' @param weighted A  Logical, should the analyzed network be weighted?
#' @param threshold A numeric value that defaults to 0. Edges with absolute weight that are not above this value are REMOVED from the network. 
#' @param shape A character containing the shape of the nodes. "circle", "square", "triangle" and "diamond" are supported.
#' @param node_color A vector with a color for each element in the groups list.
#' @param node_info A character in the rowdata to display more info for node.
#' @param vsize A value indicating the size of the nodes.
#' @param esize Size of the largest edge.
#' @param posCol Color of positive edges.
#' @param negCol Color of negative edges. 
#' @param fade if TRUE (default) and if 'edge.color' is assigned, transparency will be added to edges that are not transparent (or for which no transparency has been assigned) relative to the edge strength, similar if 'trans' is set to TRUE.
#' @param negDashed Logical, set to TRUE to make negative edges dashed (overwrites lty).
#' @param label.cex Scalar on the label size.
#' @param label.color Character containing the color of the labels, defaults to "black"
#' @param label.prop Controls the proportion of the width of the node that the label rescales to. Defaults to 0.9.
#' @param legend.cex Scalar of the legend. defaults to 0.4.
#' @param legend.mode Character string indicating the type of legend to be drawn. "groups" indicates the legend should be based on the groups object, "names" indicates the legend should be based on the nodeNames object, and style1 and style2 indicate the legend should be based on both. Defaults to "style1" if both "groups" and "nodeNames" arguments are used.
#' @param show A character string include net (default), node.
#' @param ... Additional parameters, see also \code{\link[bootnet]{plot.bootnetResult}} and \code{\link[qgraph]{qgraph}}.
#' @rdname EMP_network_plot
#' @export
#' @examples
#'data(MAE)
#'
#' k1 <- MAE |>
#'   EMP_assay_extract('taxonomy') |>
#'   EMP_collapse(estimate_group = 'Genus',collapse_by = 'row') |>
#'   EMP_diff_analysis(method='wilcox.test', estimate_group = 'Group') |>
#'   EMP_filter(feature_condition = pvalue<0.05)
#' 
#' k2 <- MAE |>
#'   EMP_collapse(experiment = 'untarget_metabol',na_string=c('NA','null','','-'),
#'                estimate_group = 'MS2kegg',method = 'sum',collapse_by = 'row') |>
#'   EMP_diff_analysis(method='DESeq2', .formula = ~Group) |>
#'   EMP_filter(feature_condition = pvalue<0.05)
#' 
#' (k1 + k2 ) |> 
#'   EMP_network_analysis() |> 
#'   EMP_network_plot(show = 'node',layout = 'spring',
#'                    vsize = 5,
#'                    node_info = c('Phylum','MS2class'),
#'                    legend.cex=0.3,label.cex = 1,label.prop = 0.9,
#'                    fade = T,shape = 'circle',font=2)
#' (k1 + k2 ) |> 
#'   EMP_network_analysis() |> 
#'   EMP_network_plot(show = 'node') # get the node importance
EMP_network_plot <- function(obj,layout='spring',weighted=TRUE,threshold=0,shape='circle',
                                node_color =c('#F0E442','#009E73','#E69F00','#CC79A7','#8491B4FF'),
                                node_info = NULL,
                                vsize=NULL, esize=NULL,
                                posCol='#B2182B',negCol='#002FA7',
                                fade=FALSE,negDashed = TRUE,
                                label.cex = 0.4, 
                                label.color = 'black', 
                                label.prop = 0.9, 
                                legend.cex = 0.3,
                                legend.mode = NULL,
                                show='net',...){
  deposit <- list()

  net_result <- obj@deposit$network_result
  network <- net_result[['net']]
  feature_info <- net_result[['feature_info']]
  
  if (is.null(vsize)) {
    vsize  <- 8*exp(-network$nNode/80)+1
  }
  if (is.null(esize)) {
    esize  <- 5*exp(-network$nNode/90)+1
  }
  
  
  group_info <- c()
  for (i in 1:length(feature_info)) {
    group_info <- append(group_info,(feature_info[[i]] |> dplyr::pull('.experiment_source')))
  }
  
  if (length(unique(group_info)) > 5) {
    stop("Please specify 'node_color' with at least 6 distinct colors!")
  }

  if (!is.null(node_info)) {
    node_name <- c()
    for (i in 1:length(feature_info)) {
      node_name <- append(node_name,(feature_info[[i]] |> dplyr::pull(node_info[i])))
    }
  }else{
    node_name <- NULL
  }
  
  if (is.null(legend.mode)) {
    if (is.null(node_info)) {
      legend.mode <- 'style3'
    }else{
      legend.mode <- 'style1'
    }
  }else{
    legend.mode <- legend.mode
  } 
    
    
  net_plot <- bootnet:::plot.bootnetResult(network,layout=layout,
                                      nodeNames=node_name,
                                      threshold=threshold,
                                      weighted=weighted,
                                      shape=shape,
                                      vsize=vsize,esize=esize,
                                      groups=group_info, negDashed = negDashed,
                                      color = node_color,
                                      posCol=posCol,negCol=negCol,
                                      label.cex = label.cex, # scalar on label size
                                      label.color = label.color, # string on label colors
                                      label.prop = label.prop, # proportion of the width of the node that the label scales
                                      legend.cex =legend.cex,
                                      legend.mode = legend.mode,fade=fade,
                                      DoNotPlot=TRUE,...) 

  centrality_plot <- qgraph::centralityPlot(network,include='all',print=FALSE) |>
               spsUtil::quiet(print_cat = TRUE, message = TRUE, warning = TRUE)

  deposit[['net']] <- net_plot
  deposit[['centrality']] <- centrality_plot

  .get.plot_deposit.EMP(obj,info = 'EMP_network_plot') <- deposit
  .get.plot_specific.EMP(obj) <- show
  .get.info.EMP(obj) <- 'EMP_network_plot'
  class(obj) <- 'EMP_network_plot'
  return(obj)
}



.show_EMP_network_plot <- function(obj,plot,...) {
  result <- .get.plot_deposit.EMP(obj,info = 'EMP_network_plot')
  switch(plot,
         "net" = plot(result[['net']]),
         "node" = print(result[['centrality']]) |>
               spsUtil::quiet(print_cat = TRUE, message = TRUE, warning = TRUE)
  )
}





