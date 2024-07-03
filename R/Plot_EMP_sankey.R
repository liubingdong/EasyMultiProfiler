#' @param positive_col Edge colour for positive relation. [Default:#CC79A7]
#' @param negtive_col  Edge colour for negitive relation. [Default:steelblue]
#' @param palette Colour palette for nodes.
#' @rdname EMP_sankey_plot
#' @importFrom networkD3 sankeyNetwork
EMP_sankey_plot.EMP_cor_analysis <- function(obj,positive_col = '#CC79A7',negtive_col = 'steelblue',palette=c("#009E73","#F0E442","#6c4f9d","#CB563f","#4DBBD5FF","#00A087FF","#b66e1e","#8491B4FF")){

  name <- NULL
  call <- match.call()
  if (inherits(obj,"EMP")) {
    EMP <- obj
  }else{
    stop('Please check the input data!')
  }  

  result <- .get.result.EMP(EMP,info = 'EMP_cor_analysis')

  data_long <- result$cor_p_filter
  feature_info <- result$feature_info

  #set the edge
  idx <- data_long$value<0
  data_long$group[idx] <- 'negetive'
  data_long$group[!idx] <- 'positive'
  data_long$value=abs(data_long$value)*100

  data_long$group=as.factor(data_long$group)
  

  colnames(data_long)=c('source','target','value','pvalue','group')

  # From these flows we need to create a node data frame: it lists every entities involved in the flow
  nodes <-  c(as.character(data_long$source), as.character(data_long$target)) %>% unique()
  nodes <- feature_info %>% dplyr::filter(name %in% {{nodes}})
  group_name <- nodes$group %>% unique()
  if(length(group_name) > 8){
    warning("Group number has achieved the max number of palette, please set more palette!")
  }  

  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  data_long$IDsource=match(data_long$source, nodes$name)-1 
  data_long$IDtarget=match(data_long$target, nodes$name)-1

  
  edge_col <- paste0(paste0('\"',c(negtive_col,positive_col),'\"'),collapse=',')
  d3_node_name <- paste0(paste0('\"',group_name,'\"'),collapse=',')
  d3_node_col <- paste0(paste0('\"',palette[1:length(group_name)],'\"'),collapse=',')


  ## origin method
  ## my_color <- 'd3.scaleOrdinal() .domain(["negitive", "positive","my_unique_group1","my_unique_group2","my_unique_group3","my_unique_group4"]) .range([ "steelblue", "darkred","#00A087FF","#CC79A7","#3C5488FF","#F0E442"])'
  ## use paste0 to copy this
  my_color <- paste0('d3.scaleOrdinal() .domain(["negetive", "positive",',d3_node_name,']) .range([ ',edge_col,',',d3_node_col,'])')




  p1 <- networkD3::sankeyNetwork(Links = data_long, Nodes = nodes, Source = "IDsource", Target = "IDtarget",
                                Value = "value", NodeID = "name",sinksRight=F,LinkGroup="group",colourScale=my_color,NodeGroup="group")


  .get.plot_deposit.EMP(EMP,info='EMP_cor_sankey') <- p1
  .get.info.EMP(EMP) <- 'EMP_cor_sankey'
  .get.history.EMP(EMP) <- call
  class(EMP) <- 'EMP_cor_sankey'
  return(EMP)
}
