CorRcpp <- function(x=x,y=NULL,type=c("pearson","spearman")) {
    stopifnot(is.data.frame(x))
	  type <- match.arg(type)
	  names_x <- names(x)
  if(is.null(y)){
    if(type=="spearman"){
      x <- apply(x, 2, rank)
    }
    corres <-  cp_cor_s(mat = as.matrix(x))
    corres$R_matrix <- as.data.frame(corres$R_matrix)
    corres$P_matrix <- as.data.frame(corres$P_matrix)
    row.names(corres$R_matrix) <- names_x
    row.names(corres$P_matrix) <- names_x
    names(corres$R_matrix) <- names_x
    names(corres$P_matrix) <- names_x
  }else{
    names_y <- names(y)
    if(type=="spearman"){
	stopifnot(is.data.frame(y))
      x <- apply(x, 2, rank)
      y <- apply(y, 2, rank)
    }
    corres <-  cp_cor_t(mat = as.matrix(x),mat2=as.matrix(y))
    corres$R_matrix <- as.data.frame(corres$R_matrix)
    corres$P_matrix <- as.data.frame(corres$P_matrix)
    row.names(corres$R_matrix) <- names_x
    row.names(corres$P_matrix) <- names_x
    names(corres$R_matrix) <- names_y
    names(corres$P_matrix) <- names_y
  }
  return(corres)
}
