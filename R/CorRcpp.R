#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom Rcpp evalCpp
#' @useDynLib EasyMultiProfiler
NULL

CorRcpp <- function(x=x,y=NULL,type=c("pearson","spearman")) {
    stopifnot(is.data.frame(x))
	  type <- match.arg(type)
  if(is.null(y)){
    if(type=="spearman"){
      x <- apply(x, 2, rank)
    }
    corres <-  cp_cor_s(mat = as.matrix(x))
    corres$R_matrix <- as.data.frame(corres$R_matrix)
    corres$P_matrix <- as.data.frame(corres$P_matrix)
    row.names(corres$R_matrix) <- names(x)
    row.names(corres$P_matrix) <- names(x)
    names(corres$R_matrix) <- names(x)
    names(corres$P_matrix) <- names(x)
  }else{
    if(type=="spearman"){
	stopifnot(is.data.frame(y))
      x <- apply(x, 2, rank)
      y <- apply(y, 2, rank)
    }
    corres <-  cp_cor_t(mat = as.matrix(x),mat2=as.matrix(y))
    corres$R_matrix <- as.data.frame(corres$R_matrix)
    corres$P_matrix <- as.data.frame(corres$P_matrix)
    row.names(corres$R_matrix) <- names(x)
    row.names(corres$P_matrix) <- names(x)
    names(corres$R_matrix) <- names(y)
    names(corres$P_matrix) <- names(y)
  }
  return(corres)
}
