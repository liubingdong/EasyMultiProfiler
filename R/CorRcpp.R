#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom Rcpp evalCpp
#' @useDynLib EasyMultiProfiler
NULL
#' Fast correlation analysis
#' @param x Dataframe.
#' @param y Dataframe.
#' @param type  A character string. Methods include pearson (default), spearman.
#' @rdname CorRcpp
#' @return list
#' @export
#'
#' @examples
#' data_df1 <- data.frame(
#'   A = c(11, 23, 3, 45, 5),
#'   B = c(62, 7, 82, 9, 10),
#'   C = c(11, 54, 13, 42, 15),
#'   D = c(16, 42, 18, 23, 20),
#'   E = c(12, 54, 23, 33, 42)
#' )
#' 
#' data_df2 <- data.frame(
#'   Tom = c(11, 27, 23, 29, 30),
#'   Anny = c(31, 16, 33, 98, 35),
#'   Jerry = c(36, 37, 45, 39, 40),
#'   Cat = c(41, 21, 43, 44, 45),
#'   White = c(46, 77, 48, 94, 12)
#' )
#' re <- CorRcpp(data_df1,data_df2,type = 'spearman')
#' re$R_matrix # correlation value
#' re$P_matrix # pvalue
CorRcpp <- function(x=x,y=NULL,type=c("pearson","spearman")) {
  stopifnot(is.data.frame(x))
	type <- match.arg(type)
  if(is.null(y)){
    if(nrow(x)<4){
      stop("Cor-analysis need more than 4 samples")
    }
    if(type=="spearman"){
      x <- apply(x, 2, rank)
    }
    corres <-  cp_cor_s(mat = as.matrix(x))
    corres$R_matrix <- as.data.frame(corres$R_matrix)
    corres$P_matrix <- as.data.frame(corres$P_matrix)
    diag(corres$P_matrix) <- 0 ## advoid the NaN value in case of the same feature
    row.names(corres$R_matrix) <- colnames(x)
    row.names(corres$P_matrix) <- colnames(x)
    colnames(corres$R_matrix) <- colnames(x)
    colnames(corres$P_matrix) <- colnames(x)
  }else{
    if(nrow(x)<4 | nrow(y)<4){
      stop("Cor-analysis need more than 4 samples")
    }
    if(type=="spearman"){
      stopifnot(is.data.frame(y))
      x <- apply(x, 2, rank)
      y <- apply(y, 2, rank)
    }

    corres <-  cp_cor_t(mat = as.matrix(x),mat2=as.matrix(y))
    corres$R_matrix <- as.data.frame(corres$R_matrix)
    corres$P_matrix <- as.data.frame(corres$P_matrix)
    if(check_xy_duplicate(x,y)){
      diag(corres$P_matrix) <- 0 ## advoid the NaN value in case of the same feature
    }
    row.names(corres$R_matrix) <- colnames(x)
    row.names(corres$P_matrix) <- colnames(x)
    colnames(corres$R_matrix) <- colnames(y)
    colnames(corres$P_matrix) <- colnames(y)
  }
  return(corres)
}


check_xy_duplicate <- function(x,y) {
  rownames(x) <- colnames(x) <- NULL
  rownames(y) <- colnames(y) <- NULL
  identical(x, y)
}