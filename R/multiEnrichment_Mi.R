#' Mutual information
#'
#' @param genelist1 genelist1
#' @param genelist2 genelist2
#' @importFrom infotheo entropy
#' @export
r <- function(genelist1, genelist2) {
## https://www.cnblogs.com/jiangyaling/p/8040024.html
  H1 <- entropy(genelist1)
  H2 <- entropy(genelist2)
  H12 <- entropy(cbind(genelist1,genelist2))
  mi <- H1 + H2 - H12
  r <- mi/max(H1, H2)
}


#' Conditional mutual information
#'
#' @param genelist1 genelist1
#' @param genelist2 genelist2
#' @param genelist3 genelist3
#' @export
cmi <- function(genelist1,genelist2,genelist3){
  H1 <- entropy(genelist1)
  H2 <- entropy(genelist2)
  H3 <- entropy(genelist3)
  H13 <- entropy(cbind(genelist1,genelist3))
  H23 <- entropy(cbind(genelist2,genelist3))
  H123 <- entropy(cbind(genelist1,genelist2,genelist3))
  (H13 + H23 - H3 - H123) / (max(H1, H2))
}

