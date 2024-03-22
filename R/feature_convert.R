
#' Title
#'
#' @param feature features
#' @param from from
#' @param to to
#' @param species species 
#' @param OrgDb OrgDb
#'
#' @return data.frame
#' @export
#'
#' @examples
#' #
.feature_convert_gene <- function(feature, from = "SYMBOL", to = "ENTREZID", species = "Human", OrgDb = NULL) {
    if (!(species %in% c("Human", "Mouse", "Pig", "Zebrafish"))) {
        if (is.null(OrgDb)) {
            stop("The species is not within the built-in species range, OrgDb needs to be provided for conversion.")
        }
        result <- AnnotationDbi::select(x = OrgDb, keys = feature, keytype = from, columns = c(to))
    }
    
    data_df <- switch(species,
                      "Human" = res_Hs,
                      "Mouse" = res_Mm,
                      "Pig" = res_Ss,
                      "Zebrafish" = res_Dr
                      )
    from <- match.arg(from, colnames(data_df))
    to <- match.arg(to, colnames(data_df))
    result <- data_df[data_df[, from] %in% feature, c(from, to)]
    na_feature <- setdiff(feature, result[, 1])
    if (length(na_feature) > 0) {
        result2 <- data.frame(from = na_feature, to = NA)
        colnames(result2) <- colnames(result)
        result <- rbind(result, result2)
    }
    result <- dplyr::distinct(result)
    return(result)  
}    



#' Title
#'
#' @param EMPT EMPT
#' @param from from
#' @param to to
#' @param method method 
#'
#' @return EMPT object
#' @export
#'
#' @examples
#' #

.EMP_feature_convert_cpd <- function(EMPT,method='mean',from,to){
  raw_rowdata <- .get.row_info.EMPT(EMPT)
  ref_data <- metaboliteIDmapping_data |> 
    dplyr::select({{from}},{{to}}) |>
    dplyr::rename(feature = {{from}}) |>
    tidyr::drop_na(feature) |>
    dplyr::distinct(feature,.keep_all = TRUE)
  new_rowdata <- dplyr::left_join(raw_rowdata,ref_data,by='feature')
  .get.row_info.EMPT(EMPT) <- new_rowdata
  deposit <- EMPT |> EMP_collapse(estimate_group = to,method = method,collapse_by='row')
  return(deposit)
}


#' Title
#'
#' @param feature features
#' @param from from
#' @param to to
#' @param species species 
#' @param OrgDb OrgDb
#'
#' @return data.frame
#' @export
#'
#' @examples
#' #
.EMP_feature_convert_gene <- function(EMPT,method='mean',from,to,species = "Human", OrgDb = NULL) {
  raw_feature <- .get.row_info.EMPT(EMPT) %>% dplyr::pull(feature)
  ref_data <- .feature_convert_gene(feature=raw_feature, 
                                       from = from, 
                                       to = to, 
                                       species = species, 
                                       OrgDb=OrgDb)
  colnames(ref_data) <- c('feature',to)
  
  raw_rowdata <- .get.row_info.EMPT(EMPT)
  new_rowdata <- dplyr::left_join(raw_rowdata,ref_data,by='feature')
  .get.row_info.EMPT(EMPT) <- new_rowdata
  deposit <- EMPT |> EMP_collapse(estimate_group = to,method = method,collapse_by='row')
  return(deposit)
}


