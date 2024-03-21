
#' Title
#'
#' @param feature features
#' @param from from
#' @param to to
#' @param species speckes 
#' @param OrgDb OrgDb
#'
#' @return data.frame
#' @export
#'
#' @examples
#' #
feature_convert <- function(feature, from = "SYMBOL", to = "ENTREZID", species = "Human", OrgDb = NULL) {
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
    return(result)  
}    
