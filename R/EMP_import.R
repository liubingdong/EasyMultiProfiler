
##' @importFrom yulab.utils yread
kegg_rest <- function(rest_url) {
    `.` <- NULL
    message('Reading KEGG annotation online: "', rest_url, '"...')
    # content <- readLines(f)
    content <- yread(rest_url)
    content %<>% strsplit(., "\t") %>% do.call('rbind', .)
    res <- data.frame(from=content[,1],
                      to=content[,2])
    return(res)
}



#' Title
#'
#' @param x wait_for_add
#' @param type wait_for_add
#' @importFrom memoise forget
#' @importFrom stringr str_detect
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
humann_function_import <- function(x,type) {
  from <- to <- feature <- Name <- `.` <- NULL
  temp <- read.table(x,header = T,sep = '\t',quote="")
  if (type== 'KO') {
    ref = kegg_rest("https://rest.kegg.jp/list/ko") %>%
  dplyr::rename(feature=from,Name = to)
  }else if (type == 'EC') {
    ref = kegg_rest("https://rest.kegg.jp/list/ec") %>%
  dplyr::rename(feature=from,Name = to)
  }
  colnames(temp)[1] <- 'feature'
  temp %<>% dplyr::filter(!stringr::str_detect(feature,'\\|') & !stringr::str_detect(feature,'UN'))
  temp2 <- dplyr::left_join(temp,ref,by = 'feature')
  temp_rowdata <- temp2 %>% dplyr::select(feature,Name)
  temp %<>% tibble::column_to_rownames('feature') %>% as.matrix()
  deposit <- SummarizedExperiment(assays=list(counts=temp),
                                  rowData = temp_rowdata)
  return(deposit)

}

#' Title
#'
#' @param x wait_for_add
#' @param sep wait_for_add
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
humann_taxonomy_import <- function(x,sep = '|') {
  feature <- `.` <- NULL
  temp <- read.table(x,header = T,sep = '\t',quote="")
  colnames(temp)[1] <- 'feature'
  temp%<>%dplyr::filter(stringr::str_detect(feature,'\\|t_'))
  temp %>% dplyr::pull(feature) %>% read.table(text = .,sep = sep) -> temp_name
  colnames(temp_name) <- c('Kindom','Phylum','Class','Order','Family','Genus','Species','Strain')[1:ncol(temp_name)]
  temp_name <- data.frame(feature = temp$feature,temp_name)
  temp %<>% tibble::column_to_rownames('feature') %>% as.matrix()
  deposit <- SummarizedExperiment(assays=list(counts=temp),
                                  rowData = temp_name)
  return(deposit)
}

#' Title
#'
#' @param x wait_for_add
#' @param humann_format wait_for_add
#' @param assay_name wait_for_add
#' @param sep wait_for_add
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_taxonomy_import <- function(x,humann_format=FALSE,assay_name=NULL,sep=if (humann_format == "FALSE") ';' else '|') {
  feature <- `.` <- NULL
  if(humann_format == TRUE){
    deposit <- humann_taxonomy_import(x,sep = sep)
  }else{
    temp <- read.table(x,header = T,sep = '\t',quote="")
    colnames(temp)[1] <- 'feature'
    temp %>% dplyr::pull(feature) %>% read.table(text = .,sep = sep) -> temp_name
    colnames(temp_name) <- c('Kindom','Phylum','Class','Order','Family','Genus','Species','Strain')[1:ncol(temp_name)]
    temp_name <- data.frame(feature = temp$feature,temp_name)
    temp  %<>% tibble::column_to_rownames('feature') %>% as.matrix()
    deposit <- SummarizedExperiment(assays=list(counts=temp),
                                    rowData = temp_name)
  }
  if (!is.null(assay_name)) {
    assayNames(deposit, 1) <- assay_name
  }
  return(deposit)
}


#' Title
#'
#' @param x wait_for_add
#' @param type wait_for_add
#' @param assay_name wait_for_add
#' @param humann_format wait_for_add
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SummarizedExperiment `assayNames<-`
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_function_import <- function(x,type,assay_name=NULL,humann_format=FALSE) {
  from <- to <- feature <- Name <- NULL
  if(humann_format == T){
    deposit <- humann_function_import(x,type)
  }else {
    temp <- read.table(x,header = T,sep = '\t',quote="")
    if (type== 'KO') {
      ref = kegg_rest("https://rest.kegg.jp/list/ko") %>%
    dplyr::rename(feature=from,Name = to)
    }else if (type == 'EC') {
      ref = kegg_rest("https://rest.kegg.jp/list/ec") %>%
    dplyr::rename(feature=from,Name = to)
    }
    colnames(temp)[1] <- 'feature'
    temp2 <- dplyr::left_join(temp,ref,by = 'feature')
    temp_rowdata <- temp2 %>% dplyr::select(feature,Name)
    temp %<>% tibble::column_to_rownames('feature') %>% as.matrix()
    deposit <- SummarizedExperiment(assays=list(counts=temp),
                                    rowData = temp_rowdata)
  }
  if (!is.null(assay_name)) {
    assayNames(deposit, 1) <- assay_name
  }
  return(deposit)
}


#' Title
#'
#' @param x wait_for_add
#' @param sampleID wait_for_add
#' @param dfmap wait_for_add
#' @param assay_name wait_for_add
#' @param assay wait_for_add
#' @importFrom dplyr all_of
#'
#' @return xx object
#' @export
#'
#' @examples
#' # example
EMP_normal_import <- function(x,sampleID=NULL,dfmap=NULL,assay_name=NULL,assay=NULL){
  colname <- feature <- NULL
  if (!is.null(sampleID)) {
    sampleID <- sampleID
  }else{
    if (is.null(dfmap) | is.null(assay) ) {
      stop("If sampleID is NULL, please input dfmap and assay!")
    }else{
      sampleID <- dfmap %>% as.data.frame() %>% dplyr::filter(assay == !!assay) %>%
        dplyr::pull(colname)
    }
  }
  data <- read.table(x,header = T,sep = '\t',quote="")
  row_data <- data %>% dplyr::select(!all_of(!!sampleID))
  assay_data <- data %>% dplyr::select(feature,all_of(!!sampleID)) %>%
    tibble::column_to_rownames('feature')
  obj <- SummarizedExperiment(assays=list(counts= as.matrix(assay_data)),
                              rowData = row_data)
  if (!is.null(assay_name)) {
    assayNames(deposit, 1) <- assay_name
  }
  return(obj)
}


