
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
#' @param file A file path.
#' @param data A dataframe.The row must be the feature and the column is the sample.
#' @param type A character string. Methods include ko and ec.
#' @importFrom stringr str_detect
#'
#' @return SummmariseExperiment object
#' @export
#'
#' @examples
#' # add example
humann_function_import <- function(file=NULL,data=NULL,type) {
  from <- to <- feature <- Name <- `.` <- NULL
  if (!is.null(data)) {
    temp <- data
  }else {
    temp <- read.table(file,header = T,sep = '\t',quote="")
  }
  if (type== 'ko') {
    ref = kegg_rest("https://rest.kegg.jp/list/ko") %>%
  dplyr::rename(feature=from,Name = to)
  }else if (type == 'ec') {
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
#' @param file A file path.
#' @param data A dataframe.The row must be the feature and the column is the sample.
#' @param sep The field separator character. Values on each line of the file are separated by this character. (defacult:'|')
#'
#' @return SummmariseExperiment object
#' @export
#'
#' @examples
#' # add example
humann_taxonomy_import <- function(file=NULL,data=NULL,sep = '|') {
  feature <- `.` <- NULL
  if (!is.null(data)) {
    temp <- data
  }else {
    temp <- read.table(file,header = T,sep = '\t',quote="")
  }  
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
#' @param file A file path.
#' @param data A dataframe.The row must be the feature and the column is the sample.
#' @param humann_format A boolean. Whether the function improt the data according to the humann format.
#' @param assay_name A character string. Indicate what kind of result the data belongs to, such as counts, relative abundance, TPM, etc.
#' @param sep The field separator character. Values on each line of the file are separated by this character. (defacult:'|')
#'
#' @return SummmariseExperiment object
#' @export
#'
#' @examples
#' # add example
EMP_taxonomy_import <- function(file=NULL,data=NULL,humann_format=FALSE,assay_name=NULL,sep=if (humann_format == "FALSE") ';' else '|') {
  feature <- `.` <- NULL
  if(humann_format == TRUE){
    deposit <- humann_taxonomy_import(file=file,data=data,sep = sep)
  }else{
    if (!is.null(data)) {
      temp <- data
    }else {
      temp <- read.table(file=file,header = T,sep = '\t',quote="")
    }     
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
#' @param file A file path.
#' @param data A dataframe.The row must be the feature and the column is the sample.
#' @param assay_name A character string. Indicate what kind of result the data belongs to, such as counts, relative abundance, TPM, etc.
#' @param humann_format A boolean. Whether the function improt the data according to the humann format.
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SummarizedExperiment `assayNames<-`
#'
#' @return SummmariseExperiment object
#' @export
#'
#' @examples
#' # add example
EMP_function_import <- function(file=NULL,data=NULL,type,assay_name=NULL,humann_format=FALSE) {
  from <- to <- feature <- Name <- NULL
  if(humann_format == T){
    deposit <- humann_function_import(file=file,data=data,type=type)
  }else {
    if (!is.null(data)) {
      temp <- data
    }else {
      temp <- read.table(file=file,header = T,sep = '\t',quote="")
    }     
    if (type== 'ko') {
      ref = kegg_rest("https://rest.kegg.jp/list/ko") %>%
    dplyr::rename(feature=from,Name = to)
    }else if (type == 'ec') {
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
#' @param file A file path.
#' @param data A dataframe.The row must be the feature and the column is the sample.
#' @param sampleID a seris of string character.
#' @param dfmap A dataframe. Indicate the experiment name, sample source and sample tube details in the omics data.
#' @param assay_name A character string. Indicate what kind of result the data belongs to, such as counts, relative abundance, TPM, etc.
#' @param assay A character string. Indicate the experiment name of the import data in the dfmap.
#' @importFrom dplyr all_of
#'
#' @return SummmariseExperiment object
#' @export
#'
#' @examples
#' # example
EMP_normal_import <- function(file=NULL,data=NULL,sampleID=NULL,dfmap=NULL,assay_name=NULL,assay=NULL){
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
  if (!is.null(data)) {
    data <- data
  }else {
    data <- read.table(file=file,header = T,sep = '\t',quote="")
  }   
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


