
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


#' @importFrom stringr str_detect
#' @importFrom SummarizedExperiment SummarizedExperiment
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
  temp <- temp[rowSums(temp[,-1]) != 0,] # filter away empty feature!
  rownames(temp) <- NULL # necessary!
  temp2 <- dplyr::left_join(temp,ref,by = 'feature')
  temp_rowdata <- temp2 %>% dplyr::select(feature,Name)
  temp %<>% tibble::column_to_rownames('feature') %>% as.matrix()
  deposit <- SummarizedExperiment(assays=list(counts=temp),
                                  rowData = temp_rowdata)
  return(deposit)

}


#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom stringr str_replace_all
humann_taxonomy_import <- function(file=NULL,data=NULL,sep = '|') {
  feature <- `.` <- NULL
  if (!is.null(data)) {
    temp <- data
  }else {
    temp <- read.table(file,header = T,sep = '\t',quote="")
  }  
  colnames(temp)[1] <- 'feature'
  temp%<>%dplyr::filter(stringr::str_detect(feature,'\\|t_'))
  temp <- temp[rowSums(temp[,-1]) != 0,] # filter away empty feature!
  rownames(temp) <- NULL # necessary!
  temp %<>%
      dplyr::mutate(feature = stringr::str_replace_all(feature, " ", "_"))  # Space in the value will lead to unexperted error  
  temp %>% dplyr::pull(feature) %>% read.table(text = .,sep = sep) -> temp_name
  colnames(temp_name) <- c('Kindom','Phylum','Class','Order','Family','Genus','Species','Strain')[1:ncol(temp_name)]
  temp_name <- data.frame(feature = temp$feature,temp_name) %>%
    .impute_tax() ## impute the NA tax
  temp %<>% tibble::column_to_rownames('feature') %>% as.matrix()
  deposit <- SummarizedExperiment(assays=list(counts=temp),
                                  rowData = temp_name)
  return(deposit)
}

#' Import microbial data into SummariseExperiment
#'
#' @param file A file path.
#' @param data A dataframe.The row must be the feature and the column is the sample.
#' @param humann_format A boolean. Whether the function improt the data according to the humann format.
#' @param assay_name A character string. Indicate what kind of result the data belongs to, such as counts, relative abundance, TPM, etc.
#' @param sep The field separator character. Values on feature column of the file are separated by this character. (defacult:';',when humann_format=T defacult:'|')
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom stringr str_replace_all
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
    temp <- temp[rowSums(temp[,-1]) != 0,] # filter away empty feature!
    rownames(temp) <- NULL # necessary!
    temp %<>%
      dplyr::mutate(feature = stringr::str_replace_all(feature, " ", "_"))  # Space in the value will lead to unexperted error 
    temp_name <- temp %>% dplyr::pull(feature) %>% read.table(text = .,sep = sep) %>% 
      dplyr::mutate_if(~ any(. == ""), ~ dplyr::na_if(., ""))  # make the empty value into NA in order to adapt .impute_tax
    colnames(temp_name) <- c('Kindom','Phylum','Class','Order','Family','Genus','Species','Strain')[1:ncol(temp_name)]
    temp_name <- data.frame(feature = temp$feature,temp_name) %>%
      .impute_tax() ## impute the NA tax
    temp  %<>% tibble::column_to_rownames('feature') %>% as.matrix()
    deposit <- SummarizedExperiment(assays=list(counts=temp),
                                    rowData = temp_name)
  }
  if (!is.null(assay_name)) {
    assayNames(deposit, 1) <- assay_name
  }
  return(deposit)
}


#' Import gene data into SummariseExperiment
#'
#' @param file A file path. The file should be deposited in txt format.
#' @param data A dataframe.The row must be the feature and the column is the sample.
#' @param type type
#' @param assay_name A character string. Indicate what kind of result the data belongs to, such as counts, relative abundance, TPM, etc.
#' @param humann_format A boolean. Whether the function improt the data according to the humann format.
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SummarizedExperiment SummarizedExperiment
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
    temp <- temp[rowSums(temp[,-1]) != 0,] # filter away empty feature!    
    rownames(temp) <- NULL # necessary!
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


#' Import combined data into SummariseExperiment
#'
#' @param file A file path. The file should be deposited in txt format.
#' @param data A dataframe.The row must be the feature and the column is the sample.
#' @param sampleID a seris of string character.
#' @param dfmap A dataframe. Indicate the experiment name, sample source and sample tube details in the omics data.
#' @param assay_name A character string. Indicate what kind of result the data belongs to, such as counts, relative abundance, TPM, etc.
#' @param assay A character string. Indicate the experiment name of the import data in the dfmap.
#' @importFrom dplyr all_of
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @return SummmariseExperiment object
#' @export
#'
#' @examples
#' # example
EMP_normal_import <- function(file=NULL,data=NULL,sampleID=NULL,dfmap=NULL,assay_name=NULL,assay=NULL){
  row_data <- assay_data <- obj <- colname <- feature <- NULL
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
  colnames(data)[1] <- 'feature'
  data <- data[rowSums(data[,sampleID]) != 0,] # filter away empty feature!
  rownames(data) <- NULL # necessary!
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

#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom MultiAssayExperiment MultiAssayExperiment
EMP_easy_normal_import <- function(file=NULL,data=NULL,assay='experiment',assay_name=NULL,coldata=NULL,output='MAE') {
  sampleID <- row_data <- assay_data <- SE_object <- feature <- Name <- NULL
  if (!is.null(data)) {
    data <- data
  }else {
    data <- read.table(file=file,header = T,sep = '\t',quote="")
  }  
  
  colnames(data)[1] <- 'feature'
  data <- data[rowSums(data[,-1]) != 0,] # filter away empty feature!
  rownames(data) <- NULL # necessary!
  row_data <- data %>% 
    dplyr::mutate(Name=feature) %>%
    dplyr::select(feature,Name)
  
  assay_data <- data %>%
    tibble::column_to_rownames('feature')
  sampleID <- colnames(assay_data)
  
  SE_object <- SummarizedExperiment(assays=list(counts= as.matrix(assay_data)),
                              rowData = row_data)
  if (!is.null(assay_name)) {
    assayNames(SE_object, 1) <- assay_name
  }
  
  if (output == 'SE') {
    return(SE_object)
  }else if(output=='MAE'){
   if (is.null(coldata)) {
     stop("If output is MultiAssayExperiment, please input the coldata contaning the group and other information!")
   }
  
   
   if (ncol(coldata) == 0) {
     stop("coldata must contain one column informantion at least!")
   }else if (ncol(coldata) == 1){
     coldata <- coldata %>% as.data.frame() %>% dplyr::mutate(colname = sampleID)
   }else{
     coldata <- coldata
   }
   
   objlist <- list(SE_object)
   names(objlist)[1] <- assay
   
   dfmap <- data.frame(assay=assay,primary=sampleID,colname=sampleID) %>% 
     dplyr::mutate(assay = factor(assay))
   
   MAE_object <- MultiAssayExperiment(objlist, coldata, dfmap)
    
   return(MAE_object)
  }else{
    stop("output only support SE and MAE")
  }
}

#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom MultiAssayExperiment MultiAssayExperiment
EMP_easy_taxonomy_import <- function(file=NULL,data=NULL,assay='experiment',assay_name=NULL,coldata=NULL,
                                     humann_format=FALSE,output='MAE',
                                     sep=if (humann_format == "FALSE") ';' else '|') {
  
  sampleID <- assay_data <- SE_object <- NULL

  SE_object <- EMP_taxonomy_import(file=file,data=data,humann_format=humann_format,assay_name=assay_name,sep=sep)                                    
  
  sampleID <-  dimnames(SE_object)[[2]]                                 
  if (output == 'SE') {
    return(SE_object)
  }else if(output=='MAE'){
    if (is.null(coldata)) {
      stop("If output is MultiAssayExperiment, please input the coldata contaning the group and other information!")
    }
    
    if (ncol(coldata) == 0) {
      stop("coldata must contain one column informantion at least!")
    }else if (ncol(coldata) == 1){
      coldata <- coldata %>% as.data.frame() %>% dplyr::mutate(colname = sampleID)
    }else{
      coldata <- coldata
    }
    
    objlist <- list(SE_object)
    names(objlist)[1] <- assay
    
    dfmap <- data.frame(assay=assay,primary=sampleID,colname=sampleID) %>% 
      dplyr::mutate(assay = factor(assay))
    
    MAE_object <- MultiAssayExperiment::MultiAssayExperiment(objlist, coldata, dfmap)
    
    return(MAE_object)
  }else{
    stop("output only support SE and MAE")
  }
}


#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom MultiAssayExperiment MultiAssayExperiment
EMP_easy_function_import <- function(file=NULL,data=NULL,type,assay='experiment',assay_name=NULL,coldata=NULL,
                                     humann_format=FALSE,output='MAE') {
  
  sampleID <- assay_data <- SE_object <- NULL
  SE_object <- EMP_function_import(file=file,data=data,type,humann_format=humann_format,assay_name=assay_name)                                    
  
  sampleID <-  dimnames(SE_object)[[2]]                                 
  if (output == 'SE') {
    return(SE_object)
  }else if(output=='MAE'){
    if (is.null(coldata)) {
      stop("If output is MultiAssayExperiment, please input the coldata contaning the group and other information!")
    }
    
    if (ncol(coldata) == 0) {
      stop("coldata must contain one column informantion at least!")
    }else if (ncol(coldata) == 1){
      coldata <- coldata %>% as.data.frame() %>% dplyr::mutate(colname = sampleID)
    }else{
      coldata <- coldata
    }
    
    objlist <- list(SE_object)
    names(objlist)[1] <- assay
    
    dfmap <- data.frame(assay=assay,primary=sampleID,colname=sampleID) %>% 
      dplyr::mutate(assay = factor(assay))
    
    MAE_object <- MultiAssayExperiment::MultiAssayExperiment(objlist, coldata, dfmap)
    
    return(MAE_object)
  }else{
    stop("output only support SE and MAE")
  }
}


#' Easily import data into SummariseExperiment or MultiAssayExperiment
#'
#' @param file A file path. The file should be deposited in txt format.
#' @param data A dataframe.The row must be the feature and the column is the sample.
#' @param type A character string. Methods include tax, ko, ec and normal.
#' @param assay A character string. Set the experiment name.(default:experiment)
#' @param assay_name A character string. Indicate what kind of result the data belongs to, such as counts, relative abundance, TPM, etc.
#' @param humann_format A boolean. Whether the function improt the data according to the humann format.
#' @param coldata A dataframe containing one column informantion at least.
#' @param sep The field separator character. Only activated when type =''tax. Values on each line of the file are separated by this character. (defacult:'|') 
#' @param output A character string. Set the output result in SE(SummariseExperiment) or MAE(MultiAssayExperiment) format.
#' @rdname EMP_easy_import
#' @return SummmariseExperiment or MultiAssayExperiment object
#' @export
#'
#' @examples
#' # example

EMP_easy_import <- function(file=NULL,data=NULL,type,assay='experiment',assay_name=NULL,coldata=NULL,
                            humann_format=FALSE,output='MAE',
                            sep=if (humann_format == "FALSE") ';' else '|'){
  obj <- NULL
  switch(type,
         "tax" = {obj <- EMP_easy_taxonomy_import(file=file,data=data,assay=assay,assay_name=assay_name,coldata=coldata,
                                                  humann_format=humann_format,output=output,sep=sep)},
         "ko" = {obj <- EMP_easy_function_import(file=file,data=data,type,assay=assay,assay_name=assay_name,coldata=coldata,
                                                  humann_format=humann_format,output=output)},
         "ec" = {obj <- EMP_easy_function_import(file=file,data=data,type,assay=assay,assay_name=assay_name,coldata=coldata,
                                                  humann_format=humann_format,output=output)},
         "normal" = {obj <- EMP_easy_normal_import(file=file,data=data,assay=assay,assay_name=assay_name,coldata=coldata,output=output)},
         {
           stop('type only support tax, ko, ec and normal!')
         }
  )
  return(obj)
}

