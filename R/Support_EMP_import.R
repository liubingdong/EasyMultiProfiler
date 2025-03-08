
##' @importFrom yulab.utils yread
kegg_rest <- function(rest_url) {
    `.` <- NULL
    #info_output <- paste0('Reading KEGG annotation online: "', rest_url, '"...')
    #EMP_message(info_output,color = 32,order = 1,show='message')
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
  
  if (any(is.na(temp[,-1]))) {
    EMP_message("The NA value has been detected in the data and changed into 0!",color = 32,order = 1,show='message')
    temp[,-1][is.na(temp[,-1])] <- 0
  }

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

  if (any(is.na(temp[,-1]))) {
    EMP_message("The NA value has been detected in the data and changed into 0!",color = 32,order = 1,show='message')
    temp[,-1][is.na(temp[,-1])] <- 0
  }

  temp <- temp[rowSums(temp[,-1]) != 0,] # filter away empty feature!
  rownames(temp) <- NULL # necessary!
  temp %<>%
      dplyr::mutate(feature = stringr::str_replace_all(feature, " ", "_"))  # Space in the value will lead to unexperted error  
  temp %>% dplyr::pull(feature) %>% read.table(text = .,sep = sep) -> temp_name
  colnames(temp_name) <- c('Kindom','Phylum','Class','Order','Family','Genus','Species','Strain')[1:ncol(temp_name)]
  
  #temp_name <- double_tax_name(temp_name,sep=';')
  
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
#' @param file_format A string including biom and gzv.
#' @param humann_format A boolean. Whether the function import the data in the Metaphlan or Humann format.
#' @param start_level A string. Specific the start level of input data from Domain,Kindom,Phylum,Class,Order,Family,Genus,Species,Strain.(default:Kindom)
#' @param duplicate_feature A boolean. Whether the feature exist the dupicated name.
#' @param assay_name A character string. Indicate what kind of result the data belongs to, such as counts, relative abundance, TPM, etc.
#' @param sep The field separator character. Values on feature column of the file are separated by this character. (defacult:';',when humann_format=T defacult:'|')
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom stringr str_replace_all
#' @importFrom biomformat biom
#' @importFrom biomformat read_hdf5_biom
#' @importFrom biomformat observation_metadata
#' @importFrom biomformat biom_data
#' @importFrom utils read.csv
#' @importFrom utils unzip
#' @importFrom dplyr any_of
#' @importFrom dplyr contains
#' @return SummmariseExperiment object
#' @details Paramter file_format and humann_format help the function import data properly. Data in humann format is usually is generated from Metaphlan and Humann. Data in biom format is usually is generated from Qiime1. Data in qzv format is usually is generated from Qiime2.
#' @export
#'
#' @examples
#' # More examples and tutorial could be found on website: 
#' # https://liubingdong.github.io/EasyMultiProfiler_tutorial/
EMP_taxonomy_import <- function(file=NULL,data=NULL,humann_format=FALSE,file_format=NULL,start_level='Kindom',assay_name=NULL,duplicate_feature=NULL,sep=if (humann_format == "TRUE") '|' else ';') {
  feature <- `.` <- check_duplicated_feature <- biom_df <- biom_data <- tax_data <- otuid <- unzipfiles <- data_file <- NULL
  if(humann_format == TRUE){
    deposit <- humann_taxonomy_import(file=file,data=data,sep = sep)
  }else{
    if (!is.null(data)) {
      temp <- data
    }else{
      if (!is.null(file_format)) {
        switch(file_format,
               biom={
                 biom_data <- biomformat::biom(biomformat::read_hdf5_biom(biom_file=file)) 
                 biom_df <- biom_data %>%
                   biomformat::biom_data() %>%
                   as.matrix() %>% 
                   as.data.frame() %>%
                   tibble::rownames_to_column('otuid')
                 tax_data <- biomformat::observation_metadata(biom_data) %>% 
                   dplyr::select(dplyr::starts_with('taxonomy')) %>%
                   tidyr::unite(col = "feature", everything(), sep = ";") %>%
                   tibble::rownames_to_column('otuid')
                 temp <- dplyr::left_join(biom_df,tax_data,by='otuid') %>%
                   dplyr::select(-otuid) %>%
                   dplyr::select(feature,everything()) 
               },
               qzv={
                 tmpdir <- tempdir()
                 unzipfiles <- unzip(file, exdir=tmpdir)
                 data_file <- unzipfiles[grep("level-7.csv", unzipfiles)]
                 temp <- read.csv(data_file,header = T,row.names=1,check.names=FALSE) %>%  # check.names is to properly read ; in column names
                   t() %>% as.data.frame() %>%
                   tibble::rownames_to_column('feature') 
               }, 
               {
                  stop('Parameter file_format must be NULL, biom or qzv')
               }
 
        )        
      }else{
        # the process below will support the level-7 microbial data and ASV/OTU raw table generated from biom
        lines <- readLines(file,warn = FALSE)
        ## if the first line with #, the line will be ignored
        if (substring(lines[1], 1, 1) == "#") {
          lines <- lines[-1]
        }
        ## use textConnection() to pass the content to the function read.table()
        temp <- read.table(text = paste(lines, collapse = "\n"), header = TRUE,sep = '\t',comment.char="") %>%
          dplyr::select(-any_of(contains(c('OTU ID','OTU.ID','OTU_ID','ASV ID','ASV.ID','ASV_ID')))) %>%
          dplyr::select(where(is.character),everything())        
      }
  }

  colnames(temp)[1] <- 'feature'

  temp[["feature"]] <- gsub(paste0(sep,' '), sep, temp[["feature"]])
  temp[["feature"]] <- gsub(paste0(' ',sep), sep, temp[["feature"]])

  if (any(is.na(temp[,-1]))) {
    EMP_message("The NA value has been detected in the data and changed into 0!",color = 32,order = 1,show='message')
    temp[,-1][is.na(temp[,-1])] <- 0
  }

  temp <- temp[rowSums(temp[,-1]) != 0,] # filter away empty feature!
  rownames(temp) <- NULL # necessary!
  
  ## Mistake-proofing
  strings_to_remove1 <- c('d__','k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__') # __ must go first than _
  strings_to_remove2 <- c('d_','k_', 'p_', 'c_', 'o_', 'f_', 'g_', 's_')

  total_tax_info <- c('Domain','Kindom','Phylum','Class','Order','Family','Genus','Species','Strain')
  idx <- which(total_tax_info %in% start_level) 
  strings_to_remove1_append <- paste0(sep,strings_to_remove1[(idx + 1) : length(strings_to_remove1)])
  strings_to_remove2_append <- paste0(sep,strings_to_remove2[(idx + 1) : length(strings_to_remove2)])

  for (i in strings_to_remove1_append) {
    temp[["feature"]] <- gsub(i, sep, temp[["feature"]])
  }
  temp[["feature"]] <- gsub(paste0("^", strings_to_remove1[idx]), "", temp[["feature"]])

  for (i in strings_to_remove2_append) {
    temp[["feature"]] <- gsub(i, sep, temp[["feature"]])
  }
  temp[["feature"]] <- gsub(paste0("^", strings_to_remove2[idx]), "", temp[["feature"]])

  # Delete the exception of tax anotation for silva in the level of Kindom
  temp[["feature"]] <- gsub(paste0("^", "d__"), "", temp[["feature"]]) # __ must go first than _
  temp[["feature"]] <- gsub(paste0("^", "d_"), "", temp[["feature"]])

  # Delete any other exceptions or empty tax annotations
  temp[["feature"]] <- gsub(paste0(sep,'__'), sep, temp[["feature"]])
  temp[["feature"]] <- gsub(paste0(sep,'_'), sep, temp[["feature"]])
  temp[["feature"]] <- gsub(paste0('__',sep), sep, temp[["feature"]])
  temp[["feature"]] <- gsub(paste0('_',sep), sep, temp[["feature"]])


  temp %<>%
    dplyr::mutate(feature = stringr::str_replace_all(feature, " ", "_")) 
  
  if (is.null(duplicate_feature)){
    check_duplicated_feature <- any(duplicated(temp$feature))
    if(check_duplicated_feature){
      duplicate_feature <- TRUE
      EMP_message("The duplicated name in the feature have been detected.\nParameter duplicate_feature is forced to be TRUE.",color = 32,order = 1,show='message')
    }else{
      duplicate_feature <- FALSE
    }
  }
  
  if(duplicate_feature==FALSE){
    
    temp_name <- temp %>% dplyr::pull(feature) %>% read.table(text = .,sep = sep,blank.lines.skip=F,quote = "",row.names = NULL,header = FALSE) %>% 
      dplyr::mutate_if(~ any(. == ""), ~ dplyr::na_if(., "")) 

    #total_tax_info <- c('Domain','Kindom','Phylum','Class','Order','Family','Genus','Species','Strain')
    real_tax_info <- total_tax_info[match(start_level, total_tax_info):length(total_tax_info)]
    colnames(temp_name) <- real_tax_info[1:ncol(temp_name)]

    ## make the dupicate anotation into NA
    #strings_to_remove3 <- c(
    #  "uncultured", "unculturebacterium", "uncultureorganism", "unidentified", "unculturemarine",
    #  "uncultureLactobacillus", "uncultureprokaryote", "uncultureBacteroidales", "uncultureDesulfovibrionaceae",
    #  "uncultureLachnospiraceae", "uncultureErysipelotrichales", "unculturePrevotella", "uncultureAnaerotruncus",
    #  "uncultureRuminococcus", "uncultureWautersiella", "uncultureCandidatus", "uncultureEubacterium",
    #  "uncultureClostridiales", "unculturecompost", "uncultureActinobacteridae", "uncultureClostridiaceae",
    #  "uncultureOdoribacter", "uncultureCoxiella", "uncultureactinobacterium", "uncultureErysipelotrichaceae",
    #  "uncultureAnaerolineae", "uncultureBacteroidaceae", "uncultureAlphaproteobacteria", "unculturerumen",
    #  "uncultureRhodospirillaceae", "unculturesoil", "unculturespirochete", "uncultureMethanobacteriales"
    #)
    #occurrences <- colSums(sapply(strings_to_remove3, function(pattern) grepl(pattern, unlist(temp_name))))
    #strings_to_remove3_real <- names(occurrences[occurrences >1]) # search the dupicate anotation which appear more than once
    #temp_name <- apply(temp_name, 2, function(x) {
    #  for (string in strings_to_remove3_real) {
    #    x <- gsub(string, NA, x)       
    #  }
    #  return(x)
    #})
    #temp_name <- as.data.frame(temp_name)    

    temp_name <- data.frame(feature = temp$feature,temp_name) %>%
      .impute_tax() ## impute the NA tax

    #temp_name <- double_tax_name(temp_name,sep=';')

    temp  %<>% tibble::column_to_rownames('feature') %>% as.matrix()
    deposit <- SummarizedExperiment(assays=list(counts=temp),
                                    rowData = temp_name)
  }else if(duplicate_feature==TRUE){
    
    temp_name <- temp %>% dplyr::pull(feature) %>% read.table(text = .,sep = sep,blank.lines.skip=F,quote = "",row.names = NULL,header = FALSE) %>%
      dplyr::mutate_if(~ any(. == ""), ~ dplyr::na_if(., ""))

    #total_tax_info <- c('Domain','Kindom','Phylum','Class','Order','Family','Genus','Species','Strain')
    real_tax_info <- total_tax_info[match(start_level, total_tax_info):length(total_tax_info)]
    colnames(temp_name) <- real_tax_info[1:ncol(temp_name)]
    
   ## make the dupicate anotation into NA
   #strings_to_remove3 <- c(
   #  "uncultured", "unculturebacterium", "uncultureorganism", "unidentified", "unculturemarine",
   #  "uncultureLactobacillus", "uncultureprokaryote", "uncultureBacteroidales", "uncultureDesulfovibrionaceae",
   #  "uncultureLachnospiraceae", "uncultureErysipelotrichales", "unculturePrevotella", "uncultureAnaerotruncus",
   #  "uncultureRuminococcus", "uncultureWautersiella", "uncultureCandidatus", "uncultureEubacterium",
   #  "uncultureClostridiales", "unculturecompost", "uncultureActinobacteridae", "uncultureClostridiaceae",
   #  "uncultureOdoribacter", "uncultureCoxiella", "uncultureactinobacterium", "uncultureErysipelotrichaceae",
   #  "uncultureAnaerolineae", "uncultureBacteroidaceae", "uncultureAlphaproteobacteria", "unculturerumen",
   #  "uncultureRhodospirillaceae", "unculturesoil", "unculturespirochete", "uncultureMethanobacteriales"
   #)
   #occurrences <- colSums(sapply(strings_to_remove3, function(pattern) grepl(pattern, unlist(temp_name))))
   #strings_to_remove3_real <- names(occurrences[occurrences >1]) # search the dupicate anotation which appear more than once
   #temp_name <- apply(temp_name, 2, function(x) {
   #  for (string in strings_to_remove3_real) {
   #    x <- gsub(string, NA, x)       
   #  }
   #  return(x)
   #})
   #temp_name <- as.data.frame(temp_name)    


    temp_name <- data.frame(feature = temp$feature,temp_name) %>%
      .impute_tax()  %>% ## impute the NA tax
      dplyr::mutate(feature = paste0('feature_',1:nrow(temp_name)))

    #temp_name <- double_tax_name(temp_name,sep=';')
      
    temp  %<>% 
      dplyr::mutate(feature = paste0('feature_',1:nrow(temp))) %>%
      tibble::column_to_rownames('feature') %>% as.matrix()
    
    deposit <- SummarizedExperiment(assays=list(counts=temp),
                                    rowData = temp_name)
  }else{
    stop("Parameter duplicate_feature must be TRUE or FALSE!")
  }
  
  if (!is.null(assay_name)) {
    assayNames(deposit, 1) <- assay_name
  }
 }
 return(deposit) 
}  
#' Import kegg data into SummariseExperiment
#'
#' @param file A file path. The file should be deposited in txt format.
#' @param data A dataframe.The row must be the feature and the column is the sample.
#' @param type A character string including ko and ec.
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
#' # More examples and tutorial could be found on website: 
#' # https://liubingdong.github.io/EasyMultiProfiler_tutorial/
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

    if (any(duplicated(data[,1]))) {
      stop("Duplicated ko or ec name, please check the data!")
    }

    if (any(is.na(temp[,-1]))) {
      EMP_message("The NA value has been detected in the data and changed into 0!",color = 32,order = 1,show='message')
      temp[,-1][is.na(temp[,-1])] <- 0
    }

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


#' Import table data into SummariseExperiment
#'
#' @param file A file path. The file should be deposited in txt format.
#' @param data A dataframe.The row must be the feature and the column is the sample.
#' @param sampleID a seris of string character to indicate the samples.
#' @param dfmap A dataframe. Indicate the experiment name, sample source and sample tube details in the omics data.
#' @param assay_name A character string. Indicate what kind of result the data belongs to, such as counts, relative abundance, TPM, etc.
#' @param assay A character string. Indicate the experiment name of the import data in the dfmap.
#' @param duplicate_feature A boolean. Whether the feature exist the dupicated name.
#' @importFrom dplyr all_of
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @return SummmariseExperiment object
#' @export
#'
#' @examples
#' # More examples and tutorial could be found on website: 
#' # https://liubingdong.github.io/EasyMultiProfiler_tutorial/
EMP_normal_import <- function(file=NULL,data=NULL,sampleID=NULL,dfmap=NULL,assay_name=NULL,assay=NULL,duplicate_feature=NULL){
  row_data <- assay_data <- obj <- colname <- feature <- NULL
  if (!is.null(data)) {
    data <- data
  }else {
    data <- read.table(file=file,header = T,sep = '\t',quote="")
  } 
  
  if (!is.null(sampleID)) {
    sampleID <- sampleID
  }else{
    if (!is.null(dfmap) & !is.null(assay) ) {
      sampleID <- dfmap %>% as.data.frame() %>% dplyr::filter(assay == !!assay) %>%
        dplyr::pull(colname)
    }else{
      EMP_message("The function will consider all column as samples.\nIf the data contain rowdata, please define sampleID!",color = 32,order = 1,show='message')
      sampleID <- colnames(data)[-1]
    }
  }
  
  if (is.null(duplicate_feature)){
    check_duplicated_feature <- any(duplicated(data[,1]))
    if(check_duplicated_feature){
      duplicate_feature <- TRUE
      EMP_message("The duplicated name in the feature have been detected.\nParameter duplicate_feature is forced to be TRUE.",color = 32,order = 1,show='message')
    }else{
      duplicate_feature <- FALSE
    }
  }
  
  if (duplicate_feature) {
    if (colnames(data)[1] == 'feature') {
      EMP_message("The original first column has been renamed to '.feature'.",color = 32,order = 1,show='message')
      colnames(data)[1] <- '.feature'
    }
    data <- data |> dplyr::mutate(feature = paste0('feature',1:nrow(data)),.before = 1)
  }
  
  colnames(data)[1] <- 'feature'

  if (any(is.na(data[,sampleID]))) {
    EMP_message("The NA value has been detected in the data and changed into 0!",color = 32,order = 1,show='message')
    data[,sampleID][is.na(data[,sampleID])] <- 0
  }

  data <- data[rowSums(data[,sampleID]) != 0,] # filter away empty feature!
  rownames(data) <- NULL # necessary!
  
  row_data <- data %>% dplyr::select(!all_of(!!sampleID))
  if (ncol(row_data) == 1) {
    row_data <- row_data %>% dplyr::mutate(Name=feature)
  }
  assay_data <- data %>% dplyr::select(feature,all_of(!!sampleID)) %>%
    tibble::column_to_rownames('feature')

  obj <- SummarizedExperiment(assays=list(counts= as.matrix(assay_data)),
                              rowData = row_data)
  if (!is.null(assay_name)) {
    assayNames(obj, 1) <- assay_name
  }
  return(obj)
}

#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom MultiAssayExperiment MultiAssayExperiment
EMP_easy_normal_import <- function(file=NULL,data=NULL,assay='experiment',sampleID=NULL,assay_name=NULL,coldata=NULL,duplicate_feature=NULL,output='MAE') {
  
  SE_object <- EMP_normal_import(file=file,data=data,assay=assay,sampleID=sampleID,assay_name=assay_name,duplicate_feature=duplicate_feature)

  sampleID <-  dimnames(SE_object)[[2]]   
  if (output == 'SE') {
    return(SE_object)
  }else if(output=='MAE'){
   if (is.null(coldata)) {
     stop("If output is MultiAssayExperiment, please input the coldata contaning the group and other information!")
   }
  
   
   if (ncol(coldata) == 0) {
     stop("coldata must contain one column informantion at least!")
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
EMP_easy_taxonomy_import <- function(file=NULL,data=NULL,assay='experiment',assay_name=NULL,coldata=NULL,start_level='Kindom',
                                     humann_format=FALSE,file_format=NULL,duplicate_feature=NULL,output='MAE',
                                     sep=if (humann_format == TRUE) '|' else ';') {
  
  sampleID <- assay_data <- SE_object <- NULL

  SE_object <- EMP_taxonomy_import(file=file,data=data,humann_format=humann_format,file_format=file_format,start_level=start_level,duplicate_feature=duplicate_feature,assay_name=assay_name,sep=sep)                                    
  
  sampleID <-  dimnames(SE_object)[[2]]                                 
  if (output == 'SE') {
    return(SE_object)
  }else if(output=='MAE'){
    if (is.null(coldata)) {
      stop("If output is MultiAssayExperiment, please input the coldata contaning the group and other information!")
    }
    
    if (ncol(coldata) == 0) {
      stop("coldata must contain one column informantion at least!")
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
#' @param sampleID A series of character strings. This parameter helps identify the assay and rowdata only when type = 'normal'.
#' @param file_format A string including biom and gzv.
#' @param humann_format A boolean. Whether the function import the data in the Metaphlan or Humann format.
#' @param duplicate_feature A boolean. Whether the feature exist the dupicated name.
#' @param coldata A dataframe containing one column informantion at least.
#' @param start_level A string. Specific the start level of input data from Domain,Kindom,Phylum,Class,Order,Family,Genus,Species,Strain.(default:Kindom)
#' @param sep The field separator character. Only activated when type =''tax. Values on each line of the file are separated by this character. (defacult:'|') 
#' @param output A character string. Set the output result in SE(SummariseExperiment) or MAE(MultiAssayExperiment) format.
#' @rdname EMP_easy_import
#' @return SummmariseExperiment or MultiAssayExperiment object
#' @export
#' @details Paramter file_format and humann_format help the function import data properly. Data in humann format is usually is generated from Metaphlan and Humann. Data in biom format is usually is generated from Qiime1. Data in qzv format is usually is generated from Qiime2.
#' @examples
#' # More examples and tutorial could be found on website: 
#' # https://liubingdong.github.io/EasyMultiProfiler_tutorial/

EMP_easy_import <- function(file=NULL,data=NULL,type,assay='experiment',assay_name=NULL,sampleID=NULL,coldata=NULL,start_level='Kindom',
                            file_format=NULL,humann_format=FALSE,duplicate_feature=NULL,output='MAE',
                            sep=if (humann_format == "TRUE") '|' else ';'){
  obj <- NULL
  switch(type,
         "tax" = {obj <- EMP_easy_taxonomy_import(file=file,data=data,assay=assay,assay_name=assay_name,coldata=coldata,start_level=start_level,
                                                  humann_format=humann_format,file_format=file_format,duplicate_feature=duplicate_feature,output=output,sep=sep)},
         "ko" = {obj <- EMP_easy_function_import(file=file,data=data,type,assay=assay,assay_name=assay_name,coldata=coldata,
                                                  humann_format=humann_format,output=output)},
         "ec" = {obj <- EMP_easy_function_import(file=file,data=data,type,assay=assay,assay_name=assay_name,coldata=coldata,
                                                  humann_format=humann_format,output=output)},
         "normal" = {obj <- EMP_easy_normal_import(file=file,data=data,assay=assay,sampleID=sampleID,assay_name=assay_name,coldata=coldata,duplicate_feature=duplicate_feature,output=output)},
         {
           stop('type only support tax, ko, ec and normal!')
         }
  )
  return(obj)
}

