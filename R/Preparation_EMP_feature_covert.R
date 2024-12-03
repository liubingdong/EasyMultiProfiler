
.feature_convert_gene <- function(feature, from = "SYMBOL", to = "ENTREZID", species = "none", OrgDb = NULL) {
  
  # Check if package is installed, otherwise install
  if (find.package("AnnotationDbi", quiet = TRUE) %>% length() == 0) {
    message("EMP_feature_convert need install package AnnotationDbi!")
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install("AnnotationDbi", ask = FALSE)
  }  


  if (!(species %in% c("Human", "Mouse", "Pig", "Zebrafish"))) {
      if (is.null(OrgDb)) {
          stop("The species is not within the built-in species range, OrgDb needs to be provided for conversion.")
      }
      result <- AnnotationDbi::select(x = OrgDb, keys = feature, keytype = from, columns = c(to)) %>% suppressMessages()
      return(result)
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
  result <- result[!duplicated(result[, 1]), ]
  return(result)  
}    


.EMP_feature_convert_gene <- function(EMPT,method='mean',from,to,species = "none", OrgDb = NULL) {
  feature <- NULL
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


.EMP_feature_convert_cpd <- function(EMPT,method='mean',from,to){
  
  feature <- NULL
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

.EMP_feature_convert_tax <- function(EMPT,sep=';',from,add){
  raw_rowdata <- .get.row_info.EMPT(EMPT)
  if ('old_feature' %in% colnames(raw_rowdata)) {
    stop("Convert tax name should perform before EMP_collapse!")
  }else{
    if (from == 'tax_single' & add == 'tax_full') {
      new_rowdata <- raw_rowdata |> .tax_to_full(sep=sep)
    }else if (from == 'tax_full' & add == 'tax_single') {
      new_rowdata <- raw_rowdata |> .tax_to_single(sep=sep)
    }else{
      stop("Parameter from and to must be tax_single or tax_full")
    }
  }
  .get.row_info.EMPT(EMPT) <- new_rowdata
  return(EMPT)
}


.build_disease_data <- function(from,add){
  ref_df <- disease <- feature <- . <- NULL
  switch(add,
             "Human_disease" = {ref_df <- human_disease
                                id_name <- 'doid'
                                disease_name <- 'Human_disease'
             },
             "Mouse_disease" = {ref_df <- mouse_disease
                                id_name <- 'mpid'
                                disease_name <- 'Mouse_disease'
             },
             {
               stop("No proper species for disease is selected!")
             }
  )
  ref_df <- ref_df %>% dplyr::select(!!sym(from),!!sym(id_name),disease) %>%
    dplyr::filter(!is.na(!!sym(from))) %>%
    dplyr::filter(!duplicated(.)) %>%
    dplyr::rename(feature=!!sym(from)) %>%
    dplyr::filter(dplyr::if_all(-feature,~ !is.na(.))) %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(
      !!id_name := paste(unique(.data[[id_name]]), collapse = '/'),
      !!disease_name := paste(unique(disease), collapse = '/'))
  return(ref_df) 
}


.EMP_feature_convert_disease <- function(EMPT,from,add){
  
  raw_rowdata <- .get.row_info.EMPT(EMPT)
  ref_data <- .build_disease_data(from = from,add=add)

  new_rowdata <- dplyr::left_join(raw_rowdata,ref_data,by='feature')
  .get.row_info.EMPT(EMPT) <- new_rowdata
  return(EMPT)
}


#' Covert feature of microbial taxonomy, gene experssion or compund abundance
#'
#' @param obj EMPT or MultiAssayExperiment object.
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object.
#' @param method A character string. Methods include mean, sum, median, min, max. When multiple annotations appear on features, merge activate.
#' @param from A character string. For metabolite include CAS,DTXSID,DTXCID,SID,CID,KEGG,ChEBI,HMDB,Drugbank. For gene include SYMBOL,ENSEMBL,ENTREZID.
#' @param to A character string. For metabolite include CAS,DTXSID,DTXCID,SID,CID,KEGG,ChEBI,HMDB,Drugbank. For gene include SYMBOL,ENSEMBL,ENTREZID.
#' @param add A character string. For microbiome include tax_single and tax_full. For disease include Human_disease and Mouse_disease based on SYMBOL,ENTREZID,ko and ec.
#' @param species A character string. Species includ Human,Mouse,Pig,Zebrafish. If converting feature from other species,please use OrgDb. 
#' @param OrgDb Supported OrgDb listed in 'https://bioconductor.org/packages/release/BiocViews.html#___OrgDb' 
#' @param action A character string. A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @rdname EMP_feature_convert
#' @return EMPT object
#' @export
#'
#' @examples
#' data(MAE)
#' ## for gene ID convert
#' MAE |>
#'   EMP_feature_convert(experiment = 'host_gene',from = 'SYMBOL',to='ENTREZID',species = 'Human')
#' 
#' ## The built-in database only supports Human, Mouse, Pig, Zebrafish
#' ## Other species could utilize OrgDb to convert
#' \dontrun{
#' library(org.Hs.eg.db)
#' MAE |>
#'   EMP_feature_convert(experiment = 'host_gene',from = 'SYMBOL',to='ENTREZID',OrgDb = org.Hs.eg.db)
#' }
#' ## for compound ID convert
#' MAE |>
#'   EMP_collapse(experiment = 'untarget_metabol',na_string=c('NA','null','','-'),
#'                estimate_group = 'MS2kegg',method = 'sum',collapse_by = 'row') |> ## get the compound in KEGG format
#'   EMP_feature_convert(from = 'KEGG',to='HMDB')
#' ## Add full name for microbial data
#' MAE |>
#'   EMP_assay_extract('taxonomy') |>
#'   EMP_feature_convert(from = 'tax_single',add = 'tax_full') |>
#'   EMP_collapse(estimate_group = 'Phylum',collapse_by = 'row')
#' ## Add disease info and select the related feature
#' MAE |> 
#'   EMP_assay_extract('host_gene') |>
#'   EMP_feature_convert(from = 'SYMBOL',add ='Human_disease') |>
#'   EMP_assay_extract(pattern = 'cancer',pattern_ref = 'Human_disease')
EMP_feature_convert <- function(obj,experiment,method='mean',from,to=NULL,add=NULL,species = "none",OrgDb = NULL,action='add'){
  call <- match.call()
  check_result <- NULL
  if (is(obj,"MultiAssayExperiment")) {
    EMPT <- .as.EMPT(obj,
                     experiment = experiment)
  }else if(is(obj,'EMPT')) {
    EMPT <- obj
  }
  
  if (from %in% c('symbol','ensembl','entrezid')) {
    from <- toupper(from)
  }
  if (!is.null(to) && to %in% c('symbol', 'ensembl', 'entrezid')) {
    to <- toupper(to)
  }

  cpd_names_total <- c("CAS", "DTXSID", "DTXCID", "SID", "CID", "KEGG", "ChEBI", "HMDB", "Drugbank")
  gene_names_total <- c('SYMBOL','ENSEMBL','ENTREZID')
  tax_names_total <- c('tax_single','tax_full')
  disease_names_from<- c('SYMBOL','ENTREZID','ko','ec')
  disease_names_add <- c('Human_disease','Mouse_disease')

  
  if (!is.null(to) & !is.null(add)) {
    stop("The parameters 'to' and 'add' cannot both be present simultaneously.")
  }
  
  if (!is.null(to)) {
    
    if (from == to) {
      stop("parameter 'from' should be different from 'to'!")
    }

    if (from %in% gene_names_total & to %in% gene_names_total) {
      EMPT <- EMPT %>% .EMP_feature_convert_gene(method=method,from=from,to=to,species=species,OrgDb=OrgDb)
    }else if(from %in% cpd_names_total & to %in% cpd_names_total){
      EMPT <- EMPT %>%.EMP_feature_convert_cpd(method=method,from=from,to=to)
    }else{
      stop('Pleast check the parameter from and to!')
    }
  }

  if (!is.null(add)) {

    if (from == add) {
      stop("parameter 'from' should be different from 'add'!")
    }

    if(from %in% tax_names_total & add %in% tax_names_total){
      EMPT <- EMPT %>%.EMP_feature_convert_tax(sep=';',from=from,add=add)
    }else if (from %in% disease_names_from & add %in% disease_names_add) {
      EMPT <- EMPT %>%.EMP_feature_convert_disease(from=from,add=add)
    }else{
      stop('Pleast check the parameter from and add!')
    }

    message_info <- list()
    message_info %<>% append(paste0('Feature information has been added: ',add))
    .get.message_info.EMPT(EMPT) <- message_info
  }
  
  # Due to the feature change, result should be removed
  # deposit2 is deposited the inherent information, do not remove
  check_result <- length(EMPT@deposit)!=0 | length(EMPT@deposit_append)!=0 | length(EMPT@plot_deposit)!=0

  if (check_result) {
    message('Due to the feature change, all results should be re-run if needed.')
    EMPT@deposit <- NULL
    EMPT@deposit_append <- NULL
    EMPT@plot_deposit <- NULL
  }

  .get.history.EMPT(EMPT) <- call
  .get.method.EMPT(EMPT) <- 'feature_covert'
  .get.algorithm.EMPT(EMPT) <- 'feature_covert'
  .get.info.EMPT(EMPT) <- 'EMP_assay_data'
  
  if (action=='add') {
    return(EMPT)
  }else if(action=='get'){
    return(.get.assay.EMPT(EMPT))
  }else{
    stop("action should be one of add or get")
  }
  
}

