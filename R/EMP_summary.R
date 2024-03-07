.getDim <- getFromNamespace(".getDim", "MultiAssayExperiment")

#' Title
#'
#' @param MAE wait_for_add
#' @importFrom MultiAssayExperiment experiments
#' @noRd
.creat_MAE_summary <- function(MAE) {
  object <- MultiAssayExperiment::experiments(MAE)
  o_class <- class(object)
  elem_cl <- vapply(object, function(o) { class(o)[[1L]] }, character(1L))
  o_len <- length(object)
  o_names <- names(object)
  o_dim <- lapply(object, dim)

  featdim <- .getDim(o_dim, 1L)
  sampdim <- .getDim(o_dim, 2L)

  assay_name_list <- c()
  sample_info_num <- c()
  feature_info_num <- c()
  assay_miss <- c()
  coldata_miss <- c()
  rowdata_miss <- c()
  for (i in o_names) {
    assay_name_list %<>% append(names(MAE@ExperimentList@listData[[i]]@assays@data@listData))
    sample_info_num %<>% append(MAE %>% EMP_coldata_extract(i) %>% ncol-1)
    feature_info_num %<>% append(rowData(MAE[[i]]) %>% ncol-1)
    assay_miss %<>% append(assay(MAE[[i]])  %>% is.na() %>% any())
    coldata_miss %<>% append(MAE %>% EMP_coldata_extract(i) %>% is.na() %>% any())
    rowdata_miss %<>% append(rowData(MAE[[i]])  %>% is.na() %>% any())
  }


  dt <- data.frame(Sample = sampdim,
                    Feature = featdim,
                    Assay = assay_name_list,
                    Sample_atrr. = sample_info_num,
                    Feature_atrr. = feature_info_num,
                    Assay_status = assay_miss,
                    Sample_status= coldata_miss,
                    Feature_status = rowdata_miss)

  return(dt)

}


.creat_EMP_summary <- function(EMP) {
  featdim <- c()
  sampdim <- c()
  experiment_name <- names(EMP@ExperimentList)
  for (i in experiment_name) {
    feat_temp <-dim(EMP@ExperimentList[[i]])[1]
    names(feat_temp) <- i
    featdim <- append(featdim,feat_temp)

    samp_temp <- dim(EMP@ExperimentList[[i]])[2]
    names(samp_temp) <- i
    sampdim <- append(sampdim,samp_temp)
  }

  assay_name_list <- c()
  sample_info_num <- c()
  feature_info_num <- c()
  assay_miss <- c()
  coldata_miss <- c()
  rowdata_miss <- c()
  for (i in experiment_name) {
    assay_name_list %<>% append(EMP@ExperimentList[[i]] %>% .get.assay_name.EMPT())
    sample_info_num %<>% append(EMP@ExperimentList[[i]] %>% EMP_coldata_extract() %>% ncol-1)
    feature_info_num %<>% append(rowData(EMP@ExperimentList[[i]]) %>% ncol-1)
    assay_miss %<>% append(assay(EMP@ExperimentList[[i]])  %>% is.na() %>% any())
    coldata_miss %<>% append(EMP@ExperimentList[[i]] %>% EMP_coldata_extract() %>% is.na() %>% any())
    rowdata_miss %<>% append(rowData(EMP@ExperimentList[[i]])  %>% is.na() %>% any())
  }

  dt <- data.frame(Sample = sampdim,
                   Feature = featdim,
                   Assay = assay_name_list,
                   Sample_atrr. = sample_info_num,
                   Feature_atrr. = feature_info_num,
                   Assay_status = assay_miss,
                   Sample_status= coldata_miss,
                   Feature_status = rowdata_miss)
  return(dt)
}



#' Title
#'
#' @param obj wait_for_add
#' @importFrom formattable proportion_bar
#' @importFrom kableExtra cell_spec
#' @importFrom kableExtra kable_styling
#' @importFrom kableExtra column_spec
#' @importFrom kableExtra add_header_above
#' @importFrom kableExtra scroll_box
#' @importFrom kableExtra footnote
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_summary <- function(obj) {
  if (inherits(obj,"MultiAssayExperiment")) {
    dt <- .creat_MAE_summary(obj)
  }else if (inherits(obj,"EMP")) {
    dt <- .creat_EMP_summary(obj)
  }else{
    stop("EMP_summary only support MultiAssayExperiment and EMP class!")
  }

  dt %>%
    dplyr::mutate(
      Sample = formattable::proportion_bar("#4DBBD5FF")(Sample),
      Feature = formattable::proportion_bar("#CC79A7")(Feature),
      #Sample_atrr. = color_tile("white", "#CDDC39")(Sample_atrr.),
      #Feature_atrr. = color_tile("white", "#FF9800")(Feature_atrr.),
      Sample_atrr. = formattable::proportion_bar("#CDDC39")(Sample_atrr.),
      Feature_atrr. = formattable::proportion_bar("#FF9800")(Feature_atrr.),
      Assay = ifelse(Assay == 'counts',
                     yes=kableExtra::cell_spec(Assay, "html", color = "green", bold = T),
                     no=ifelse(Assay == 'relative',kableExtra::cell_spec(Assay, "html", color = "#E69F00", bold = T),
                                                                      kableExtra::cell_spec(Assay, "html", color = "red", bold = T))),
      Assay_status = ifelse(Assay_status == 'FALSE',
                            kableExtra::cell_spec(Assay_status, "html", background = "green", color = "white", bold = T),
                            kableExtra::cell_spec(Assay_status, "html",  background = "red", color = "white", bold = T)),
      Sample_status = ifelse(Sample_status == 'FALSE',
                             kableExtra::cell_spec(Sample_status, "html", background = "green", color = "white",align='justify', bold = T),
                             kableExtra::cell_spec(Sample_status, "html",  background = "red", color = "white",align='justify', bold = T)),
      Feature_status = ifelse(Feature_status == 'FALSE',
                              kableExtra::cell_spec(Feature_status, "html",  background = "green", color = "white",align='justify', bold = T),
                              kableExtra::cell_spec(Feature_status, "html",  background = "red", color = "white",align='justify', bold = T))



    ) %>%
    kableExtra::kable("html", escape = F,caption = "EasyMultiPlot Summary") %>%
    kableExtra::kable_styling("hover", full_width = F,position = 'center') %>%
    kableExtra::column_spec(9, width = "10cm") %>%
    kableExtra::column_spec(4, width = "3cm")  %>%
    kableExtra::column_spec(1, bold = T, color = "white", background = "#3C5488FF") %>%
    kableExtra::add_header_above(c(" ", "Data dimension" = 2, "Data information" = 3,"Data miss" = 3)) %>%
    kableExtra::scroll_box(width = "100%", box_css = '
    padding: 15px; border: 15px solid transparent;
    background: linear-gradient(white,white), repeating-linear-gradient(45deg, #0072b2, #d9230f 10px, #f96352 10px, #f96352 20px);
    background-clip: padding-box, border-box;') %>%
    kableExtra::footnote(general = "Typing for more information: help (EasyMultiPlot)",
             general_title = "Note: ", number_title = "",
             alphabet_title = "Type II: ", symbol_title = "Type III: ",
             footnote_as_chunk = T, title_format = c("italic", "underline")
    )
}
