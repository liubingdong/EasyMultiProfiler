

#' @importFrom Boruta Boruta

.EMP_Boruta_analysis <- function(obj,seed=123,estimate_group,...) {
  assay_data <- coldata <- tran_data <- Boruta_model <- feature_importance <- primary <-  NULL
  assay_data <- obj %>% .get.assay.EMPT()
  coldata <- obj %>% .get.mapping.EMPT() %>%
    dplyr::select(primary,{{estimate_group}}) %>%
    dplyr::rename(Group = {{estimate_group}})
  
  tran_data  <- dplyr::left_join(assay_data,coldata,by='primary') %>% 
    tibble::column_to_rownames('primary')
  tran_data$Group <- factor(tran_data$Group)
  set.seed(seed)
  Boruta_model <- Boruta(Group~.,data=tran_data,...)
  feature_importance <- Boruta_model$finalDecision %>% 
    as.data.frame() %>%
    tibble::rownames_to_column('feature') %>%
    tibble::as_tibble()
  colnames(feature_importance) <- c('feature','Boruta_decision')
  
  obj@deposit[['Boruta_model']] <- Boruta_model
  obj@deposit[['Boruta_feature_importance']] <- feature_importance
  
  .get.estimate_group.EMPT(obj) <- estimate_group
  .get.method.EMPT(obj) <- 'boruta'
  .get.algorithm.EMPT(obj) <- 'boruta'
  .get.info.EMPT(obj) <- 'EMP_marker_analysis'
  return(obj)
}

.EMP_Boruta_analysis_m <- memoise::memoise(.EMP_Boruta_analysis,cache = cachem::cache_mem(max_size = 2048 * 1024^2))

#' @importFrom randomForest randomForest
.EMP_rf_analysis <- function(obj,seed=123,estimate_group,...) {

  assay_data <- coldata <- tran_data <- randomforest_model <- feature_importance <- check_y <- primary <- NULL
  IncNodePurity <- MeanDecreaseAccuracy <- MeanDecreaseGini<- NULL
  assay_data <- obj %>% .get.assay.EMPT()
  coldata <- obj %>% .get.mapping.EMPT() %>%
    dplyr::select(primary,{{estimate_group}}) %>%
    dplyr::rename(group = {{estimate_group}})
  feature_name <- colnames(assay_data)[-1]
  tran_data  <- dplyr::left_join(assay_data,coldata,by='primary') %>% 
    tibble::column_to_rownames('primary')
  
  
  if (class(tran_data$group) == 'numeric' | class(tran_data$group) == 'integer') {
    check_y <- 'regression'
  }else if (class(tran_data$group) == 'character' | class(tran_data$group) == 'factor') {
    check_y <- 'classify'
    tran_data$group <- factor(tran_data$group)
  }
  
  rlang::check_installed(c('BiocManager'), reason = 'for .EMP_rf_analysis().', action = install.packages) 
  rlang::check_installed(c('janitor'), reason = 'for .EMP_rf_analysis().', action = BiocManager::install)

  #colnames(tran_data) <- c(paste0('V',1:(ncol(tran_data)-1)),'Group')
  tran_data <- janitor::clean_names(tran_data) # clean bad names for rf, or error.
  set.seed(seed)
  randomforest_model <- randomForest(group~.,data=tran_data,importance=TRUE,...)
  if (check_y == 'classify') {
    feature_importance <- randomforest_model$importance %>% as.data.frame() %>%
      dplyr::select(MeanDecreaseAccuracy,MeanDecreaseGini) %>%
      dplyr::rename(rf_MDA = MeanDecreaseAccuracy,rf_MDG = MeanDecreaseGini) %>%
      tibble::rownames_to_column('feature') %>%
      dplyr::mutate(feature = {{feature_name}}) %>% tibble::as_tibble()
  }else if (check_y == 'regression') {
    feature_importance <- randomforest_model$importance %>% as.data.frame() %>%
      dplyr::select('%IncMSE','IncNodePurity') %>%
      dplyr::rename(rf_IncMSE = '%IncMSE',rf_IncNodePurity = IncNodePurity) %>%
      tibble::rownames_to_column('feature') %>%
      dplyr::mutate(feature = {{feature_name}}) %>% tibble::as_tibble()
  }else {
    stop('check_y is wrong!')
  }

  obj@deposit[['rf_model']] <- randomforest_model
  obj@deposit[['rf_feature_importance']] <- feature_importance
  
  .get.estimate_group.EMPT(obj) <- estimate_group
  .get.method.EMPT(obj) <- 'randomforest'
  .get.algorithm.EMPT(obj) <- check_y
  .get.info.EMPT(obj) <- 'EMP_marker_analysis'
  return(obj)
}

.EMP_rf_analysis_m <- memoise::memoise(.EMP_rf_analysis,cache = cachem::cache_mem(max_size = 2048 * 1024^2))


#' @importFrom xgboost xgb.DMatrix
#' @importFrom xgboost xgboost
#' @importFrom xgboost xgb.importance

.EMP_xgb_analysis <- function(obj,seed=123,estimate_group,max.depth=6,eta=0.3,nrounds=50,objective=NULL,xgboost_run=NULL,verbose=0,...) {
  assay_data <- row_data <- coldata <- tran_data <- xgb_model <- feature_importance <- primary <- feature <- NULL
  Feature <- Gain <- Cover <- Frequency <-  Importance <- xgb_Gain <- xgb_Cover <- xgb_Frequency <- xgb_Importance <- nthread <- NULL
  assay_data <- obj %>% 
    .get.assay.EMPT() %>% 
    tibble::column_to_rownames('primary') %>% 
    as.matrix()
  row_data <- obj %>% .get.row_info.EMPT() %>% dplyr::select(1:2)


  if (is.null(objective) | is.null(xgboost_run)) {
    stop('Parameter xgboost_run need specify classify or regression and select the suitable parameter objective!')
  }

  if (xgboost_run == 'classify') {
    coldata <- obj %>% .get.mapping.EMPT() %>%
      dplyr::pull({{estimate_group}}) %>%
      as.factor() %>% 
      as.numeric() %>% -1
  }else if (xgboost_run == 'regression') {
    coldata <- obj %>% .get.mapping.EMPT() %>%
      dplyr::pull({{estimate_group}})
  }else {
    stop('Parameter xgboost_run must be classify or regression!')
  }

  nthread <- parallel::detectCores() - 1

  traindata <- xgb.DMatrix(data = as.matrix(assay_data),label= coldata)
  set.seed(seed)
  xbg_model <- xgboost(data = traindata, 
                       max.depth = max.depth, 
                       eta = eta, 
                       nthread = nthread, 
                       nrounds = nrounds, 
                       objective = objective,
                       verbose=verbose,...)
  feature_importance <- xgb.importance(colnames(traindata), model = xbg_model)
  xgboost::xgb.plot.importance(feature_importance,plot=F) ## This step is necessary to add importance in the feature_importance
  feature_importance <- feature_importance %>%
    dplyr::rename(feature=Feature,xgb_Gain=Gain,xgb_Cover=Cover,xgb_Frequency=Frequency,xgb_Importance=Importance)
  feature_importance <- dplyr::left_join(row_data,feature_importance,by='feature') %>%
    dplyr::select(feature,xgb_Gain,xgb_Cover,xgb_Frequency,xgb_Importance) %>%
    dplyr::mutate_all(~ifelse(is.na(.), 0, .)) 
  
  obj@deposit[['xgb_model']] <- xbg_model
  obj@deposit[['xgb_feature_importance']] <- feature_importance
  
  .get.estimate_group.EMPT(obj) <- estimate_group
  .get.method.EMPT(obj) <- 'xgboost'
  .get.algorithm.EMPT(obj) <- xgboost_run
  .get.info.EMPT(obj) <- 'EMP_marker_analysis'
  return(obj)
}

.EMP_xgb_analysis_m <- memoise::memoise(.EMP_xgb_analysis,cache = cachem::cache_mem(max_size = 2048 * 1024^2))


#' @importFrom stats coef

.EMP_lasso_analysis <- function(obj,estimate_group,seed=123,nfolds=5,lambda_select='lambda.min',...) {
  rlang::check_installed(c('BiocManager'), reason = 'for .EMP_lasso_analysis().', action = install.packages) 
  rlang::check_installed(c('glmnet'), reason = 'for .EMP_lasso_analysis().', action = BiocManager::install)
  assay_data <- coldata <- tran_data <- lasso_model <- feature_importance <- primary <- NULL
  
  assay_data <- obj %>% 
    .get.assay.EMPT() %>% 
    tibble::column_to_rownames('primary') %>% 
    as.matrix()
  
  coldata <- obj %>% .get.mapping.EMPT() %>%
    dplyr::select(primary,{{estimate_group}}) %>%
    dplyr::rename(Group = {{estimate_group}}) %>% 
    tibble::column_to_rownames('primary') %>% 
    as.matrix()
  set.seed(seed)
  lasso_model <- glmnet::cv.glmnet(assay_data,coldata,type.measure='mse', family = "gaussian",nfolds = nfolds,alpha = 1,...)
  
  feature_importance <- coef(lasso_model,s=lasso_model[[lambda_select]]) %>% as.matrix()
  feature_importance <- feature_importance[-1,] %>% 
    as.data.frame() %>%
    tibble::rownames_to_column('feature') %>% tibble::as_tibble()
  colnames(feature_importance) <- c('feature','lasso_coe')
  
  obj@deposit[['lasso_model']] <- lasso_model
  obj@deposit[['lasso_feature_importance']] <- feature_importance
  
  .get.estimate_group.EMPT(obj) <- estimate_group
  .get.method.EMPT(obj) <- 'lasso'
  .get.algorithm.EMPT(obj) <- 'regression'
  .get.info.EMPT(obj) <- 'EMP_marker_analysis'
  return(obj)
}

.EMP_lasso_analysis_m <- memoise::memoise(.EMP_lasso_analysis,cache = cachem::cache_mem(max_size = 2048 * 1024^2))


#' Marker discover based on classify or regression model
#'
#' @param obj Object in EMPT or MultiAssayExperiment format.
#' @param experiment A character string. Experiment name in the MultiAssayExperiment object.
#' @param method A character string. Method includs boruta, randomforest, xgboost, lasso.
#' @param estimate_group A character string. Select the group name in the coldata to be calculated.
#' @param seed An interger. Set the random seed to the plot.
#' @param nfolds An interger. Only actived when method = 'lasso'. More imformation in glmnet::cv.glmnet.
#' @param lambda_select A character string including lambda.min or lambda.1se. Only actived when method = 'lasso'. More imformation in glmnet::cv.glmnet.
#' @param max.depth An interger (default:6). Only actived when method = 'xgboost'. More imformation in xgboost::xgboost.
#' @param eta A number (0.3). Only actived when method = 'xgboost'. More imformation in xgboost::xgboost.
#' @param nrounds An interger (default:50). Only actived when method = 'xgboost'. More imformation in xgboost::xgboost.
#' @param xgboost_run An character string.Parameter xgboost_run need specify classify or regression and select the suitable parameter objective. More imformation in xgboost::xgboost.
#' @param objective An character string. Only actived when method = 'xgboost'. More imformation in xgboost::xgboost. eg. binary:logistic for two categories classify,multi:softmax for multible categories classify and reg:squarederror for linear regression.                                 
#' @param verbose An interger (default:0). Only actived when method = 'xgboost'. More imformation in xgboost::xgboost.
#' @param use_cached A boolean. Whether the function use the results in cache or re-compute.
#' @param action A character string. Whether to join the new information to the EMPT (add), or just get the detailed result generated here (get).
#' @param ... Further parameters passed to the function Boruta::Boruta, randomForest, xgboost::xgboost, glmnet::cv.glmnet.
#' @importFrom memoise forget
#' @return EMPT object
#' @export
#'
#' @examples
#' data(MAE)
#' ## To estimate the improtance of feature by Boruta algorithm
#' MAE |>
#'   EMP_marker_analysis(experiment = 'geno_ec',method = 'boruta',
#'                       estimate_group = 'Group') |>
#'   EMP_filter(feature_condition = Boruta_decision!= 'Rejected') ## select the Confirmed and Tentative feature
#' 
#' ## regression or classify by randomforest
#' MAE |>
#'   EMP_marker_analysis(experiment = 'geno_ec',method = 'randomforest',
#'                       estimate_group = 'Group') 
#' 
#' MAE |>
#'   EMP_marker_analysis(experiment = 'geno_ec',method = 'randomforest',
#'                       estimate_group = 'Education_Years') 
#' 
#' ## regression or classify by xgboost
#' ### For regression
#' MAE |>
#'   EMP_marker_analysis(experiment = 'geno_ec',method = 'xgboost',xgboost_run='regression',
#'                       estimate_group = 'Education_Years',objective = 'reg:squarederror')
#' ### For two categories classify
#' MAE |>
#'   EMP_marker_analysis(experiment = 'geno_ec',method = 'xgboost',xgboost_run='classify',
#'                       estimate_group = 'Group',objective = 'binary:logistic')
#' ### For multible categories classify
#' MAE |>
#'   EMP_marker_analysis(experiment = 'geno_ec',method = 'xgboost',xgboost_run='classify',
#'                       estimate_group = 'Status',objective = 'multi:softmax',
#'                       num_class=3) ## num_class is necessary
#' ## Lasso regression
#' MAE |>
#'   EMP_marker_analysis(experiment = 'geno_ko',method = 'lasso',estimate_group = 'Education_Years') |>
#'   EMP_filter(feature_condition = lasso_coe >0) # Select the imprortant feature
EMP_marker_analysis <- function(obj,experiment,method,estimate_group=NULL,seed=123,nfolds=5,lambda_select='lambda.min',
                                  max.depth=6,eta=0.3,nrounds=50,xgboost_run=NULL,objective=NULL,verbose=0,use_cached=TRUE,action='add',...){

  call <- match.call()
  if (inherits(obj,"MultiAssayExperiment")) {
    EMPT <- .as.EMPT(obj,
                     experiment = experiment)
    .get.method.EMPT(EMPT) <- method
  }else if(inherits(obj,'EMPT')) {
    EMPT <-obj
    .get.method.EMPT(EMPT) <- method
  }
  
  estimate_group <- .check_estimate_group.EMPT(EMPT,estimate_group)

  if (use_cached == FALSE) {
    memoise::forget(.EMP_Boruta_analysis_m) %>% invisible()
    memoise::forget(.EMP_rf_analysis_m) %>% invisible()
    memoise::forget(.EMP_xgb_analysis_m) %>% invisible()
    memoise::forget(.EMP_lasso_analysis_m) %>% invisible()
  }

  switch(method,
         "boruta" = {EMPT <- .EMP_Boruta_analysis_m(obj=EMPT,seed=seed,estimate_group=estimate_group,...)},
         "randomforest"  = {EMPT <- .EMP_rf_analysis_m(obj=EMPT,seed=seed,estimate_group=estimate_group,...)},
         "xgboost"  = {EMPT <- .EMP_xgb_analysis_m(obj=EMPT,seed=seed,estimate_group=estimate_group,max.depth=max.depth,
                                                      eta=eta,nrounds=nrounds,xgboost_run=xgboost_run,objective=objective,verbose=verbose,...)},
         "lasso"  = {EMPT <- .EMP_lasso_analysis_m(obj=EMPT,estimate_group=estimate_group,seed=seed,nfolds=nfolds,lambda_select=lambda_select,...)},
         )
  class(EMPT) <- 'EMP_marker_analysis'
  .get.history.EMPT(EMPT) <- call
  if (action == 'add') {
    return(EMPT)
  }else if(action == 'get') {
    return(.get.result.EMPT(EMPT))
  }else{
    warning('action should be one of add or get!')
  }
}


