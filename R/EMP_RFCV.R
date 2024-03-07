#' Title
#'
#' @param RF wait_for_add
#' @param seed_start wait_for_add
#' @param ntree wait_for_add
#' @param core wait_for_add
#' @param kfold wait_for_add
#' @param rep wait_for_add
#' @param RF_importance wait_for_add
#' @param step wait_for_add
#' @param each_ouput wait_for_add
#' @import ggplot2
#' @importFrom randomForest randomForest
#' @importFrom randomForest varImpPlot
#' @importFrom randomForest importance
#' @importFrom pROC roc
#' @return xx object
#' @export
#'
#' @examples
#' # add example
RFCVSEED <- function(RF,seed_start=123,ntree=1000,core=1,kfold=5,rep=10,RF_importance=1,step=1,each_ouput=F){

  value=c("MeanDecreaseAccuracy","MeanDecreaseGini")
  value=value[RF_importance]

  result=list()
  RF_var=list()
  RF_seed_plot=list()
  # R升级4.0后需要指明因子
  RF$Group=as.factor(RF$Group)
  RFCV_single= function(i){
    # 删除id为i的行，创建训练集
    # 选id为i的行，创建训练集
    trainingset <- subset(data, id !=i)
    testset <- subset(data, id==i)
    #运行一个随机森林模型
    set.seed(seed)
    mymodel <- randomForest::randomForest(trainingset$Group ~ ., data = trainingset[,-ncol(trainingset)], ntree = ntree,keep.forest=T,importance=TRUE)
    #去掉回应列1
    testdata=subset(testset,select =-c(Group,id) )
    temp <- as.data.frame(predict(mymodel, testdata))
    idx=(temp[,1]==testset$Group)
    num=rep(1,length(idx))
    sum(num[idx])/length(idx)
  }

  seed_confirm=c()
  seed_check=seed_start
  while (length(seed_confirm)<rep) {
    set.seed(seed_check)
    num=length(unique(sample(1:kfold, nrow(RF), replace = TRUE)))
    if (num==kfold) {
      seed_confirm=append(seed_confirm,seed_check)
    }
    seed_check=seed_check+1
  }

  for (seed in seed_confirm) {
    set.seed(seed)
    mymodel <- randomForest::randomForest(RF$Group ~ ., data = RF, ntree = ntree,keep.forest=T,importance=TRUE)
    importance=as.data.frame(randomForest::importance(mymodel))
    importance_order=importance[order(importance[,value]),]
    if (each_ouput == T){
    pdf(paste0("seed",seed,"_varimplot.pdf"),width = 20,height = 10)
    randomForest::varImpPlot(mymodel)
    dev.off()
    }
    # 获取倒序剔除特征值顺序
    del_order=row.names(importance_order)
    RF_var[[which(seed_confirm==seed)]]=rev(del_order)


    accuracy_step=list()

    # 决定step_order的，注意不要削减至空数据集 将无法训练
    data_filter_decide<-c( (ncol(RF)-1) /step)==round((ncol(RF)-1)/step)
    if (data_filter_decide&&step!=1) {
      step_order=seq(0,ncol(RF)-3,by=step)
    }else{
      step_order=seq(0,ncol(RF)-2,by=step)
    }


    progress.bar <- create_progress_bar("text")
    progress.bar$init(length(step_order))

    for (j in 1:length(step_order)) {
      # 过滤数据
      if (step_order[j]==0) { data=RF
      }else{
        id_retain=!(colnames(RF)%in%del_order[1:step_order[j]])
        data=RF[id_retain]}
      # 并行开始
      set.seed(seed)
      k=kfold
      data$id <- sample(1:k, nrow(data), replace = TRUE)
      #####
      #cl <- makeCluster(core)
      #registerDoParallel(cl)
      #accuracy=foreach(i=1:k,.packages="randomForest") %dopar% RFCV(i)
      #print(mean(unlist(accuracy)))
      #accuracy_step[[j]]=mean(unlist(accuracy))
      #stopCluster(cl)
      #####
      accuracy=c()
      for (i in 1:2) {
        accuracy=append(accuracy,RFCV_single(i))
      }
      accuracy_step[[j]]= mean(accuracy)
      progress.bar$step()
    }
    result[[which(seed_confirm==seed)]]=accuracy_step


    print(paste0("RF_CV_",seed," is done!"))

  }
  names(result)=paste( "Seed_", seed_confirm, sep="")
  names(RF_var)=paste( "Seed_", seed_confirm, sep="")
  deposit=list()
  deposit$CV_accuracy=result
  deposit$RF_var=RF_var
  deposit$Seed_num=seed_confirm
  deposit$Trails_num=length(seed_confirm)
  deposit$step_order=step_order
  deposit$Feature_num=ncol(RF)-1
  return(deposit)
}





#' Title
#'
#' @param data wait_for_add
#' @param x_break wait_for_add
#' @param cutoff_colour wait_for_add
#' @param palette wait_for_add
#' @return xx object
#' @export
#'
#' @examples
#' # add example
RFCV_plot= function(data,x_break=1,cutoff_colour=c("red"),palette=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF",
                                                                    "#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666")){
  # 确定色板颜色
  if (length(palette)==1) {
    palette=rep(palette,data$Trails_num)
  }
  deposit=list()
  Test=list()
  for (i in 1:data$Trails_num) {
    Test[[i]]=rep(paste0("Trail_",i),length(data$step_order))
  }

  Test=unlist(Test)

  xlab=rep(data$Feature_num-data$step_order,data$Trails_num)
  y=data$CV_accuracy

  for (i in 1:data$Trails_num) {y[[i]]=unlist(y[[i]])
  }

  ylab=unlist(y)
  f=data.frame(xlab,ylab,Test)

  y2=as.data.frame(y)
  total_mean=apply(y2, 1, mean)
  total_sd=apply(y2, 1, sd)
  total_points=1-total_mean+total_sd

  max_mean=max(unlist(total_mean))
  cut_off=data$Feature_num - data$step_order[which(unlist(total_points)%in%min(unlist(total_points)))]
  cut_off_min = min(cut_off)
  opt_mean=mean(as.numeric(y2[cut_off_min,]))

  p=ggplot(f,aes(x=xlab,y=(1-ylab),fill=Test,color=Test))+geom_line()+xlab( "Number of var")+ylab("CV Error")

  curve_plot=p+geom_vline(xintercept = cut_off,colour=cutoff_colour)+scale_x_reverse(breaks=seq(0, data$Feature_num, x_break))+
    labs(title=paste0("Max mean accuracy is ",max_mean," and ",
                      "Optional mean accuracy is ",opt_mean,"  ; ","Optional number is ",cut_off_min))+ theme_bw()+scale_colour_manual(values=palette)

  opt_random_num=colnames(y2)[as.logical(y2[cut_off_min,]==max(y2[cut_off_min,]))]

  # 根据cut_off 确认变量交并集
  var_select=list()
  for (i in names(data$RF_var)) {
    var_select[[i]]=data$RF_var[[i]][1:cut_off_min]
  }

  intersect_num=Reduce(intersect,var_select)
  union_num=Reduce(union,var_select)

  deposit$curve_plot=curve_plot
  deposit$opt_num=cut_off_min
  deposit$opt_random_num=opt_random_num
  deposit$intersect_num=intersect_num
  deposit$union_num=union_num
  return(deposit)
}



#' Title
#'
#' @param name wait_for_add
#' @param ... wait_for_add
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
create_progress_bar <- function(name = "none", ...) {
  if (!is.character(name)) return(name)
  name <- paste("progress", name, sep="_")

  if (!exists(name, mode = "function")) {
    warning("Cannot find progress bar ", name, call. = FALSE)
    progress_none()
  } else {
    match.fun(name)(...)
  }
}
#' Title
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
progress_none <- function() {
  list(
    init = function(x) NULL,
    step = function()  NULL,
    term = function()  NULL
  )
}
#' Title
#'
#' @param style wait_for_add
#' @param ... wait_for_add
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
progress_text <- function(style = 3, ...) {
  n <- 0
  txt <- NULL

  list(
    init = function(x) {
      txt <<- utils::txtProgressBar(max = x, style = style, ...)
      utils::setTxtProgressBar(txt, 0)
    },
    step = function() {
      n <<- n + 1
      utils::setTxtProgressBar(txt, n)
    },
    term = function() close(txt)
  )
}



#' Title
#'
#' @param RF wait_for_add
#' @param seed_start wait_for_add
#' @param ntree wait_for_add
#' @param kfold wait_for_add
#' @param rep wait_for_add
#' @param RF_importance wait_for_add
#' @param step wait_for_add
#' @param each_ouput wait_for_add
#' @param x_break wait_for_add
#' @param cutoff_colour wait_for_add
#' @param palette wait_for_add
#' @importFrom dplyr rename_all
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
EMP_RFCV <- function(RF,seed_start=123,ntree=1000,kfold=5,rep=10,RF_importance=1,step=1,each_ouput=F,
                 x_break=1,cutoff_colour=c("red"),palette=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF",
                                                            "#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666")){
  deposit <- list()
  try(RF <- subset(RF,select=-c(SampleID)),silent=T)
  try(RF %>% dplyr::rename_all(~gsub("\\|", ".", .)))
  RF$Group <- as.factor(RF$Group)
  result <- RFCVSEED(rep = rep,RF,seed_start = seed_start,ntree = ntree,
                   kfold = kfold,RF_importance = RF_importance,step = step,each_ouput=each_ouput)
  result_plot <- RFCV_plot(data = result,x_break = x_break,
                           cutoff_colour = cutoff_colour,palette = palette)
  deposit$RFCV_data <- RF
  deposit$RFCV_result <- result
  deposit$RFCV_result_plot <- result_plot
  return(deposit)
}

#' Title
#'
#' @param x wait_for_add
#' @param rf_estimate_group wait_for_add
#' @param rf_tax_select wait_for_add
#' @param ntree wait_for_add
#' @importFrom pROC roc
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
RFCV_roc <- function(x,rf_estimate_group,rf_tax_select,ntree=1000){
  if (is.list(x)) {
    data <- x$RFCV_data
    seed <- strsplit(x$RFCV_result_plot$opt_random_num,'_')[[1]][2]
    if (rf_tax_select[1]=='union') {
      rf_tax_select <- x$RFCV_result_plot$union_num
    }else if(rf_tax_select[1]=='intersect'){
      rf_tax_select <- x$RFCV_result_plot$intersect_num
    }else{
      rf_tax_select <- rf_tax_select
    }
    data_sub <- data[,c('Group',rf_tax_select)]
  }else if(is.data.frame(x)){
  data_sub <- x[,c('Group',rf_tax_select)]
  }
  set.seed(seed)
  best_model <- randomForest::randomForest(Group ~ ., data = data_sub, ntree = ntree,keep.forest=T,importance=TRUE)
  pred.prob=predict(best_model,type = "prob")
  rocobj1 <- pROC::roc(data_sub$Group, pred.prob[,rf_estimate_group],quiet=T)
  plot(rocobj1, print.auc=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE,
       auc.polygon.col="skyblue", print.thres="best",legacy.axes=T,grid=c(0.1, 0.2),main= paste0(rf_estimate_group))
}

#' Title
#'
#' @param data wait_for_add
#' @param rf_estimate_group wait_for_add
#' @param id_not wait_for_add
#'
#' @return xx object
#' @export
#'
#' @examples
#' # add example
RFCV_data_binary <-function(data,rf_estimate_group,id_not='NOT'){
  new_group <- ifelse(data$Group == rf_estimate_group, rf_estimate_group,ifelse(data$Group != rf_estimate_group, id_not,"NA"))
  data$Group <- as.factor(new_group)
  return(data)
}
