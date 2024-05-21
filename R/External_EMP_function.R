# Due to some unexpected conflicts, we extract some source code from the open-source R package and modifed them to fit EMP package.
# We strongly appreciate their contribution and if any issues, please contact us with Email: 382983280@qq.com


##### correlation analysis code is from agricolae package(https://cran.r-project.org/web/packages/agricolae/index.html)
####' Title
####'
####' @param x wait_for_add
####' @param y wait_for_add
####' @param method wait_for_add
####' @param alternative wait_for_add
####'
####' @return xx object
####' @noRd
####'
####' @examples
####' # add example
###agricolae_correlation <- function (x, y = NULL,
###                                   method = c("pearson","kendall","spearman","lin"),
###                                   alternative = "two.sided")
###{
###x1<-x
###y1<-y
###    method <- match.arg(method)
###    if (is.data.frame(y)) y1 <- as.matrix(y)
###    if (is.data.frame(x)) x1 <- as.matrix(x)
###
###if (!is.null(y1)) {
###   if (!is.matrix(x1) &  !is.matrix(y1)) {
####-------------
#### cor.vector
###name.xy<-paste(deparse(substitute(x)), "and", deparse(substitute(y)))
###xx<-cbind(x,y)
###yy<-na.omit(xx)
###nn<-length(yy[,1])
###x<-yy[,1]
###y<-yy[,2]
###corr<-correl(x,y,method=method,alternative=alternative)
###stat<-corr$stat
###coef<-corr$rho
###pvalue<-corr$pvalue
###if(method=="pearson") {
###gl<-nn-2
###cat("\nPearson's product-moment correlation\n\n")
###cat("data:",name.xy,"\n")
###cat("t =",stat,", df =",gl,", p-value =",pvalue,"\n")
###cat("alternative hypothesis: true rho is ")
###if(alternative == "two.sided" ) cat("not equal to 0")
###if(alternative == "less" ) cat("less than 0")
###if(alternative == "greater") cat("greater than 0")
###cat("\nsample estimates:\ncor\n",coef,"\n")
####list(t=t,df=gl,p.value=pvalue,rho=coef)
###}
###if(method=="spearman"){
###cat("\nSpearman's rank correlation rho\n\n")
###cat("data:",name.xy,"\n")
###cat("p-value =",pvalue,"\n")
###cat("alternative hypothesis: true rho is ")
###if(alternative == "two.sided" ) cat("not equal to 0")
###if(alternative == "less" ) cat("less than 0")
###if(alternative == "greater") cat("greater than 0")
###cat("\nsample estimates:\nrho\n",coef,"\n")
####list(S=t,p.value=pvalue,rho=coef)
###}
###if(method=="lin") {
###gl<-nn-2
###cat("\nLin's concordance index\n\n")
###cat("data:",name.xy,"\n")
###cat("t =",stat,", df =",gl,", p-value =",pvalue,"\n")
###cat("alternative hypothesis: true rho is ")
###if(alternative == "two.sided" ) cat("not equal to 0")
###if(alternative == "less" ) cat("less than 0")
###if(alternative == "greater") cat("greater than 0")
###cat("\nsample estimates:\ncor\n",coef,"\n")
####list(t=t,df=gl,p.value=pvalue,rho=coef)
###}
###if(method=="kendall"){
###cat("\nKendall's rank correlation tau\n\n")
###cat("data:",name.xy,"\n")
###cat("z-norm = ",stat,"p-value =",pvalue,"\n")
###cat("alternative hypothesis: true rho is ")
###if(alternative == "two.sided" ) cat("not equal to 0")
###if(alternative == "less" ) cat("less than 0")
###if(alternative == "greater") cat("greater than 0")
###cat("\nsample estimates:\ntau\n",coef,"\n")
####list(z.norm=stat,p.value=pvalue,tau=coef)
###}
####-------------
###}
###   else {
#### realizar cor.mv
####-------------
###     if (is.data.frame(x) | is.matrix(x)) {
###        nvarx <- ncol(x)
###        if (is.matrix(x)) x<-data.frame(x)
###        nombrex <- names(x)
###        if (is.matrix(x)) nombrex<-colnames(x)
###        x <- as.matrix(x)
###    }
###    else {
###        nvarx <- 1
###        nombrex<- deparse(substitute(x))
###        x <- as.matrix(x)
###    }
###    if (is.data.frame(y) | is.matrix(y)) {
###        nvary <- ncol(y)
###        nombrey <- names(y)
###        if (is.matrix(y)) nombrey<-colnames(y)
###        y <- as.matrix(y)
###    }
###    else {
###        nvary <- 1
###        nombrey<- deparse(substitute(y))
###        y <- as.matrix(y)
###    }
###    estimate <- rep(0, nvarx*nvary)
###    dim(estimate) <- c(nvarx, nvary)
###    dimnames(estimate) <- list(nombrex, nombrey)
###    pvalue <- estimate
###    nn <- round(estimate, 0)
###    for (i in 1:nvarx) {
###        for (j in 1:nvary) {
###            xx <- cbind(x[, i], y[, j])
###            yy <- na.omit(xx)
###            nn[i, j] <- length(yy[, 1])
###            xa <- yy[, 1]
###            yb <- yy[, 2]
###            corr <- correl(xa, yb, method = method, alternative=alternative)
###            estimate[i, j] <- corr$rho
###            pvalue[i, j] <- corr$pvalue
###         }
###    }
###    names(method) = ""
###    estimate <- round(estimate, 2)
###    pvalue <- round(pvalue, 4)
###    n1<-unique(c(nn))
###    if(length(n1)==1)nn<-n1
####cat("\nCorrelation Analysis\n\nMethod     :",method)
####cat("\nAlternative:",alternative,"\n\n")
###    lista <- list(correlation = estimate, pvalue = pvalue,
###        n.obs = nn)
###    return(lista)
####-------------
###}
###}
###else {
####-------------
#### "cor.matrix"
###nvar<-ncol(x)
###estimate<-rep(0,nvar*nvar)
###nombres<-names(x)
###if (is.matrix(x)) nombres<-colnames(x)
###dim(estimate)<-c(nvar,nvar)
###dimnames(estimate)<-list(nombres,nombres)
###pvalue<-estimate
###nn <- round(estimate,0)
###x <-as.matrix(x)
###for(i in 1:nvar){
###for(j in 1:nvar){
###xx<-cbind(x[,i],x[,j])
###yy<-na.omit(xx)
###nn[i,j]<-length(yy[,1])
###x0<-yy[,1]
###y0<-yy[,2]
###if (i==j) {pvalue[i,j]=0 ; estimate[i,j]=1}
###else {
###corr<-correl(x0,y0,method=method,alternative=alternative)
###estimate[i,j]<-corr$rho
###pvalue[i,j]<-corr$pvalue
###}
###}
###}
###names(method)=""
###estimate<-round(estimate,2)
###diag(pvalue)<-1
####cat("\nCorrelation Analysis\n\nMethod     :",method)
####cat("\nAlternative:",alternative,"\n\n")
###n1<-unique(c(nn))
###if(length(n1)==1)nn<-n1
###lista<-list(correlation=estimate,pvalue=pvalue, n.obs=nn)
###return(lista)
####-------------
###}
###}
###vark <- function(x, y){
###  ties.x <- rle(sort(x))$lengths
###  ties.y <- rle(sort(y))$lengths
###  n <- length(x)
###  t1 <- n * (n - 1) * (2 * n + 5)
###  t2 <- sum(ties.x * (ties.x - 1) * (2 * ties.x + 5))
###  t3 <- sum(ties.y * (ties.y - 1) * (2 * ties.y + 5))
###  v1 <- (t1 - t2 - t3)/18
###  t1 <- sum(ties.x * (ties.x - 1) * (ties.x - 2))
###  t2 <- sum(ties.y * (ties.y - 1) * (ties.y - 2))
###  v2 <- (t1 * t2)/(9 * n * (n - 1) * (n - 2))
###  t1 <- sum(ties.x * (ties.x - 1)) * sum(ties.y * (ties.y - 1))
###  v3 <- t1/(2 * n * (n - 1))
###  v1 + v2 + v3
###}
###
### #' @importFrom stats pnorm
###kendall <- function(data1,data2) {
###  n<-length(data1)
###  n2<-0;n1<-0;is<-0
###  n0<-n-1
###  for (j in 1:n0) { 
###    jj<-j+1
###    for (k in jj:n) {
###    a1<-data1[j]-data1[k]
###    a2<-data2[j]-data2[k]
###    aa<-a1*a2
###    if (! is.na(aa)) {
###      if (aa) {
###        n1<-n1+1
###        n2<-n2+1
###        if(aa > 0.0) is<-is+1
###        else is=is-1
###      } else { 
###        if (a1) n1<-n1+1
###        if (a2) n2<-n2+1
###        }
###      }
###    }
###  }
###  tau<-is/(sqrt(n1)*sqrt(n2))
###  #
###  #svar<-(4.0*n+10.0)/(9.0*n*(n-1.0))
###  #z<-tau/sqrt(svar)
###  #prob<-erfcc(abs(z)/1.4142136)
###  z<-is/sqrt(vark(data1,data2))
###  prob<- 2*pnorm(-abs(z))
###  return(list(stat=z,tau=tau,pvalue=prob))
###}


## WGCNA analysis code is from agricolae package(https://cran.r-project.org/web/packages/WGCNA/index.html)
#' Title
#'
#' @param x wait_for_add
#' @param y wait_for_add
#' @param method wait_for_add
#' @param alternative wait_for_add
#'
#' @return xx object
#' @noRd
#'
#' @examples
#' # add example
correl <- function(x,y,method = "pearson", alternative = "two.sided"){
  x<-1.0*x;y<-1.0*y
  n<-length(x)
  if(method=="kendall"){
    # where is kendall?
  corr <- kendall(x,y)
  stat<-corr$stat
  rho<-corr$tau
  if(alternative == "two.sided" ) pvalue<-corr$pvalue
  if(alternative == "less" ) pvalue<-1-corr$pvalue/2
  if(alternative == "greater") pvalue<-corr$pvalue/2
  }
  if(method=="spearman" ){
  a<-rank(x)
  b<-rank(y)
  x<-a
  y<-b
  }
  if ((method =="pearson") | (method=="spearman")) {
  sumx<-sum(x^2)-sum(x)^2/n
  sumy<-sum(y^2)-sum(y)^2/n
  sumxy<-sum(x*y)-sum(x)*sum(y)/n
  rho<-sumxy/sqrt(sumx*sumy)
  gl<-n-2
  stat<-rho*sqrt(gl)/(sqrt(1-rho^2))
  if(alternative == "two.sided" ) pvalue<-2*(1-pt(abs(stat),gl))
  if(alternative == "less" ) pvalue<-pt(abs(stat),gl)
  if(alternative == "greater") pvalue<-1-pt(abs(stat),gl)
  }
  if (method =="lin") {
  mx<-mean(x)
  my<-mean(y)
  sumx<-(sum(x^2)-sum(x)^2/n)/n
  sumy<-(sum(y^2)-sum(y)^2/n)/n
  sumxy<-(sum(x*y)-sum(x)*sum(y)/n)/n
  r<-sumxy/sqrt(sumx*sumy)
  rho<-2*sumxy/(sumx+sumy+(mx-my)^2)
  gl<-n-2
  sdlin<-sqrt((1/gl)*((1-r^2)*rho^2*(1-rho^2)/r^2+2*rho^3*(1-rho)*(mx-my)^2/(r*sqrt(sumx*sumy))-rho^4*(mx-my)^4/(2*sumx*sumy*r^2)))
  stat<-rho/sdlin
  if(alternative == "two.sided" ) pvalue<-2*(1-pt(abs(stat),gl))
  if(alternative == "less" ) pvalue<-pt(abs(stat),gl)
  if(alternative == "greater") pvalue<-1-pt(abs(stat),gl)
  }
  list(stat=stat,rho=rho,pvalue=pvalue)
}

#' @import dynamicTreeCut
#' @importFrom impute impute.knn
WGCNA_blockwiseModules <- function (datExpr, weights = NULL, checkMissingData = TRUE, blocks = NULL, 
    maxBlockSize = 5000, blockSizePenaltyPower = 5, nPreclusteringCenters = as.integer(min(ncol(datExpr)/20, 
        100 * ncol(datExpr)/maxBlockSize)), randomSeed = 54321, 
    loadTOM = FALSE, corType = "pearson", maxPOutliers = 1, quickCor = 0, 
    pearsonFallback = "individual", cosineCorrelation = FALSE, 
    power = 6, networkType = "unsigned", replaceMissingAdjacencies = FALSE, 
    TOMType = "signed", TOMDenom = "min", suppressTOMForZeroAdjacencies = FALSE, 
    suppressNegativeTOM = FALSE, getTOMs = NULL, saveTOMs = FALSE, 
    saveTOMFileBase = "blockwiseTOM", deepSplit = 2, detectCutHeight = 0.995, 
    minModuleSize = min(20, ncol(datExpr)/2), maxCoreScatter = NULL, 
    minGap = NULL, maxAbsCoreScatter = NULL, minAbsGap = NULL, 
    minSplitHeight = NULL, minAbsSplitHeight = NULL, useBranchEigennodeDissim = FALSE, 
    minBranchEigennodeDissim = mergeCutHeight, stabilityLabels = NULL, 
    stabilityCriterion = c("Individual fraction", "Common fraction"), 
    minStabilityDissim = NULL, pamStage = TRUE, pamRespectsDendro = TRUE, 
    reassignThreshold = 1e-06, minCoreKME = 0.5, minCoreKMESize = minModuleSize/3, 
    minKMEtoStay = 0.3, mergeCutHeight = 0.15, impute = TRUE, 
    trapErrors = FALSE, numericLabels = FALSE, nThreads = 0, 
    useInternalMatrixAlgebra = FALSE, useCorOptionsThroughout = TRUE, 
    verbose = 0, indent = 0, ...) 
{
    .corFnc <- getFromNamespace(".corFnc", "WGCNA")
    .TOMDenoms <- getFromNamespace(".TOMDenoms", "WGCNA")
    .TOMTypes <- getFromNamespace(".TOMTypes", "WGCNA")
    .checkAndScaleWeights <- getFromNamespace(".checkAndScaleWeights", "WGCNA")
    .corOptionList <- getFromNamespace(".corOptionList", "WGCNA")
    .corTypes <- getFromNamespace(".corTypes", "WGCNA")
    .networkTypes <- getFromNamespace(".networkTypes", "WGCNA")
    .orderLabelsBySize <- getFromNamespace(".orderLabelsBySize", "WGCNA")
    .pearsonFallbacks <- getFromNamespace(".pearsonFallbacks", "WGCNA")
    .useNThreads <- getFromNamespace(".useNThreads", "WGCNA")
    collectGarbage <- getFromNamespace("collectGarbage", "WGCNA")
    corPvalueFisher <- getFromNamespace("corPvalueFisher", "WGCNA")
    goodSamplesGenes <- getFromNamespace("goodSamplesGenes", "WGCNA")
    mergeCloseModules <- getFromNamespace("mergeCloseModules", "WGCNA")
    moduleEigengenes <- getFromNamespace("moduleEigengenes", "WGCNA")
    plotDendroAndColors <- getFromNamespace("plotDendroAndColors", "WGCNA")
    projectiveKMeans <- getFromNamespace("projectiveKMeans", "WGCNA")
    spaste <- getFromNamespace("spaste", "WGCNA")


    spaces = dynamicTreeCut::indentSpaces(indent)
    if (verbose > 0) 
        printFlush(paste(spaces, "Calculating module eigengenes block-wise from all genes"))
    seedSaved = FALSE
    if (!is.null(randomSeed)) {
        if (exists(".Random.seed")) {
            seedSaved = TRUE
            savedSeed = .Random.seed
        }
        set.seed(randomSeed)
    }
    intCorType = pmatch(corType, .corTypes)
    if (is.na(intCorType)) 
        stop(paste("Invalid 'corType'. Recognized values are", 
            paste(.corTypes, collapse = ", ")))
    intTOMType = pmatch(TOMType, .TOMTypes)
    if (is.na(intTOMType)) 
        stop(paste("Invalid 'TOMType'. Recognized values are", 
            paste(.TOMTypes, collapse = ", ")))
    TOMDenomC = pmatch(TOMDenom, .TOMDenoms) - 1
    if (is.na(TOMDenomC)) 
        stop(paste("Invalid 'TOMDenom'. Recognized values are", 
            paste(.TOMDenoms, collapse = ", ")))
    if ((maxPOutliers < 0) | (maxPOutliers > 1)) 
        stop("maxPOutliers must be between 0 and 1.")
    if (quickCor < 0) 
        stop("quickCor must be positive.")
    if (nThreads < 0) 
        stop("nThreads must be positive.")
    if (is.null(nThreads) || (nThreads == 0)) 
        nThreads = .useNThreads()
    if ((power < 1) | (power > 30)) 
        stop("power must be between 1 and 30.")
    intNetworkType = charmatch(networkType, .networkTypes)
    if (is.na(intNetworkType)) 
        stop(paste("Unrecognized networkType argument.", "Recognized values are (unique abbreviations of)", 
            paste(.networkTypes, collapse = ", ")))
    fallback = pmatch(pearsonFallback, .pearsonFallbacks)
    if (is.na(fallback)) 
        stop(spaste("Unrecognized value '", pearsonFallback, 
            "' of argument 'pearsonFallback'.", "Recognized values are (unique abbreviations of)\n", 
            paste(.pearsonFallbacks, collapse = ", ")))
    datExpr = as.matrix(datExpr)
    dimEx = dim(datExpr)
    if (length(dimEx) != 2) 
        stop("datExpr has incorrect dimensions.")
    weights = .checkAndScaleWeights(weights, datExpr, scaleByMax = FALSE)
    haveWeights = length(weights) > 0
    nGenes = dimEx[2]
    nSamples = dimEx[1]
    allLabels = rep(0, nGenes)
    AllMEs = NULL
    allLabelIndex = NULL
    originalSampleNames = rownames(datExpr)
    if (is.null(originalSampleNames)) 
        originalSampleNames = spaste("Row.", 1:nrow(datExpr))
    originalGeneNames = colnames(datExpr)
    if (is.null(originalGeneNames)) 
        originalGeneNames = spaste("Column.", 1:ncol(datExpr))
    if (!is.null(blocks) && (length(blocks) != nGenes)) 
        stop("Input error: the length of 'geneRank' does not equal the number of genes in given 'datExpr'.")
    if (!is.null(getTOMs)) 
        warning("getTOMs is deprecated, please use saveTOMs instead.")
    if (checkMissingData) {
        gsg = goodSamplesGenes(datExpr, weights = weights, verbose = verbose - 
            1, indent = indent + 1)
        if (!gsg$allOK) {
            datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
            if (haveWeights) 
                weights = weights[gsg$goodSamples, gsg$goodGenes]
        }
        nGGenes = sum(gsg$goodGenes)
        nGSamples = sum(gsg$goodSamples)
    }
    else {
        nGGenes = nGenes
        nGSamples = nSamples
        gsg = list(goodSamples = rep(TRUE, nSamples), goodGenes = rep(TRUE, 
            nGenes), allOK = TRUE)
    }
    if (any(is.na(datExpr))) {
        datExpr.scaled.imputed = t(impute.knn(t(scale(datExpr)))$data)
    }
    else datExpr.scaled.imputed = scale(datExpr)
    # corFnc = .corFnc[intCorType]
    corFnc = WGCNA::cor
    corOptions = .corOptionList[[intCorType]]
    corFncAcceptsWeights = intCorType == 1
    if (useCorOptionsThroughout) 
        corOptions = c(corOptions, list(cosine = cosineCorrelation))
    if (intCorType == 2 && useCorOptionsThroughout) 
        corOptions = c(corOptions, list(maxPOutliers = maxPOutliers, 
            pearsonFallback = pearsonFallback))
    signed = networkType %in% c("signed", "signed hybrid")
    otherArgs = list(...)
    if (useBranchEigennodeDissim) {
        branchSplitFnc = list("branchEigengeneDissim")
        externalSplitOptions = list(list(corFnc = corFnc, corOptions = corOptions, 
            signed = signed))
        nExternalBranchSplitFnc = 1
        externalSplitFncNeedsDistance = FALSE
        minExternalSplit = minBranchEigennodeDissim
    }
    else {
        branchSplitFnc = list()
        externalSplitOptions = list()
        externalSplitFncNeedsDistance = logical(0)
        nExternalBranchSplitFnc = 0
        minExternalSplit = numeric(0)
    }
    if (!is.null(stabilityLabels)) {
        stabilityCriterion = match.arg(stabilityCriterion)
        branchSplitFnc = c(branchSplitFnc, if (stabilityCriterion == 
            "Individual fraction") "branchSplitFromStabilityLabels.individualFraction" else "branchSplitFromStabilityLabels")
        minExternalSplit = c(minExternalSplit, minStabilityDissim)
        externalSplitFncNeedsDistance = c(externalSplitFncNeedsDistance, 
            FALSE)
        print(dim(stabilityLabels))
        externalSplitOptions = c(externalSplitOptions, list(list(stabilityLabels = stabilityLabels)))
    }
    if ("useBranchSplit" %in% names(otherArgs)) {
        if (otherArgs$useBranchSplit) {
            nExternalBranchSplitFnc = nExternalBranchSplitFnc + 
                1
            branchSplitFnc[[nExternalBranchSplitFnc]] = "branchSplit"
            externalSplitOptions[[nExternalBranchSplitFnc]] = list(discardProp = 0.08, 
                minCentralProp = 0.75, nConsideredPCs = 3, signed = signed, 
                getDetails = FALSE)
            externalSplitFncNeedsDistance[nExternalBranchSplitFnc] = FALSE
            minExternalSplit[nExternalBranchSplitFnc] = otherArgs$minBranchSplit
        }
    }
    if (is.null(blocks)) {
        if (nGGenes > maxBlockSize) {
            if (verbose > 1) 
                printFlush(paste(spaces, "....pre-clustering genes to determine blocks.."))
            clustering = projectiveKMeans(datExpr, preferredSize = maxBlockSize, 
                checkData = FALSE, sizePenaltyPower = blockSizePenaltyPower, 
                nCenters = nPreclusteringCenters, randomSeed = randomSeed, 
                verbose = verbose - 2, indent = indent + 1)
            gBlocks = .orderLabelsBySize(clustering$clusters)
            if (verbose > 2) {
                printFlush("Block sizes:")
                print(table(gBlocks))
            }
        }
        else gBlocks = rep(1, nGGenes)
        blocks = rep(NA, nGenes)
        blocks[gsg$goodGenes] = gBlocks
    }
    else {
        gBlocks = blocks[gsg$goodGenes]
    }
    blockLevels = as.numeric(levels(factor(gBlocks)))
    blockSizes = table(gBlocks)
    if (any(blockSizes > sqrt(2^31) - 1)) 
        printFlush(spaste(spaces, "Found block(s) with size(s) larger than limit of 'int' indexing.\n", 
            spaces, " Support for such large blocks is experimental; please report\n", 
            spaces, " any issues to Peter.Langfelder@gmail.com."))
    nBlocks = length(blockLevels)
    dendros = list()
    TOMFiles = rep("", nBlocks)
    blockGenes = list()
    maxUsedLabel = 0
    for (blockNo in 1:nBlocks) {
        if (verbose > 1) 
            printFlush(paste(spaces, "..Working on block", blockNo, 
                "."))
        blockGenes[[blockNo]] = c(1:nGenes)[gsg$goodGenes][gBlocks == 
            blockLevels[blockNo]]
        block = c(1:nGGenes)[gBlocks == blockLevels[blockNo]]
        selExpr = as.matrix(datExpr[, block])
        if (haveWeights) 
            selWeights = weights[, block]
        nBlockGenes = length(block)
        TOMFiles[blockNo] = spaste(saveTOMFileBase, "-block.", 
            blockNo, ".RData")
        if (loadTOM) {
            if (verbose > 2) 
                printFlush(paste(spaces, "  ..loading TOM for block", 
                  blockNo, "from file", TOMFiles[blockNo]))
            x = try(load(file = TOMFiles[blockNo]), silent = TRUE)
            if (x != "TOM") {
                loadTOM = FALSE
                printFlush(spaste("Loading of TOM in block ", 
                  blockNo, " failed:\n file ", TOMFiles[blockNo], 
                  "\n  either does not exist or does not contain the object 'TOM'.\n", 
                  "  Will recalculate TOM."))
            }
            else if (!inherits(TOM, "dist")) {
                printFlush(spaste("TOM file ", TOMFiles[blockNo], 
                  " does not contain object of the right type or size.\n", 
                  " Will recalculate TOM."))
            }
            else {
                size.1 = attr(TOM, "Size")
                if (length(size.1) != 1 || size.1 != nBlockGenes) {
                  printFlush(spaste("TOM file ", TOMFiles[blockNo], 
                    " does not contain object of the right type or size.\n", 
                    " Will recalculate TOM."))
                  loadTOM = FALSE
                }
                else {
                  tom = as.matrix(TOM)
                  rm(TOM)
                  collectGarbage()
                }
            }
        }
        if (!loadTOM) {
            callVerb = max(0, verbose - 1)
            callInd = indent + 2
            CcorType = intCorType - 1
            CnetworkType = intNetworkType - 1
            CTOMType = intTOMType - 1
            warn = as.integer(0)
            tom = .Call("tomSimilarity_call", selExpr, weights, 
                as.integer(CcorType), as.integer(CnetworkType), 
                as.double(power), as.integer(CTOMType), as.integer(TOMDenomC), 
                as.double(maxPOutliers), as.double(quickCor), 
                as.integer(fallback), as.integer(cosineCorrelation), 
                as.integer(replaceMissingAdjacencies), as.integer(suppressTOMForZeroAdjacencies), 
                as.integer(suppressNegativeTOM), as.integer(useInternalMatrixAlgebra), 
                warn, as.integer(nThreads), as.integer(callVerb), 
                as.integer(callInd), PACKAGE = "WGCNA")
            if (saveTOMs) {
                TOM = as.dist(tom)
                TOMFiles[blockNo] = paste(saveTOMFileBase, "-block.", 
                  blockNo, ".RData", sep = "")
                if (verbose > 2) 
                  printFlush(paste(spaces, "  ..saving TOM for block", 
                    blockNo, "into file", TOMFiles[blockNo]))
                save(TOM, file = TOMFiles[blockNo])
                rm(TOM)
                collectGarbage()
            }
        }
        dissTom = 1 - tom
        rm(tom)
        collectGarbage()
        if (verbose > 2) 
            printFlush(paste(spaces, "....clustering.."))
        dendros[[blockNo]] = fastcluster::hclust(as.dist(dissTom), 
            method = "average")
        if (verbose > 2) 
            printFlush(paste(spaces, "....detecting modules.."))
        datExpr.scaled.imputed.block = datExpr.scaled.imputed[, 
            block]
        if (nExternalBranchSplitFnc > 0) 
            for (extBSFnc in 1:nExternalBranchSplitFnc) externalSplitOptions[[extBSFnc]]$expr = datExpr.scaled.imputed.block
        collectGarbage()
        blockLabels = try(cutreeDynamic(dendro = dendros[[blockNo]], 
            deepSplit = deepSplit, cutHeight = detectCutHeight, 
            minClusterSize = minModuleSize, method = "hybrid", 
            distM = dissTom, maxCoreScatter = maxCoreScatter, 
            minGap = minGap, maxAbsCoreScatter = maxAbsCoreScatter, 
            minAbsGap = minAbsGap, minSplitHeight = minSplitHeight, 
            minAbsSplitHeight = minAbsSplitHeight, externalBranchSplitFnc = branchSplitFnc, 
            minExternalSplit = minExternalSplit, externalSplitOptions = externalSplitOptions, 
            externalSplitFncNeedsDistance = externalSplitFncNeedsDistance, 
            assumeSimpleExternalSpecification = FALSE, pamStage = pamStage, 
            pamRespectsDendro = pamRespectsDendro, verbose = verbose - 
                3, indent = indent + 2), silent = FALSE)
        collectGarbage()
        if (verbose > 8) {
            labels0 = blockLabels
            if (interactive()) 
                plotDendroAndColors(dendros[[blockNo]], labels2colors(blockLabels), 
                  dendroLabels = FALSE, main = paste("Block", 
                    blockNo), rowText = blockLabels, textPositions = 1, 
                  rowTextAlignment = "center")
            if (FALSE) 
                plotDendroAndColors(dendros[[blockNo]], labels2colors(allLabels), 
                  dendroLabels = FALSE, main = paste("Block", 
                    blockNo))
        }
        if (inherits(blockLabels, "try-error")) {
            if (verbose > 0) {
                printFlush(paste(spaces, "*** cutreeDynamic returned the following error:\n", 
                  spaces, blockLabels, spaces, "Stopping the module detection here."))
            }
            else warning(paste("blockwiseModules: cutreeDynamic returned the following error:\n", 
                "      ", blockLabels, "---> Continuing with next block. "))
            next
        }
        if (sum(blockLabels > 0) == 0) {
            if (verbose > 1) {
                printFlush(paste(spaces, "No modules detected in block", 
                  blockNo))
            }
            blockNo = blockNo + 1
            next
        }
        blockLabels[blockLabels > 0] = blockLabels[blockLabels > 
            0] + maxUsedLabel
        maxUsedLabel = max(blockLabels)
        if (verbose > 2) 
            printFlush(paste(spaces, "....calculating module eigengenes.."))
        MEs = try(moduleEigengenes(selExpr[, blockLabels != 0], 
            blockLabels[blockLabels != 0], impute = impute, verbose = verbose - 
                3, indent = indent + 2), silent = TRUE)
        if (inherits(MEs, "try-error")) {
            if (trapErrors) {
                if (verbose > 0) {
                  printFlush(paste(spaces, "*** moduleEigengenes failed with the following message:"))
                  printFlush(paste(spaces, "       ", MEs))
                  printFlush(paste(spaces, "    ---> Stopping module detection here."))
                }
                else warning(paste("blockwiseModules: moduleEigengenes failed with the following message:", 
                  "\n     ", MEs, "---> Continuing with next block. "))
                next
            }
            else stop(MEs)
        }
        propMEs = MEs$eigengenes
        blockLabelIndex = as.numeric(substring(names(propMEs), 
            3))
        deleteModules = NULL
        changedModules = NULL
        if (verbose > 2) 
            printFlush(paste(spaces, "....checking kME in modules.."))
        for (mod in 1:ncol(propMEs)) {
            modGenes = (blockLabels == blockLabelIndex[mod])
            KME = do.call(corFnc, c(list(selExpr[, modGenes], 
                propMEs[, mod]), if (corFncAcceptsWeights) list(weights.x = if (haveWeights) weights[, 
                modGenes] else NULL, weights.y = NULL) else NULL, 
                corOptions))
            if (intNetworkType == 1) 
                KME = abs(KME)
            if (sum(KME > minCoreKME) < minCoreKMESize) {
                blockLabels[modGenes] = 0
                deleteModules = c(deleteModules, mod)
                if (verbose > 3) 
                  printFlush(paste(spaces, "    ..deleting module ", 
                    mod, ": of ", sum(modGenes), " total genes in the module\n       only ", 
                    sum(KME > minCoreKME), " have the requisite high correlation with the eigengene.", 
                    sep = ""))
            }
            else if (sum(KME < minKMEtoStay) > 0) {
                if (verbose > 2) 
                  printFlush(paste(spaces, "    ..removing", 
                    sum(KME < minKMEtoStay), "genes from module", 
                    mod, "because their KME is too low."))
                blockLabels[modGenes][KME < minKMEtoStay] = 0
                if (sum(blockLabels[modGenes] > 0) < minModuleSize) {
                  deleteModules = c(deleteModules, mod)
                  blockLabels[modGenes] = 0
                  if (verbose > 3) 
                    printFlush(paste(spaces, "    ..deleting module ", 
                      blockLabelIndex[mod], ": not enough genes in the module after removal of low KME genes.", 
                      sep = ""))
                }
                else {
                  changedModules = union(changedModules, blockLabelIndex[mod])
                }
            }
        }
        if (!is.null(deleteModules)) {
            propMEs = propMEs[, -deleteModules, drop = FALSE]
            modGenes = is.finite(match(blockLabels, blockLabelIndex[deleteModules]))
            blockLabels[modGenes] = 0
            modAllGenes = is.finite(match(allLabels, blockLabelIndex[deleteModules]))
            allLabels[modAllGenes] = 0
            blockLabelIndex = blockLabelIndex[-deleteModules]
        }
        if (sum(blockLabels > 0) == 0) {
            if (verbose > 1) {
                printFlush(paste(spaces, "No significant modules detected in block", 
                  blockNo))
            }
            blockNo = blockNo + 1
            next
        }
        if (is.null(AllMEs)) {
            AllMEs = propMEs
        }
        else AllMEs = cbind(AllMEs, propMEs)
        allLabelIndex = c(allLabelIndex, blockLabelIndex)
        assigned = block[blockLabels != 0]
        allLabels[gsg$goodGenes][assigned] = blockLabels[blockLabels != 
            0]
        rm(dissTom)
        collectGarbage()
        blockNo = blockNo + 1
    }
    deleteModules = NULL
    goodLabels = allLabels[gsg$goodGenes]
    reassignIndex = rep(FALSE, length(goodLabels))
    if (sum(goodLabels != 0) > 0) {
        propLabels = goodLabels[goodLabels != 0]
        assGenes = (c(1:nGenes)[gsg$goodGenes])[goodLabels != 
            0]
        KME = do.call(match.fun(corFnc), c(list(datExpr[, goodLabels != 
            0], AllMEs), if (corFncAcceptsWeights) list(weights.x = if (haveWeights) weights[, 
            goodLabels != 0] else NULL, weights.y = NULL) else NULL, 
            corOptions))
        if (intNetworkType == 1) 
            KME = abs(KME)
        nMods = ncol(AllMEs)
        for (mod in 1:nMods) {
            modGenes = c(1:length(propLabels))[propLabels == 
                allLabelIndex[mod]]
            KMEmodule = KME[modGenes, mod]
            KMEbest = apply(KME[modGenes, , drop = FALSE], 1, 
                max)
            candidates = (KMEmodule < KMEbest)
            candidates[!is.finite(candidates)] = FALSE
            if (FALSE) {
                modDiss = dissTom[goodLabels == allLabelIndex[mod], 
                  goodLabels == allLabelIndex[mod]]
                mod.k = colSums(modDiss)
                boxplot(mod.k ~ candidates)
            }
            if (sum(candidates) > 0) {
                pModule = corPvalueFisher(KMEmodule[candidates], 
                  nSamples)
                whichBest = apply(KME[modGenes[candidates], , 
                  drop = FALSE], 1, which.max)
                pBest = corPvalueFisher(KMEbest[candidates], 
                  nSamples)
                reassign = ifelse(is.finite(pBest/pModule), (pBest/pModule < 
                  reassignThreshold), FALSE)
                if (sum(reassign) > 0) {
                  if (verbose > 2) 
                    printFlush(paste(spaces, " ..reassigning", 
                      sum(reassign), "genes from module", mod, 
                      "to modules with higher KME."))
                  allLabels[assGenes[modGenes[candidates][reassign]]] = whichBest[reassign]
                  changedModules = union(changedModules, whichBest[reassign])
                  if (length(modGenes) - sum(reassign) < minModuleSize) {
                    deleteModules = c(deleteModules, mod)
                  }
                  else changedModules = union(changedModules, 
                    mod)
                }
            }
        }
    }
    if (!is.null(deleteModules)) {
        AllMEs = AllMEs[, -deleteModules, drop = FALSE]
        genes = is.finite(match(allLabels, allLabelIndex[deleteModules]))
        allLabels[genes] = 0
        allLabelIndex = allLabelIndex[-deleteModules]
        goodLabels = allLabels[gsg$goodGenes]
    }
    if (verbose > 1) 
        printFlush(paste(spaces, "..merging modules that are too close.."))
    if (numericLabels) {
        colors = allLabels
    }
    else {
        colors = labels2colors(allLabels)
    }
    mergedAllColors = colors
    MEsOK = TRUE
    mergedMods = try(mergeCloseModules(datExpr, colors[gsg$goodGenes], 
        cutHeight = mergeCutHeight, relabel = TRUE, corFnc = corFnc, 
        corOptions = corOptions, impute = impute, verbose = verbose - 
            2, indent = indent + 2), silent = TRUE)
    if (inherits(mergedMods, "try-error")) {
        warning(paste("blockwiseModules: mergeCloseModules failed with the following error message:\n    ", 
            mergedMods, "\n--> returning unmerged colors.\n"))
        MEs = try(moduleEigengenes(datExpr, colors[gsg$goodGenes], 
            impute = impute, verbose = verbose - 3, indent = indent + 
                3), silent = TRUE)
        if (inherits(MEs, "try-error")) {
            if (!trapErrors) 
                stop(MEs)
            if (verbose > 0) {
                printFlush(paste(spaces, "*** moduleEigengenes failed with the following error message:"))
                printFlush(paste(spaces, "     ", MEs))
                printFlush(paste(spaces, "*** returning no module eigengenes.\n"))
            }
            else warning(paste("blockwiseModules: moduleEigengenes failed with the following error message:\n    ", 
                MEs, "\n--> returning no module eigengenes.\n"))
            allSampleMEs = NULL
            MEsOK = FALSE
        }
        else {
            if (sum(!MEs$validMEs) > 0) {
                mergedAllColors[gsg$goodGenes] = MEs$validColors
                MEs = MEs$eigengenes[, MEs$validMEs]
            }
            else MEs = MEs$eigengenes
            allSampleMEs = as.data.frame(matrix(NA, nrow = nSamples, 
                ncol = ncol(MEs)))
            allSampleMEs[gsg$goodSamples, ] = MEs[, ]
            names(allSampleMEs) = names(MEs)
            rownames(allSampleMEs) = make.unique(originalSampleNames)
        }
    }
    else {
        mergedAllColors[gsg$goodGenes] = mergedMods$colors
        allSampleMEs = as.data.frame(matrix(NA, nrow = nSamples, 
            ncol = ncol(mergedMods$newMEs)))
        allSampleMEs[gsg$goodSamples, ] = mergedMods$newMEs[, 
            ]
        names(allSampleMEs) = names(mergedMods$newMEs)
        rownames(allSampleMEs) = make.unique(originalSampleNames)
    }
    if (seedSaved) 
        .Random.seed <<- savedSeed
    if (!saveTOMs) 
        TOMFiles = NULL
    names(colors) = originalGeneNames
    names(mergedAllColors) = originalGeneNames
    list(colors = mergedAllColors, unmergedColors = colors, MEs = allSampleMEs, 
        goodSamples = gsg$goodSamples, goodGenes = gsg$goodGenes, 
        dendrograms = dendros, TOMFiles = TOMFiles, blockGenes = blockGenes, 
        blocks = blocks, MEsOK = MEsOK)
}

### ## This code 'chectCores' is from the function 'detectCores' in the package 'parallel' for core detect
### chectCores <- function (all.tests = FALSE, logical = TRUE){
###   systems <- list(linux = "grep \"^processor\" /proc/cpuinfo 2>/dev/null | wc -l",
###                   darwin = if (logical) "/usr/sbin/sysctl -n hw.logicalcpu 2>/dev/null" else "/usr/sbin/sysctl -n hw.physicalcpu 2>/dev/null",
###                   solaris = if (logical) "/usr/sbin/psrinfo -v | grep 'Status of.*processor' | wc -l" else "/bin/kstat -p -m cpu_info | grep :core_id | cut -f2 | uniq | wc -l",
###                   freebsd = "/sbin/sysctl -n hw.ncpu 2>/dev/null", openbsd = "/sbin/sysctl -n hw.ncpuonline 2>/dev/null")
###   nm <- names(systems)
###   m <- pmatch(nm, R.version$os)
###   m <- nm[!is.na(m)]
###   if (length(m)) {
###     cmd <- systems[[m]]
###     if (!is.null(a <- tryCatch(suppressWarnings(system(cmd,
###                                                        TRUE)), error = function(e) NULL))) {
###       a <- gsub("^ +", "", a[1])
###       if (grepl("^[1-9]", a))
###         return(as.integer(a))
###     }
###   }
###   if (all.tests) {
###     for (i in seq(systems)) for (cmd in systems[i]) {
###       if (is.null(a <- tryCatch(suppressWarnings(system(cmd,
###                                                         TRUE)), error = function(e) NULL)))
###         next
###       a <- gsub("^ +", "", a[1])
###       if (grepl("^[1-9]", a))
###         return(as.integer(a))
###     }
###   }
###   NA_integer_
### }
