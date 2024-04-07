#' NetPEA
#'
#' @param genelist genelist.
#' @param network network.
#' @param p restart probability.
#' @param threshold threshold.
#' @param TERM2GENE TERM2GENE.
#' @param TERM2NAME TERM2NAME.
#' @param nperm Number of permutations to do.
#' @param keepAll If FALSE, just keep pathways that have overlap with input genes.
#' @importFrom stats pnorm
#' @importFrom stats sd
#' @export
NetPEA <- function(genelist, network, p = 0.5, TERM2GENE = NULL,
                      TERM2NAME = NULL, threshold = 1e-9, nperm = 100,
                      keepAll = FALSE){
    u <- RWR(genelist, network, p, threshold)
    TERM2GENE <- as.data.frame(TERM2GENE)
    TERM2GENE <- TERM2GENE[!is.na(TERM2GENE[,1]), ]
    TERM2GENE <- TERM2GENE[!is.na(TERM2GENE[,2]), ]
    TERM2GENE <- unique(TERM2GENE)
    PATHID2EXTID <- split(as.character(TERM2GENE[,2]), as.character(TERM2GENE[,1]))

    if (!keepAll) {
        pathLen <- unlist(lapply(PATHID2EXTID, function(x) {length(intersect(x, names(genelist)))}))
        PATHID2EXTID <- PATHID2EXTID[which(pathLen > 0)]
    }

    if ( missing(TERM2NAME) || is.null(TERM2NAME) || all(is.na(TERM2NAME))) {
      # nothing
    } else {
        TERM2NAME <- as.data.frame(TERM2NAME)
        TERM2NAME <- TERM2NAME[!is.na(TERM2NAME[,1]), ]
        TERM2NAME <- TERM2NAME[!is.na(TERM2NAME[,2]), ]
        TERM2NAME <- unique(TERM2NAME)
    }

    distance_scores <- rep(0, length(PATHID2EXTID))
    for (i in seq_len(length(PATHID2EXTID))) {
        uu <- u[PATHID2EXTID[[i]]]
        uu <- uu[!is.na(uu)]
        distance_scores[i] <- mean(1 - uu, na.rm = TRUE)    
        # distance_scores[i] <- mean(1 / uu, na.rm = TRUE)  
    }
    names(distance_scores) <- names(PATHID2EXTID)

    if (nperm == 1) {
        # do nothing
     } else {
        ## random
        distance_scores_random <- matrix(0, nrow = length(PATHID2EXTID), ncol = nperm)
        rownames(distance_scores_random) <- names(distance_scores)
        if (ncol(network) == nrow(network)) {
            netGenes <- rownames(network)
        } else {
            netGenes <- union(network[, 1], network[, 2])
        }
        sampleGene <- genelist[names(genelist) %in% netGenes] 
        for (i in seq_len(nperm)) {
            names(sampleGene) <- sample(netGenes, length(sampleGene))
            u <- quiet(RWR(sampleGene, network, p, threshold))
            for (j in seq_len(length(PATHID2EXTID))) {
                uu <- u[PATHID2EXTID[[j]]]
                uu <- uu[!is.na(uu)]
                distance_scores_random[j, i] <- mean(1 - uu, na.rm = TRUE) 
                # distance_scores_random[j, i] <- mean(1 / uu, na.rm = TRUE)      
            }
        }
        means <- rowMeans(distance_scores_random)
        sds <- apply(distance_scores_random, 1, FUN = sd, na.rm = TRUE)
        zscore <- (distance_scores - means) / sds
        pvalue <- pnorm(zscore, lower.tail=FALSE)
    
    
        result <- data.frame(id = names(distance_scores), 
                             distanceScore = distance_scores,
                             zscore = zscore,
                             pvalue = pvalue)
        if (sum(is.na(result$pvalue)) > 0) {
            result <- result[!is.na(result$pvalue), ]
        }
    
        if (sum(is.na(result$distanceScore)) > 0) {
            result <- result[!is.na(result$distanceScore), ]
        }
    
        result <- result[order(result[, "pvalue"], decreasing = FALSE), ]
        result <- result[!is.na(result[, 2]), ]
        result$name <- TERM2NAME[match(result$id, TERM2NAME[, 1]), 2]
    }

    
    result
}




#' NetPEA2 weighted NetPEA(EnrichNet2)
#'
#' @param genelist genelist.
#' @param network network.
#' @param p restart probability.
#' @param threshold threshold.
#' @param TERM2GENE TERM2GENE.
#' @param TERM2NAME TERM2NAME.
#' @param nperm Number of permutations to do.
#' @param keepAll If FALSE, just keep pathways that have overlap with input genes.
#' @param nperm_method If 1, random genelist, if 2, random genelist and network.
#' @importFrom stats pnorm
#' @importFrom stats quantile
#' @importFrom stats sd
#' @export
NetPEA2 <- function(genelist, network, p = 0.5, TERM2GENE = NULL,
                      TERM2NAME = NULL, threshold = 1e-9, nperm = 100,
                      keepAll = FALSE, nperm_method = 1){
    u <- RWR(genelist, network, p, threshold)
    TERM2GENE <- as.data.frame(TERM2GENE)
    TERM2GENE <- TERM2GENE[!is.na(TERM2GENE[,1]), ]
    TERM2GENE <- TERM2GENE[!is.na(TERM2GENE[,2]), ]
    TERM2GENE <- unique(TERM2GENE)
    PATHID2EXTID <- split(as.character(TERM2GENE[,2]), as.character(TERM2GENE[,1]))

    if (!keepAll) {
        pathLen <- unlist(lapply(PATHID2EXTID, function(x) {length(intersect(x, names(genelist)))}))
        PATHID2EXTID <- PATHID2EXTID[which(pathLen > 0)]
    }

    if ( missing(TERM2NAME) || is.null(TERM2NAME) || all(is.na(TERM2NAME))) {
      # nothing
    } else {
        TERM2NAME <- as.data.frame(TERM2NAME)
        TERM2NAME <- TERM2NAME[!is.na(TERM2NAME[,1]), ]
        TERM2NAME <- TERM2NAME[!is.na(TERM2NAME[,2]), ]
        TERM2NAME <- unique(TERM2NAME)
    }
    ###############################
          # compute gene frequencies accross genesets
    gf = table(TERM2GENE[,2])
 
    if (quantile(gf, 0.99) > mean(gf) + 3 * sd(gf)) {
      gf[gf > quantile(gf, 0.99)] <- quantile(gf, 0.99)
    }
    gff <- function(x) {
      1 + ((max(x) - x)/(max(x) - min(x)))^0.5
    }
    # compute weights
    gf = gff(gf)

    restg = setdiff(names(genelist), names(gf))
    appendd = rep(1, length(restg))
    names(appendd) <- restg
    gf = c(gf, appendd)
  
############################
    distance_scores <- rep(0, length(PATHID2EXTID))
    for (i in seq_len(length(PATHID2EXTID))) {
        uu <- u[PATHID2EXTID[[i]]]
        uu <- uu[!is.na(uu)]
        # weight
        uu_new <- (1 - uu) * gf[names(uu)]
        distance_scores[i] <- mean(uu_new, na.rm = TRUE)    
    }

    

    names(distance_scores) <- names(PATHID2EXTID)

    if (nperm == 1) {
        # do nothing
     } else if (nperm_method == 1) {
        ## random
        distance_scores_random <- matrix(0, nrow = length(PATHID2EXTID), ncol = nperm)
        rownames(distance_scores_random) <- names(distance_scores)
        if (ncol(network) == nrow(network)) {
            netGenes <- rownames(network)
        } else {
            netGenes <- union(network[, 1], network[, 2])
        }
        sampleGene <- genelist[names(genelist) %in% netGenes] 
        for (i in seq_len(nperm)) {
            names(sampleGene) <- sample(netGenes, length(sampleGene))
            u <- quiet(RWR(sampleGene, network, p, threshold))
            for (j in seq_len(length(PATHID2EXTID))) {
                uu <- u[PATHID2EXTID[[j]]]
                uu <- uu[!is.na(uu)]
                # weight
                uu_new <- (1 - uu) * gf[names(uu)]
                distance_scores_random[j, i] <- mean(uu_new, na.rm = TRUE)  
            }
        }
        means <- rowMeans(distance_scores_random)
        sds <- apply(distance_scores_random, 1, FUN = sd, na.rm = TRUE)
        zscore <- (distance_scores - means) / sds
        pvalue <- pnorm(zscore, lower.tail=FALSE)
    
    
        result <- data.frame(id = names(distance_scores), 
                             distanceScore = distance_scores,
                             zscore = zscore,
                             pvalue = pvalue)
        if (sum(is.na(result$pvalue)) > 0) {
            result <- result[!is.na(result$pvalue), ]
        }
    
        if (sum(is.na(result$distanceScore)) > 0) {
            result <- result[!is.na(result$distanceScore), ]
        }
    
        result <- result[order(result[, "pvalue"], decreasing = FALSE), ]
        result <- result[!is.na(result[, 2]), ]
        result$name <- TERM2NAME[match(result$id, TERM2NAME[, 1]), 2]
    } else {
        NULL
    }
 
    return(result)
}
