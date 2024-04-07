#' wide2long
#'
#' @param network network
#' @param direct if TRUE, build a directed network
#' @importFrom Matrix Matrix
#' @export
wide2long <- function(network, direct = FALSE) {
    netGenes <- union(network[, 1], network[, 2])
    # A <- matrix(0, length(netGenes), length(netGenes))
    A <- Matrix(0, length(netGenes), length(netGenes))
    rownames(A) <- colnames(A) <- netGenes
    for (i in seq_len(nrow(network))) {
        A[network[i, 1], network[i, 2]] <- network[i, 3]
    }

    if (!direct) {
        for (i in seq_len(nrow(network))) {
            A[network[i, 2], network[i, 1]] <- network[i, 3]
        }
    }
    A
}


#' wide2long undirected spMatrix network
#'
#' @param network network
#' @importFrom Matrix spMatrix
#' @export
wide2long_spMatrix <- function(network) {
    netGenes <- union(network[, 1], network[, 2])
    ii <- match(network[, 1], netGenes)
    jj <- match(network[, 2], netGenes)
    i <- c(ii, jj)
    j <- c(jj, ii)
    x <- as.numeric(c(network[, 3], network[, 3]))
    A <- spMatrix(length(netGenes), length(netGenes), i = i, j = j, x = x)
    rownames(A) <- colnames(A) <- netGenes
    A
}


#' Random Walk with Restart
#'
#' @param genelist a vector of -log(pvalue), with the names of gene id.
#' @param network network
#' @param p restart probability
#' @param threshold threshold
#' @export
RWR <- function(genelist, network, p = 0.5, threshold = 1e-9) {
    if (ncol(network) == 2) {
        network[, 3] <- 1
    }

    if (ncol(network) == nrow(network)) {
        A <- network
        netGenes <- rownames(A)
    } else {
        A <- wide2long_spMatrix(network)
        netGenes <- union(network[, 1], network[, 2])
    }
    
    u_old <- v <- rep(0, length(netGenes))
    names(v) <- netGenes
    genes <- intersect(names(v), names(genelist))
    v[genes] <- genelist[genes]
    # v[netGenes %in% genelist] <- 1
    u <- v
    ## colSums(A) may have 0 value.
    colsum <- Matrix::colSums(A)
    colsum[colsum == 0] <- 1
    A <- Matrix::t(A / colsum)
    # A <- A / colSums(A)
    kk <- 1
    while(sum(abs(u-u_old)) > threshold & kk < 100) {
        u_old <- u
        u <- (1-p) * A %*% u_old + p * v
        kk <- kk + 1
    }
    print(kk)
    result <- as.numeric(u)
    names(result) <- netGenes
    result
}


#' NetEnrich
#'
#' @param genelist a vector of gene weights, which names are genes.
#' @param network network
#' @param p restart probability
#' @param threshold threshold
#' @param TERM2GENE TERM2GENE
#' @param TERM2NAME TERM2NAME
#' @param pvalueCutoff Cutoff value of pvalue.
#' @param qvalueCutoff Cutoff of qvalue.
#' @importFrom qvalue qvalue
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", 
#' "bonferroni", "BH", "BY", "fdr", "none"
#' @param n n
#' @param nperm Number of permutations to do.
#' @export
NetEnrich <- function(genelist, network, p = 0, TERM2GENE = NULL,
                      TERM2NAME = NULL, threshold = 1e-9, n = 10, 
                      nperm = 100,  pvalueCutoff = 0.05, 
                      pAdjustMethod = "BH", qvalueCutoff = 0.2){
    if (inherits(genelist, "character")) {
        genename <- genelist
        genelist <- rep(1, length(genename))
        names(genelist) <- genename
    }
    u <- RWR(genelist, network, p, threshold)
    TERM2GENE <- as.data.frame(TERM2GENE)
    TERM2GENE <- TERM2GENE[!is.na(TERM2GENE[,1]), ]
    TERM2GENE <- TERM2GENE[!is.na(TERM2GENE[,2]), ]
    TERM2GENE <- unique(TERM2GENE)
    PATHID2EXTID <- split(as.character(TERM2GENE[,2]), as.character(TERM2GENE[,1]))
    if (missing(TERM2NAME) || is.null(TERM2NAME) || all(is.na(TERM2NAME))) {
      # nothing
    } else {
        TERM2NAME <- as.data.frame(TERM2NAME)
        TERM2NAME <- TERM2NAME[!is.na(TERM2NAME[,1]), ]
        TERM2NAME <- TERM2NAME[!is.na(TERM2NAME[,2]), ]
        TERM2NAME <- unique(TERM2NAME)
    }
    distance_scores <- vector("list", length(PATHID2EXTID))
    for (i in seq_len(length(PATHID2EXTID))) {
        uu <- u[PATHID2EXTID[[i]]]
        uu <- uu[!is.na(uu)]
        distance_scores[[i]] <- sort(1 - uu)  
        # distance_scores[[i]] <- sort(1 / uu)   
    }
    path_scores <- distance_scores
    gene2path_hash <- matrix(0, nrow=length(path_scores), ncol=n)
    for(j in 1:length(path_scores)) {			
		splitscores <- cut(path_scores[[j]], breaks=seq(0.0,1,1/n))	
		freqs <- table(splitscores)
		gene2path_hash[j,] <- as.numeric(freqs)	
	}
    # path2path_hash <- apply(gene2path_hash, 2, sum)/nrow(gene2path_hash)
    path2path_hash <- colMeans(gene2path_hash)
    set_indices <- PATHID2EXTID
    xd_vec <- numeric(length(set_indices))	
	# number of distance bins minus 1
	diam <- n-1	   	
	for(j in 1:length(set_indices)) {
			xd_dist <- xd_distance(gene2path_hash[j,], path2path_hash, diam)
			xd_vec[j] <- xd_dist						
	}

    if (nperm > 1) {
        
        xd_vec_random <- matrix(0, nrow = length(PATHID2EXTID), ncol = nperm)
        rownames(xd_vec_random) <- names(PATHID2EXTID)
        if (ncol(network) == nrow(network)) {
            netGenes <- rownames(network)
        } else {
            netGenes <- union(network[, 1], network[, 2])
        }
        sampleGene <- genelist[names(genelist) %in% netGenes] 
        for (i in seq_len(nperm)) {
            distance_scores <- vector("list", length(PATHID2EXTID))
            names(sampleGene) <- sample(netGenes, length(sampleGene))
            u <- quiet(RWR(sampleGene, network, p, threshold))
            for (j in seq_len(length(PATHID2EXTID))) {
                uu <- u[PATHID2EXTID[[j]]]
                uu <- uu[!is.na(uu)]
                distance_scores[[j]] <- sort(1 - uu)  
                # distance_scores[[j]] <- sort(1 / uu)    
            }
            path_scores <- distance_scores
            gene2path_hash <- matrix(0, nrow=length(path_scores), ncol=n)
            for(j in 1:length(path_scores)) {			
	        	splitscores <- cut(path_scores[[j]], breaks=seq(0.0,1,1/n))	
	        	freqs <- table(splitscores)
	        	gene2path_hash[j,] <- as.numeric(freqs)	
	        }
            path2path_hash <- apply(gene2path_hash, 2, sum)/nrow(gene2path_hash)
            set_indices <- PATHID2EXTID
	        # number of distance bins minus 1
	        diam <- n-1	   	
	        for(j in 1:length(set_indices)) {
	        		xd_dist <- xd_distance(gene2path_hash[j,], path2path_hash, diam)
	        		xd_vec_random[j, i] <- xd_dist						
	        }
        }
        means <- rowMeans(xd_vec_random)
        sds <- apply(xd_vec_random, 1, FUN = sd, na.rm = TRUE)
        zscore <- (xd_vec - means) / sds
        pvalue <- pnorm(zscore, lower.tail=FALSE)
    }

    genes <- names(genelist)

    nodelabels <- names(u)
    # overlap_ids <- sapply(set_indices, function(x) nodelabels[intersect(x, genes)])
    overlap_ids <- sapply(set_indices, function(x) intersect(x, genes))
    overlaps <- sapply(overlap_ids, function(x) length(x))    
    path_lengths <- sapply(set_indices, length)	
    set_names <- names(PATHID2EXTID)
    if (nperm > 1) {
	    resmat <- data.frame(set_names, xd_vec, rep(length(genes),length(set_indices)),
             path_lengths, overlaps, pvalue) 	
	    names(resmat) <- c("path_names", "xd_scores", "upload_sizes", 
            "pathway_sizes", "overlap_sizes", "pvalue")
	    o <- order(pvalue, decreasing=FALSE)
    } else {
        resmat <- data.frame(set_names, xd_vec, rep(length(genes),length(set_indices)),
             path_lengths, overlaps) 	
	    names(resmat) <- c("path_names", "xd_scores", "upload_sizes", 
            "pathway_sizes", "overlap_sizes")
	    o <- order(xd_vec, decreasing=FALSE)
    }

    # optional: output including overlapping gene ids
	# return (list(resmat[o,], overlap_ids[o]))  
	# compact output (no overlapping gene ids)
	
    ##################
    if (!is.null(TERM2NAME)) {
        resmat$Description <- TERM2NAME[match(resmat$path_names, TERM2NAME[, 1]), 2]
    } else {
        resmat$Description <- resmat$path_names
    }
    resmat$adj <- p.adjust(resmat$pvalue, method = pAdjustMethod)
    resmat$qval <- tryCatch(qvalue(resmat$pvalue, lambda=0.05, pi0.method="bootstrap")$qvalue, 
        error=function(e) 0)
    resmat$GeneRatio <- apply(data.frame(a=resmat$overlap_sizes, b=resmat$upload_sizes), 1, function(x)
                       paste(x[1], "/", x[2], sep="", collapse=""))
    nTermGene <- resmat$pathway_sizes
    N <- length(unique(TERM2GENE[, 2]))
    resmat$BgRatio <- apply(data.frame(a=nTermGene, b=N), 1, function(x)
                     paste(x[1], "/", x[2], sep="", collapse=""))   
    resmat$geneID <- vapply(overlap_ids[resmat$path_names], function(i) paste(i, collapse="/"), FUN.VALUE = "1")

    # resmat <- resmat[o,]
    result <- data.frame(ID          = resmat$path_names,
                         Description = resmat$Description,
                         GeneRatio   = resmat$GeneRatio,  
                         BgRatio     = resmat$BgRatio,  
                         pvalue      = resmat$pvalue,
                         p.adjust    = resmat$adj,
                         qvalue      = resmat$qval,
                         geneID      = resmat$geneID,  
                         Count       = resmat$overlap_sizes)   
    
    
    result <- result[result$pvalue < pvalueCutoff, ]
    result <- result[result$qvalue < qvalueCutoff, ]
    result <- result[order(result$pvalue), ]
    background <- TERM2GENE[ ,2] |> as.character() |> unique()
    geneSets <- split(TERM2GENE[, 2], TERM2GENE[, 1])
    x <- new("enrichResult",
             result         = result,
             pvalueCutoff   = pvalueCutoff,
             pAdjustMethod  = pAdjustMethod,
             qvalueCutoff   = qvalueCutoff,
             gene           = names(genelist),
             universe       = background,
             geneSets       = geneSets,
             organism       = "UNKNOWN",
             keytype        = "UNKNOWN",
             ontology       = "UNKNOWN",
             readable       = FALSE
            )
    return(x)
}


#' get_Xd
#'
#' @param distance_score distance_score
#' @param n n
#' @param genelist genelist
#' @param TERM2GENE TERM2GENE
#' @param pia pia
#' @importFrom dplyr ntile
#' @export
get_Xd <- function(distance_score, n, pia, genelist, TERM2GENE) {

    distance <- data.frame(
        distance_score_sort = sort(distance_score),
        ranks = ntile(seq_len(length(distance_score)), n))

    distance_list <- split(rownames(distance), distance[, 2])
    pic <- lapply(distance_list, function(x) {
        length(intersect(x, names(genelist)))}) 
    pic <- unlist(pic)
    if (sum(pic) > 0) {
        pic <- pic / sum(pic) * 100
    }
    sum((pic - pia) / n / seq_len(n))
}


#' get_Xd2
#'
#' @param distance_score distance_score
#' @param n n
#' @param genelist genelist
#' @param u u
#' @importFrom dplyr ntile
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr across
#' @importFrom magrittr %>%
#' @export
get_Xd2 <- function(distance_score, n, genelist, u) {

    distance_u <- data.frame(
        distance = as.numeric(sort(1 - u)),
        ranks = ntile(seq_len(length(u)), n)
    )
    #haha <- split(distance_u[, 1], distance_u[, 2])
    ranks <- distance <- NULL
    distance_min <- distance_u %>% group_by(ranks) %>% summarize(across(distance, min))
    distance_max <- distance_u %>% group_by(ranks) %>% summarize(across(distance, max))
    distance_length <- distance_u %>% group_by(ranks) %>% summarize(across(distance, length))
    
    distance_min_max <- data.frame(min = distance_min[, 2], max = distance_max[, 2], length = distance_max[,2])
    
    distance_pathway <- rep(0, n)
    for (i in 1:n) {
        distance_pathway[i] <- sum(distance_score >= distance_min_max[i, 1] & distance_score < distance_min_max[i, 2])
    }

    pic <- distance_pathway / sum(distance_pathway) * 100 
    pia <- distance_min_max[, 3] / sum(distance_min_max[, 3]) * 100 
    sum((pic - pia) / n / seq_len(n))
}







# from hadley wickham in "https://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html"
#' Suppressing output
#'
#' @param x some code
#' @noRd
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}


xd_distance <- function(distrib1, distrib2, diam) {

			# compute Xd-distance			
			sum_xd <- 0.0
			
			total_real <- sum(distrib1)
			total_rand <- sum(distrib2)
			
			for(xd_iter in 1:(diam+1))
			{
				# prevent division by zero				
				term1 <- distrib1[xd_iter]/total_real
				if(distrib1[xd_iter] == 0)
					term1 <- 0
								
				term2 <- distrib2[xd_iter]/total_rand
				if(distrib2[xd_iter] == 0)
					term2 <- 0
				
				# convert counts to percentages
				sum_xd <- sum_xd + (100*term1 -	100*term2) / (xd_iter * (diam+1))
				
			}
		
			return (sum_xd)
}


#' NetEnrich_rank 
#'
#' Imitate ActivePathways
#' @param genelist genelist
#' @param network network
#' @param p restart probability
#' @param threshold threshold
#' @param TERM2GENE TERM2GENE
#' @param TERM2NAME TERM2NAME
#' @param n n
#' @param nperm Number of permutations to do.
#' @export
NetEnrich_rank <- function(genelist, network, p = 0, TERM2GENE = NULL,
                      TERM2NAME = NULL, threshold = 1e-9, n = 10, nperm = 100){ 
    ranks <- seq(1, length(genelist), 5)                  
    results <- vector("list", length(ranks))
    for (k in 1:length(ranks)) {
        results[[k]] <- NetEnrich(genelist = genelist[1:ranks[k]], 
            network = network, p = p, TERM2GENE = TERM2GENE,
            TERM2NAME = TERM2NAME, threshold = threshold, n = n, nperm = nperm)
        
    }

    pathways <- unique(TERM2GENE[, 1])
    result_mat <- matrix(1, length(ranks), length(pathways))
    colnames(result_mat) <- pathways
    for (k in 1:length(ranks)) {
        result <- results[[k]] 
        result_mat[k, result[, "path_names"]] <- result[, "pvalue"]
    }
    
    resultdf <- data.frame(path_names = pathways, pvalue = rep(1, length(pathways)))
    rownames(resultdf) <- pathways
    for (i in 1:nrow(resultdf)) {
        resultdf[i, 2] <- min(result_mat[, i], na.rm = TRUE)
    }
    return(resultdf)
}


