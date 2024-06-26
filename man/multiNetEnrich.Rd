% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiEnrichment_MultiNetEnrich.R
\name{multiNetEnrich}
\alias{multiNetEnrich}
\title{multiNetEnrich}
\usage{
multiNetEnrich(
  multiGene,
  network,
  p = 0,
  TERM2GENE = NULL,
  TERM2NAME = NULL,
  threshold = 1e-09,
  n = 10,
  nperm = 100,
  pvalueCutoff = 0.05,
  cutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  combineMethod = "fisher",
  stoufferWeights = NULL,
  output = "enrichResult",
  combineLevel = "enrichResult"
)
}
\arguments{
\item{multiGene}{a data.frame of multi-omics gene difference analysis results (pvalue).
Each row is a gene, and each column represents an omics dataset.}

\item{network}{network}

\item{p}{restart probability}

\item{TERM2GENE}{user input annotation of TERM TO GENE mapping, 
a data.frame of 2 column with term and gene}

\item{TERM2NAME}{user input of TERM TO NAME mapping, 
a data.frame of 2 column with term and name}

\item{threshold}{threshold}

\item{n}{number of bins}

\item{nperm}{Number of permutations to do.}

\item{pvalueCutoff}{Cutoff value of pvalue.}

\item{cutoff}{Pvalue threshold of differentially expressed genes.}

\item{pAdjustMethod}{one of "holm", "hochberg", "hommel", 
"bonferroni", "BH", "BY", "fdr", "none"}

\item{qvalueCutoff}{Cutoff of qvalue.}

\item{combineMethod}{The method of combining pvalues, one of 
"fisher", "edgington", "stouffer" and "Brown"(only used in ActivePathways method).}

\item{stoufferWeights}{weights of stouffer combine method.}

\item{output}{output class, one of "enrichResult", "compareClusterResult" and "list".}

\item{combineLevel}{one of "gene" and "enrichResult"}
}
\description{
multiNetEnrich
}
