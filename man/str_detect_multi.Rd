% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Support_EMP_function.R
\name{str_detect_multi}
\alias{str_detect_multi}
\title{Multi-detect the presence/absence of a match}
\usage{
str_detect_multi(string, pattern, exact = FALSE)
}
\arguments{
\item{string}{A character vector.}

\item{pattern}{Pattern to look for.}

\item{exact}{A boolean. Whether search the pattern matched completely.}
}
\value{
Boolean vector
}
\description{
Multi-detect the presence/absence of a match
}
\examples{
text <- c('Bacilli_unclassfiled','Bacteroidia_uncuture','Other')
str_detect_multi(text,c('Bacilli','bacteroidia'),exact=FALSE) # Ignore the capital letter
str_detect_multi(text,c('Bacilli','Bacteroidia'),exact=TRUE) # Set the matched completely
}
\author{
Bingdong Liu
}
