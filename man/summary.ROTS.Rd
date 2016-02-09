\name{summary.ROTS}
\alias{summary.ROTS}
\title{ Summary of a ROTS object }
\description{
  Summarizes the differentisl testing resutls from ROTS package.
}
\usage{
  \method{summary}{ROTS}(object,
        fdr = NULL,
        num.genes = NULL,
        verbose = TRUE,
        ...)
}

\arguments{
  \item{object}{a \code{ROTS} object created from differential expression testing run by ROTS.}
  \item{fdr}{ selected cutoff for FDR value.}
  \item{num.genes}{ selected cutoff number for number of differentially detected features. }
  \item{verbose}{ If TRUE (default), summary function will print out 10 first detections fulfill the cutoff criteria. }
  \item{...}{ other arguments passed to the summary function. }
}

\details{
  This fucntion returns the summary information (including row number, test-statistic, pvalue and FDR value)
  for the selected features.}

\value{
Returns a matrix where the rows are the selected features and columns are the Row number, ROTS-statistic, pvalue and FDR.
}

\references{
  L. L. Elo, S. Filen, R. Lahesmaa and T. Aittokallio:
    Reproducibility-optimized test statistic for ranking genes in
  microarray studies. IEEE/ACM Transactions on Computational Biology and
  Bioinformatics 5: 423--431, 2008.
}

\author{ Fatemeh Seyednasrollah, Tomi Suomi, Laura L. Elo
  
  Maintainer: Fatemeh Seyednasrollah <fatsey@utu.fi>
}

\seealso{ \code{\link{affySpikeIn}} }
\examples{
  ## ROTS-statistic for the Affymetrix spike-in data. 
  rots.out <- ROTS(data = affySpikeIn, groups = c(rep(0,5), rep(1,5)),
      B = 100, K = 500 , seed = 1234)
  ## Summary of the ROTS results.
  rots.summary <- summary(rots.out, fdr = 0.05)
}
\keyword{ math }


