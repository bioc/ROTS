\name{summary.ROTS}
\alias{summary.ROTS}
\title{ Summary of a ROTS object }

\description{
  Summarizes the differential testing results from ROTS package.
}

\usage{
  \method{summary}{ROTS}(object, fdr=NULL, num.genes=NULL, verbose=TRUE, ...)
}

\arguments{
  \item{object}{a \code{ROTS} object created from differential expression testing run by ROTS.}
  \item{fdr}{ selected cutoff for FDR value.}
  \item{num.genes}{ selected cutoff number for number of differentially detected features. }
  \item{verbose}{ If TRUE (default), summary function will print out 10 first detections which fulfill the cutoff criteria. }
  \item{...}{ other arguments passed to the summary function. }
}

\details{
  This function returns the summary information (including row number, test-statistic, pvalue and FDR value) for the selected features.
  }

\value{
	Returns a matrix where the rows are the selected features and columns are the Row number, ROTS-statistic, pvalue and FDR.
}

\references{
	Suomi T, Seyednasrollah F, Jaakkola MK, Faux T, Elo LL.\cr
	ROTS: An R package for reproducibility-optimized statistical testing.\cr
	PLoS Comput Biol 2017; 13: e1005562.
}

\author{
	Fatemeh Seyednasrollah, Tomi Suomi, Laura L. Elo\cr
	Maintainer: Tomi Suomi <tomi.suomi@utu.fi>
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


