\name{plot.ROTS}
\alias{plot.ROTS}
\title{ Plotting of a ROTS object }
\description{
  Plots the ROTS objects created with the ROTS package.
}
\usage{
  \method{plot}{ROTS}(x, fdr=0.05, type="volcano", ...)
}
\arguments{
  \item{x}{A \code{ROTS} object created from differential expression testing run by ROTS.}
  \item{fdr}{ Selected cutoff for FDR value.}
  \item{type}{ Type of plot to be created. Options are volcano (default) and heatmap. }
  \item{...}{ Other arguments passed to the plot function. }
}

\details{
  This function plots the results from a ROTS object using given false discovery rate threshold.}

\value{
Plots the results from ROTS object.
}

\author{ Fatemeh Seyednasrollah, Tomi Suomi, Laura L. Elo
  Maintainer: Fatemeh Seyednasrollah <fatsey@utu.fi>
}

\seealso{ \code{\link{affySpikeIn}} }
\examples{
  ## ROTS-statistic for the Affymetrix spike-in data. 
  rots.out <- ROTS(data = affySpikeIn, groups = c(rep(0,5), rep(1,5)),
      B = 100, K = 500 , seed = 1234)
  ## Plotting of the ROTS results.
  plot(rots.out)
}
\keyword{ math }


