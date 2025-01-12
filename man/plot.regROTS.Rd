\name{plot.regROTS}
\alias{plot.regROTS}
\title{Plotting of a regROTS object}

\description{
  Plots the regROTS objects created with the ROTS package.
}

\usage{
  \method{plot}{regROTS}(x, fdr=0.05, type=NULL, features=NULL, labels=FALSE, ...)
}

\arguments{
  \item{x}{A \code{regROTS} object created from differential expression testing run by lmeROTS or lmROTS.}
  \item{fdr}{Selected cutoff for FDR value.}
  \item{type}{Type of plot to be created. Options are: 'volcano', 'reproducibility', 'pvalue'.}
  \item{features}{Option to plot only selected variables from the output.}
  \item{labels}{Option to print labels for differentially expressed features.}
  \item{...}{Other arguments passed to the plot function.}
}

\details{
  This function plots the results from a ROTS object using given false discovery rate threshold.
}

\value{
	Plots the results from ROTS object.
}

\author{
	Tomi Suomi, Laura L. Elo\cr
	Maintainer: Tomi Suomi <tomi.suomi@utu.fi>
}

\examples{
  ## Run and plot lmeROTS-statistic 
  \dontrun{rots.out <- lmeROTS("Group + Age + (1|Individual)", data, metadata)}
  \dontrun{plot(rots.out, type="reproducibility")}
}
\keyword{math}


