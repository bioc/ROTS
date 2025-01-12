\name{summary.regROTS}
\alias{summary.regROTS}
\title{Summary of a regROTS object}

\description{
  Summarizes the differential testing results from ROTS package.
}

\usage{
  \method{summary}{regROTS}(object, fdr=0.05, ...)
}

\arguments{
  \item{object}{A \code{regROTS} object created from differential expression testing run by ROTS.}
  \item{fdr}{Selected cutoff for FDR value.}
  \item{...}{Other arguments passed to the summary function.}
}

\details{
  This function returns the summary information of \code{regROTS} object.
  }

\value{
	Summary of reproducibility-optimisation.
}

\author{
	Tomi Suomi, Laura L. Elo\cr
	Maintainer: Tomi Suomi <tomi.suomi@utu.fi>
}

\examples{
  ## Run and get summary of lmeROTS-statistic 
  \dontrun{rots.out <- lmeROTS("Group + Age + (1|Individual)", data, metadata)}
  \dontrun{summary(rots.out)}
}

\keyword{ math }


