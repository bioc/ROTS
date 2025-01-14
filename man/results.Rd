\name{results}
\alias{results}
\title{Gathering results from a regROTS object}

\description{
  Gathers the results from regROTS objects created with the ROTS package.
}

\usage{
  results(object, order=FALSE, pvalue=NULL, FDR=NULL, coef=NULL, logfc=NULL)
}

\arguments{
  \item{object}{A \code{regROTS} object created from differential expression testing run by lmeROTS or lmROTS.}
  \item{order}{Ordering of results based on significance (default FALSE).}
  \item{pvalue}{Optional p-value cutoff for filtering the results.}
  \item{FDR}{Optional false discovery rate cutoff for filtering the results.}
  \item{coef}{Optional model coefficient cutoff for filtering the results (\code{regROTS}).}
  \item{logfc}{Optional log fold change cutoff for filtering the results (\code{ROTS}).}
}

\details{
  This is a convenience function to gather results from regROTS objects.
}

\value{
	A list of reproducibility-optimized results for each variable in the model.
}

\author{
	Tomi Suomi, Laura L. Elo\cr
	Maintainer: Tomi Suomi <tomi.suomi@utu.fi>
}

\examples{
  ## Run and plot lmeROTS-statistic 
  \dontrun{rots.out <- lmeROTS("Group + Age + (1|Individual)", data, metadata)}
  \dontrun{results <- results(rots.out)}
}
\keyword{math}
