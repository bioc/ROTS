\name{lmROTS}
\alias{lmROTS}
\title{Reproducibility-Optimized Test Statistic (ROTS)}

\description{
	Calculates the reproducibility-optimized test statistic (ROTS) for
	ranking genes in order of evidence for differential expression using
	linear model
}

\usage{
	lmROTS(formula, data, metadata, B = 100, K = NULL, seed = NULL,
	BPPARAM = bpparam())
}

\arguments{
  \item{formula}{A string representing a formula that follows the \code{lm}
    syntax and where the dependent variable (e.g. expression ~ ...) is omitted.
    \code{lmROTS} uses each row of \code{data} as the dependent variable
    and looks for explanatory variables from \code{metadata}.}
	\item{data}{A numeric data matrix in which rows correspond
		to genes and columns correspond to samples.}
	\item{metadata}{A numeric data matrix in which rows correspond
		to samples and columns correspond to explanatory variables.}
	\item{B}{An integer specifying the number of bootstrap and permutation
		resamplings (default 100).}
	\item{K}{An integer indicating the largest top list size considered.
	  If no value is given, 1/4 of the features are used.}
	\item{seed}{An integer seed for the random number generator.}
	\item{BPPARAM}{BiocParallel back-end if other than default.}
	
}

\details{The reproducibility-optimization procedure ROTS enables the
	selection of a suitable gene ranking statistic directly from the given
	dataset. This function allows the use of linear models as a basis
	for reproducibility-optimization procedure.
}

\value{
  \code{lmROTS} returns an object of class \code{regROTS}, which is a list
  of ROTS results for each of the explanatory variables in the formula. They
  contain the following components:
	
	\item{d}{the value of the optimized ROTS-statistic for each gene.}
  \item{coef}{coefficients of the model for each gene.}
	\item{pvalue}{the corresponding pvalues.}
	\item{FDR}{the corresponding FDR-values.}
	\item{a1}{the optimized parameter a1.}
	\item{a2}{the optimized parameter a2.}
	\item{k}{the optimized top list size.}
	\item{R}{the optimized reproducibility value.}
	\item{Z}{the optimized reproducibility Z-score.}
	\item{ztable}{the full reproducibility table.}
}

\author{
	Tomi Suomi, Laura L. Elo\cr
	Maintainer: Tomi Suomi <tomi.suomi@utu.fi>
}

\examples{
## lmROTS-statistic
rots.out <- lmeROTS("Group + Age", data, metadata)
}

\keyword{math}
