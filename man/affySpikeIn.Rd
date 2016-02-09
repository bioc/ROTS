\name{affySpikeIn}
\alias{affySpikeIn}
\alias{affySpikeIn.gnames}
\alias{affySpikeIn.L}
\alias{spikedgene}

\docType{data}
\title{Gene expression data from the Affymetrix spike-in experiment}
\description{
	The Affymetrix HG-U95Av2 spike-in data contains two groups
	of 12 samples.  In this carefully controlled experiment, it is known
	that 14 genes are differentially expressed between the two sample
	groups, while the rest of the genes are equally expressed. In order to reduce the running time of the vignette example, this package uses a subset of 1000 genes of the original dataset including all the spiked genes with 5 samples from each group. The preprocessed RMA-normalized data were obtained from the Bioconductor
	2.2 package \code{DEDS}.  For more details about the spike-in
	experiment, see http://www.affymetrix.com.
}

\usage{affySpikeIn}
\format{
	\describe{
		\item{\code{affySpikeIn}}{a numeric matrix of gene expression levels with 12626 rows (genes) and 24 columns (samples).}
		\item{\code{affySpikeIn.L}}{a vector indicating the sample groups.}
		\item{\code{affySpikeIn.gnames}}{a character vector containing the names of the genes.}
		\item{\code{spikedgene}}{a numeric vector indicating the locations of the 14 spiked genes.}
	}
}
\keyword{datasets}
