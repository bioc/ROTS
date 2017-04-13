plot.ROTS <- function(x, fdr=0.05, type=NULL, labels=FALSE, ...) {
  
  # Check for plot type
  if (!is.null(type)) {
    if(!(type %in% c("volcano","heatmap","ma","reproducibility","pvalue","pca"))) {
      stop("Plot type not available. The options are: 'volcano', 'heatmap', 'ma', 'reproducibility', 'pvalue', 'pca'")
    }
  } else {
    stop("Plot type not selected. The options are: 'volcano', 'heatmap', 'ma', 'reproducibility', 'pvalue', 'pca'")
  }
  
  
  # Differentially expressed features
  de <- which(x$FDR<fdr)
  
  # Volcano plot
  if(type=="volcano") {
    par(xpd=TRUE)
    plot(x$logfc, -log10(x$pvalue), xlab="log2 fold change", ylab="-log10 p-value", pch=20, cex=0.5, bty="l")
    if(length(de)>0) {
      points(x$logfc[de], -log10(x$pvalue)[de], pch=21, col="red")
      if (labels==TRUE) {
        text(x$logfc[de], -log10(x$pvalue)[de], labels=row.names(x$data)[de], pos=3, cex=0.5)
      }
    }
  }
  
  # Heatmap
  if(type=="heatmap" & length(de)>0) {
    try({
      heatmap(as.matrix(x$data[de,]), col=colorRampPalette(c("blue","white","red"))(512), ...)
    },silent=TRUE)
  }
  
  # MA plot
  if(type=="ma") {
    M <- x$logfc
    A <- 0.5*(rowMeans(x$data[,x$cl==1])+rowMeans(x$data[,x$cl==2]))
    par(xpd=TRUE)
    plot(A, M, xlab="A", ylab="M", pch=20, cex=0.5, panel.first=abline(h=0, col="grey", lty=2), bty="l")
    lines(loess.smooth(A, M), col="blue", lty=2)
    if(length(de)>0) {
      points(A[de], M[de], pch=21, col="red")
      if (labels==TRUE) {
        text(A[de], M[de], labels=row.names(x$data)[de], pos=3, cex=0.5)
      }
    }
  }
  
  # Reproducibility
  if(type=="reproducibility") {
	z <- c(0,x$ztable[row.names(x$ztable)==x$a1,])
	k <- c(0,as.numeric(colnames(x$ztable)))
    plot(k, z, pch=20, xlab="Top list size", ylab="Reproducibility Z-score", cex=0.5, panel.first=lines(k, z, col="grey"), bty="l")
    points(k[which(z==max(z))], z[which(z==max(z))], pch=21, col="red")
    text(k[which(z==max(z))], z[which(z==max(z))], labels=round(max(z),digits=3), pos=4)
    legend("bottomleft",c(paste("a1 =",x$a1),paste("a2 =",x$a2)), bty="n", inset=0.1, y.intersp=2)
  }
  
  # Histogram of p-values
  if(type=="pvalue") {
    hist(x$pvalue, breaks=20, col="grey", main="", xlab="p-value", ylab="Frequency")
  }
  
  # PCA
  if(type=="pca") {
    if(length(de)>0) {
      pca <- prcomp(t(x$data[de,]), center=TRUE, scale.=TRUE)
      par(xpd=TRUE)
      plot(pca$x[,1], pca$x[,2], xlab="Principal component 1", ylab="Principal component 2", pch=20, cex=2, col= x$cl, bty="l")
      if (labels==TRUE) {
        text(pca$x[,1], pca$x[,2], labels=colnames(x$data), pos=3, cex=1)
      }
    }
  }
  
}
