plot.regROTS <- function(x, fdr=0.05, type=NULL, features=NULL, labels=FALSE, ...) {
  
  # Check for plot type
  if (!is.null(type)) {
    if(!(type %in% c("volcano","reproducibility","pvalue"))) {
      stop("Plot type not available. The options are: 'volcano', 'reproducibility', 'pvalue'")
    }
  } else {
    stop("Plot type not selected. The options are: 'volcano', 'reproducibility', 'pvalue'")
  }
  
  # Differentially expressed features
  de <- lapply(x, function(y) which(y$FDR<0.05))
  
  # Which indices to plot
  sel.plot <- seq_along(x) 
  if(!is.null(features)) {
    if(!all(features %in% names(x))) {
      stop(paste("Available features are:", paste(names(x), collapse=", ")))
    }
    sel.plot <- match(features,names(x))
  }
  
  # Volcano plot
  if(type=="volcano") {
    for(s in sel.plot) {
      plot(x[[s]]$coef, -log10(x[[s]]$pvalue), xlab="model coefficient", ylab="-log10 p-value", pch=20, cex=0.5, bty="l")
      mtext(names(x)[s])
      if(length(de[[s]])>0) {
        points(x[[s]]$coef[de[[s]]], -log10(x[[s]]$pvalue[de[[s]]]), pch=21, col="red")
        if (labels==TRUE) {
          par(xpd=TRUE)
          text(x[[s]]$coef[de[[s]]], -log10(x[[s]]$pvalue[de[[s]]]), labels=names(x[[s]]$pvalue[de[[s]]]), pos=3, cex=0.5)
        }
      }
    }
  }
  
  # Reproducibility
  if(type=="reproducibility") {
    for(s in sel.plot) {
      sel <- which(x[[s]]$ztable == max(x[[s]]$ztable[is.finite(x[[s]]$ztable)]), arr.ind=TRUE)
      k <- as.numeric(colnames(x[[s]]$ztable))
      z <- x[[s]]$ztable[sel[1],]
      plot(k, z, pch=20, xlab="Top list size", ylab="Reproducibility Z-score", cex=0.5, panel.first=lines(k, z, col="grey"), bty="l")
      mtext(names(x)[s])
      points(k[which(z==max(z))], z[which(z==max(z))], pch=21, col="red")
      text(k[which(z==max(z))], z[which(z==max(z))], labels=round(max(z),digits=3), pos=4)
      legend("topright",c(paste("a1 =",x[[s]]$a1),paste("a2 =",x[[s]]$a2)), bty="n")
    }
  }
  
  # Histogram of p-values
  if(type=="pvalue") {
    for(s in sel.plot) {
      hist(x[[s]]$pvalue, breaks=20, col="grey", main="", xlab="p-value", ylab="Frequency")
      mtext(names(x)[s])
    }
  }
  
}
