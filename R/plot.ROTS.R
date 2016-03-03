plot.ROTS <- function(x, fdr=0.05, type="volcano", ...) {
  
  # Check for plot type
  if(!(type %in% c("volcano","heatmap"))) {
    stop("Plot type not available")
  }
  
  # Differentially expressed features
  de <- which(x$FDR<fdr)
  
  # Volcano plot
  if(type=="volcano") {
    plot(x$logfc, -log10(x$pvalue), xlab="log2 fold change", ylab="-log10 p-value")
    if(length(de)>0) {
      points(x$logfc[de], -log10(x$pvalue)[de], pch=16, col="darkgreen")
    }
  }
  
  # Heatmap
  if(type=="heatmap" & length(de)>0) {
    try({
      heatmap(as.matrix(x$data[de,]), col=colorRampPalette(c("blue","white","red"))(512))
    },silent=TRUE)
  }
  
}
