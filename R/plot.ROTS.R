plot.ROTS <- function(x, fdr=0.05, ...) {
  
  # Differentially expressed features
  de <- which(x$FDR<fdr)
  
  # Volcano plot
  plot(x$logfc, -log10(x$pvalue), xlab="log2 fold change", ylab="-log10 p-value")
  if(length(de)>0) {
    points(x$logfc[de], -log10(x$pvalue)[de], pch=16, col="darkgreen")
  }
  
  # Heatmap
  try({
    if(length(de)>0) {
      heatmap(as.matrix(x$data[de,]), col=colorRampPalette(c("blue","white","red"))(512))
    }
  },silent=TRUE)
  
}
