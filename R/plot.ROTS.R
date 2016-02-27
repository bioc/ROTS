plot.ROTS <- function(obj, fdr=0.05, ...) {
  
  # Differentially expressed features
  de <- which(obj$FDR<fdr)
  
  # Volcano plot
  plot.default(obj$logfc, -log10(obj$pvalue), xlab="log2 fold change", ylab="-log10 p-value")
  if(length(de)>0) {
    points(obj$logfc[de], -log10(rots.out$pvalue)[de], pch=16, col="darkgreen")
  }
  
  # Heatmap
  if(length(de)>0) {
    heatmap(as.matrix(obj$data[de,]), col=colorRampPalette(c("blue","white","red"))(512))
  }
  
}
