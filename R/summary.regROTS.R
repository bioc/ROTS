`summary.regROTS` <-
function(object, fdr=0.05, ...) {
  cat("Reproducibility statistics:\n")
  print(do.call("rbind", lapply(object, function(y) unlist(y[c("a1","a2","k","R","Z")]))))
  
  n <- sapply(object, function(y) sum(y$FDR<fdr))
  cat(paste("\nA total of ", sum(n), " rows satisfy condition FDR < ", fdr, "\n", sep=""))
  for(s in seq_along(n)){
    cat(paste(names(n)[s], ": ", n[s], "\n", sep=""))
  }
  
  for(s in seq_along(n)){
    cat(paste("\n",names(n)[s],"\n",sep=""))
    res <- data.frame(coef=object[[s]]$coef, d=object[[s]]$d, p=object[[s]]$pvalue, fdr=object[[s]]$FDR)
    res <- res[order(res$p, -abs(res$d), -abs(res$coef), decreasing=FALSE),]
    sel <- which(res$fdr<fdr)
    if(length(sel)==0) {
      cat(paste("No rows satisfy condition FDR < ", fdr, "\n", sep=""))
    } else {
      if (length(sel)>10) {
        sel <- sel[1:10]
        print(res[sel,])
        cat("(truncated to 10 rows)\n")
      } else {
        print(res[sel,])
      }
    }
  }
}
