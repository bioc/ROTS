# Function to run the optimization
`optimizeModel` <- function(data, model.original, model.boot, model.null, B, K, seed, BPPARAM) {
  
  # Parameters to test
  a.test <- c(-1, (0:20) / 100, (11:50) / 50, (6:25) / 5)
  k.test <- c( (1:20) * 5, (11:50) * 10, (21:40) * 25, (11:1000) * 100)
  
  if (is.null(K)) {
    K <- floor(nrow(data)/4)
    message(paste("No max top list size K given, using",K))
  }
  if (K>nrow(data)*0.9) {
    warning("Top list size K is more than 90% of the data. Be cautious with reproducibility estimates.")
  }
  K <- min(K,nrow(data))
  k.test <- k.test[k.test < K]
  
  n <- (ncol(model.original)/2)
  names <- gsub("coef\\.","",colnames(model.original)[1:n])
  
  results <- lapply(1:n, function(v) {
    message(paste("Optimizing parameters:",names[v]))
    
    # Test parameters
    ztable <- do.call("rbind", bplapply(a.test, function(a) {
      sapply(k.test, function(k) {
        if(a==-1) {
          d.original <- model.original[,v+(0*n)]
          d.boot <- sapply(model.boot, function(x) x[,v+(0*n)])
          d.null <- sapply(model.null, function(x) x[,v+(0*n)])
        } else {
          d.original <- model.original[,v+(0*n)]/(a+model.original[,v+(1*n)])
          d.boot <- sapply(model.boot, function(x) x[,v+(0*n)]/(a+x[,v+(1*n)]))
          d.null <- sapply(model.null, function(x) x[,v+(0*n)]/(a+x[,v+(1*n)]))
        }
        r.boot <- mean(apply(d.boot, 2, function(x) sum(rank(-abs(d.original))<=k & rank(-abs(x))<=k)/k))
        r.null <- mean(apply(d.null, 2, function(x) sum(rank(-abs(d.original))<=k & rank(-abs(x))<=k)/k))
        sd.boot <- sd(apply(d.boot, 2, function(x) sum(rank(-abs(d.original))<=k & rank(-abs(x))<=k)/k))
        (r.boot-r.null)/sd.boot
      })
    }, BPPARAM=BPPARAM))
    rownames(ztable) <- a.test
    colnames(ztable) <- k.test
    ztable[is.infinite(ztable)] <- NA
    ztable[is.na(ztable)] <- NA
    
    # Produce final results
    a1 <- a.test[which(ztable == max(ztable[is.finite(ztable)]), arr.ind=TRUE)[1]]
    a2 <- 1
    if(a1==-1) {
      a1 <- 1
      a2 <- 0
    }
    k <- k.test[which(ztable == max(ztable[is.finite(ztable)]), arr.ind=TRUE)[2]]
    coef <- model.original[,v]
    d <- model.original[,v] / (a1+a2*model.original[,v+n])
    d.boot <- sapply(model.boot, function(x) x[,v]/(a1+a2*x[,v+n]))
    d.null <- sapply(model.null, function(x) x[,v]/(a1+a2*x[,v+n]))
    r.boot <- mean(apply(d.boot, 2, function(x) sum(rank(-abs(d))<=k & rank(-abs(x))<=k)/k))
    r.null <- mean(apply(d.null, 2, function(x) sum(rank(-abs(d))<=k & rank(-abs(x))<=k)/k))
    sd.boot <- sd(apply(d.boot, 2, function(x) sum(rank(-abs(d))<=k & rank(-abs(x))<=k)/k))
    R <- r.boot
    Z <- (r.boot-r.null)/sd.boot
    pvalue <- unlist(bplapply(d, function(x) {
      1-(sum(abs(d.null)<abs(x),na.rm=TRUE)/(length(d.null)+1))
    }, BPPARAM=BPPARAM))
    FDR <- p.adjust(pvalue, method="BH")
    
    names(d) <- names(coef) <- names(pvalue) <- names(FDR) <- rownames(data)
    out <- list(d=d, coef=coef, pvalue=pvalue, FDR=FDR, a1=a1, a2=a2, k=k, R=R, Z=Z, ztable=ztable)
    return(out)
    
  })
  names(results) <- names
  return(results)
}
