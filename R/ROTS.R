`ROTS` <-
  function(data, groups, B=1000, K=NULL, paired=FALSE, seed=NULL, a1=NULL, a2=NULL, log=TRUE, progress=FALSE) {
    if (is(data, "ExpressionSet"))
           data <- Biobase::exprs(data)  
    ## Set random number generator seed for reproducibility
    if(!is.null(seed))
      set.seed(seed, kind="default")
    
    ## Check groups
    if(length(groups)!=ncol(data)) {
      stop(paste("Number of samples in the data does not match the groups."))
    }
	
    ## If rownames(data) == NULL, use integers as rownames. Summary function requires that
    ## the data matrix has rownames.
    if(is.null(rownames(data)))
      rownames(data) <- 1:nrow(data)
    
    ##  The reproducibility values in a dense lattice
    ssq <- c( (0:20) / 100, (11:50) / 50, (6:25) / 5)
    N <- c( (1:20) * 5, (11:50) * 10, (21:40) * 25, (11:1000) * 100)
    
    ## Check for top list size
    if (is.null(K)) {
      K <- floor(nrow(data)/4)
      message(paste("No top list size K given, using",K))
    }
    
    ## The top list size cannot be larger than the total number of genes
    K <- min(K,nrow(data))
    N <- N[N < K]
    
    ## Reorder the data according to the group labels and check for NA rows.
    data <- data[,order(groups)]
    groups <- sort(groups)
    
    for (i in unique(groups)) {
      if(any(rowSums(is.na(data[,which(groups==i)])) >= length(which(groups==i))-1)) {
        stop(paste("The data matrix of group",i,"contains rows with less than two non-missing values, please remove these rows."))
      }
    }
    
    cl <- groups+(1-min(groups))
	
    ## Check number of samples for paired test
    if (paired) {
      for (i in unique(cl)[-1]) {
        if(length(which(cl==1))!=length(which(cl==i))) stop("Uneven number of samples for paired test.")
      }
    }
	
    ## Calculate fold change
    if (length(unique(cl))==2) {
      if(log) {
        logfc <- rowMeans(data[,which(cl==1)], na.rm=TRUE) - rowMeans(data[,which(cl==2)], na.rm=TRUE)
      } else {
        logfc <- rowMeans(log2(data[,which(cl==1)]+1), na.rm=TRUE) - rowMeans(log2(data[,which(cl==2)]+1), na.rm=TRUE)
      }
    } else {
      logfc <- rep(NA,nrow(data))
    }
    
    ## ---------------------------------------------------------------------------
 
    ## Bootstrap samples
    message("Bootstrapping samples")
    samples <- bootstrapSamples(data, 2*B, cl, paired)
    ## Permutated samples
    pSamples <- permutatedSamples(data, nrow(samples), cl)
    
    ## Test statistics in the bootstrap (permutate) datasets.
    ## For each bootstrap (permutated) dataset, determine the signal log-ratio
    ## (D) and the standard error (S).
    D <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    S <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    pD <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    pS <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    
    if (progress) pb <- txtProgressBar(min=0, max=nrow(samples), style=3)
    for (i in seq_len(nrow(samples))) {
      samples.R <- split(samples[i,], cl)
      pSamples.R <- split(pSamples[i,], cl)
      
      ## If a1 and a2 parameters are given, we don't need the bootstrap
      ## dataset
      if( is.null(a1) | is.null(a2) ){
        fit <- testStatistic(paired, lapply(samples.R, function(x) data[,x]))
        D[,i] <- fit$d
        S[,i] <- fit$s
      }
      
      pFit <- testStatistic(paired, lapply(pSamples.R, function(x) data[,x]))
      pD[,i] <- pFit$d
      pS[,i] <- pFit$s
      
      if (progress) setTxtProgressBar(pb, i)
    }
    if (progress) close(pb)
    
    ## Free up memory
    rm(samples, pSamples)
    gc()
    
    ## ---------------------------------------------------------------------------
    
    if( is.null(a1) | is.null(a2) ){
      ## Optimize the parameters
      message("Optimizing parameters")
      
      ## Calculate the reproducibility values for all the given a1-values and top
      ## list sizes in both bootstrap and permuted data and their standard
      ## deviations in the bootstrap case
      
      ## Reproducibility matrix (bootstrap data):
      ## the rows correspond to the a1-values given in ssq (+ 1 for signal
      ## log-ratio only: a1=1, a2=0), the columns correspond to the different top
      ## list sizes
      reprotable <- matrix(nrow=length(ssq) + 1, ncol=length(N))
      colnames(reprotable) <- N
      row.names(reprotable) <- c(ssq, "slr")
      
      ## Reproducibility matrix (permuted data):
      ## the rows correspond to the a1-values given in ssq (+ 1 for signal
      ## log-ratio only: a1=1, a2=0), the columns correspond to the different top
      ## list sizes
      reprotable.P <- matrix(nrow=length(ssq) + 1, ncol=length(N))
      colnames(reprotable.P) <- N
      row.names(reprotable.P) <- c(ssq, "slr")
      
      ## Standard deviation matrix for the reproducibility values (bootstrap
      ## data): the rows correspond to the a1-values given in ssq (+ 1 for signal
      ## log-ratio only: a1=1, a2=0), the columns correspond to the different
      ## top list sizes
      reprotable.sd <- matrix(nrow=length(ssq) + 1, ncol=length(N))
      colnames(reprotable.sd) <- N
      row.names(reprotable.sd) <- c(ssq, "slr")
      
      if (progress) pb <- txtProgressBar(min=0, max=length(ssq), style=3)
      for(i in 1 : length(ssq)){
        ## The overlaps between bootstrap samples. Rows correspond to different
        ## bootrsrap samples and columns corrospond to different top list size.
        ## Repeat for each parameter combination a1 and a2.
        overlaps <- matrix(0, nrow=B, ncol=length(N))
        overlaps.P <- matrix(0, nrow=B, ncol=length(N))
        
        ## Call the custom c++-loop 1.
        cResults =  calculateOverlaps1 (D, S, pD, pS, nrow(D), as.integer(N), length(N),
                          ssq[i], as.integer(B), overlaps, overlaps.P)
        
        ## Colmeans & rowMeans are a lot faster than apply
        #         reprotable[i, ] <- colMeans(overlaps)
        #         reprotable.P[i, ] <- colMeans(overlaps.P)
        
        reprotable[i, ] <- colMeans(cResults[["overlaps"]])
        reprotable.P[i, ] <- colMeans(cResults[["overlaps_P"]])
        
        ## Standard deviation for each column
        ## same as reprotable.sd[i, ] <- apply(overlaps, 2, sd)
        ## or just sd(overlaps), but a lot faster.
        #         reprotable.sd[i,] <- sqrt(rowSums((t(overlaps) - reprotable[i,])^2) /
        #                                     (nrow(overlaps) - 1))
        reprotable.sd[i,] <- sqrt(rowSums((t(cResults[["overlaps"]]) - reprotable[i,])^2) /
                                    (nrow(cResults[["overlaps"]]) - 1)) 
        
        if (progress) setTxtProgressBar(pb, i)
      }
      if (progress) close(pb)
      
      i <- length(ssq) + 1
      
      overlaps <- matrix(0, nrow=B, ncol=length(N))
      overlaps.P <- matrix(0, nrow=B, ncol=length(N))
      
      ## Call the custom c++-loop 2.
      cResults = calculateOverlaps2 (D, pD, nrow(D), as.integer(N), length(N),
                      as.integer(B), overlaps, overlaps.P)
      
      ## Free up memory
      rm(D, S)
      gc()
      
      reprotable[i, ] <- colMeans(cResults[["overlaps"]])
      reprotable.P[i, ] <- colMeans(cResults[["overlaps_P"]])
      ## Standard deviation for each column
      reprotable.sd[i,] <- sqrt(rowSums((t(cResults[["overlaps"]]) - reprotable[i,])^2) /
                                  (nrow(cResults[["overlaps"]]) - 1))
      
      ## Free up memory
      rm(overlaps, overlaps.P , cResults)
      gc()
      
      ## -------------------------------------------------------------------------
      
      ## Calculate the Z-statistic and find the top list size and the
      ## (a1,a2)-combination giving the maximal Z-value
      ztable <- (reprotable - reprotable.P) / reprotable.sd
      ## Rownames of ztable are c(ssq, "slr") and colnames are N
      
      ## Free up memory
      rm(reprotable.P, reprotable.sd)
      gc()
      
	  sel <- which(ztable == max(ztable[is.finite(ztable)]), arr.ind=TRUE)
      ## Sel is a matrix containing the location(s) of the largest value (row,
      ## col). If the location of the largest value is not unique then nrow(sel)
      ## > 2 (length(sel) > 2)
      
      if(length(sel)>2) sel <- sel[1,]
      
      if(sel[1] < nrow(reprotable)) {
        a1 <- as.numeric(row.names(reprotable)[sel[1]])
        a2 <- 1
      }
      if(sel[1] == nrow(reprotable)) {
        a1 <- 1
        a2 <- 0
      }
      k <- as.numeric(colnames(reprotable)[sel[2]])
      R <- reprotable[sel[1],sel[2]]
      Z <- ztable[sel[1],sel[2]]
      
      ## Free up memory
      rm(reprotable)
      gc()
      
      ## Calculate the reproducibility-optimized test statistic based on the
      ## reproducibility-maximizing a1, a2 and k values and the corresponding FDR
      fit <- testStatistic(paired, lapply(split(1:length(cl), cl), function(x) data[,x]))
      d <- fit$d / (a1 + a2 * fit$s)
      pD <- pD/(a1 + a2 * pS)
      
      ## Free up memory
      rm(pS)
      gc()
      
      message("Calculating p-values")
      p <- calculateP(d, pD)
      
      message("Calculating FDR")
      FDR <- calculateFDR(d, pD, progress)
      
      ## Free up memory
      rm(pD)
      gc()
      
      ROTS.output <- list(
        data = data,     # Original data
        B = B,           # Number of resamplings
        d = d,           # ROTS-statistic
        logfc = logfc,   # Fold change
        pvalue = p,	     # p-value
        FDR = FDR,       # FDR
        a1 = a1,         # a1
        a2 = a2,         # a2
        k = k,           # top list size
        R = R,	         # reproducibility value
        Z = Z,	         # z-score
        ztable = ztable, # z-score table
        cl = cl)         # classes
    }
    
    else{ # !is.null(a1 & !is.null(a2)
      ## Calculate statistic based on the given parameter values
      ## and the corresponding FDR
      fit <- testStatistic(paired, lapply(split(1:length(cl), cl), function(x) data[,x]))
      d <- fit$d / (a1 + a2 * fit$s)
      message("Calculating p-values")
      p <- calculateP(d, pD/(a1 + a2 * pS))
      message("Calculating FDR")
      FDR <- calculateFDR(d, pD/(a1 + a2 * pS), progress)
      
      ROTS.output <- list(
        data = data,   # Original data
        B = B,		     # Number of resamplings
        d = d,		     # ROTS-statistic
        logfc = logfc, # Fold change
        pvalue = p,	   # p-value
        FDR = FDR,	   # FDR
        a1 = a1,	     # a1
        a2 = a2,	     # a2
        k = NULL,	     # top list size
        R = NULL,	     # reproducibility value
        Z = NULL,	     # z-score
        cl = cl)       # classes
    }
    
    ## Define the class
    class(ROTS.output) <- "ROTS"
    return(ROTS.output)
  }
