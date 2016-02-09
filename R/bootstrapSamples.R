`bootstrapSamples` <-
function(data, B, labels){
  l1 <- which(labels==1)
  l2 <- which(labels==2)
  
  samples <- matrix(nrow=B, ncol=length(l1)+length(l2))
  for(i in seq_len(B)){
    j <- 1
    while(j==1){
      a1 <- sample(l1, length(l1), replace=TRUE)
      j <- length(unique(a1))
    }
    
    j <- 1
    while(j==1){
      a2 <- sample(l2, length(l2), replace=TRUE)
      j <- length(unique(a2))
    }
    samples[i,] <- c(a1,a2)
  }
  
  return(samples)
}

