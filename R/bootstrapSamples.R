`bootstrapSamples` <-
function(data, B, labels, paired){
  
  samples <- matrix(nrow=B, ncol=length(labels))
  
  for(i in 1:B){
    for (label in unique(labels)) {
      pos <- which(labels==label)
      if (length(pos)>1) {
        samples[i,pos] <- sample(pos, length(pos), replace=TRUE)
      } else {
        samples[i,pos] <- pos
      }
    }
  }
  
  if (paired) {
    for(i in 1:B){
      for (label in unique(labels)[-1]) {
        pos <- which(labels==label)
        samples[i,pos] <- samples[i,which(labels==1)]+pos[1]-1
      }
    }
  }
  
  return(samples)
}

`bootstrapSamples.surv` <-
function(data, B){
  
  samples <- matrix(nrow=B, ncol=ncol(data))
  for(i in seq_len(B)){
    samples[i,] <- sort(sample(seq_len(ncol(data)), replace=TRUE), decreasing=FALSE)
  }
  
  return(samples)
}