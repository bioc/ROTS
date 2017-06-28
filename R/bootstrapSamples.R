`bootstrapSamples` <-
function(data, B, labels, paired){
  
  samples <- matrix(nrow=B, ncol=length(labels))
  
  for(i in 1:B){
    for (label in unique(labels)) {
      pos <- which(labels==label)
      samples[i,pos] <- sample(pos, length(pos), replace=TRUE)
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
