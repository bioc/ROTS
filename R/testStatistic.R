`testStatistic` <-
 function(paired, samples){
   
   # Two groups
   if (length(samples)==2) {
     X <- samples[[1]]
     Y <- samples[[2]]
     
     ## Calculates the test statistic for each row.
     ## X and Y are the data matrices of the two groups.
     ## Each row of these two matrices must contain at least TWO not NA values.
     ## Thus the "variance" always exists.  
     
     ## Row means
     mX <- rowMeans(X, na.rm=TRUE)
     mY <- rowMeans(Y, na.rm=TRUE)
     
     ## Pooled standard deviations for each row
     sX <- rowSums((X - mX)^2, na.rm=TRUE)
     sY <- rowSums((Y - mY)^2, na.rm=TRUE)
     
     if(!paired) {
       ## Number of not NA values in each row
       nX <- rowSums(!is.na(X))
       nY <- rowSums(!is.na(Y))
       
       ## d == difference between the group means for each row (==gene)
       ## s == pooled standard deviation for each row (==gene)        
       d <- mY - mX
       s <- sqrt(((sX + sY) / (nX + nY - 2)) * (1 / nX + 1 / nY))
       
       ## Cases with less than two non-missing values.
       ## Set d = 0, s = 1
       ind <- which( nY < 2 | nX < 2 )
       d[ind] <- 0
       s[ind] <- 1
     }
     
     if(paired) {
       ## Add for paired
       sXY <- rowSums((X - mX)*(Y - mY), na.rm=TRUE)
       
       ## Number of not NA values in each row
       n <- rowSums(!is.na(X*Y))
       
       ## d == difference between the group means for each row (==gene)
       ## s == pooled standard deviation for each row (==gene)        
       d <- mY - mX
       s <- sqrt(((sX + sY) / (n + n - 2)) * (2 / n) - 2/(n*n-n)*sXY)
       
       ## Cases with less than two non-missing values.
       ## Set d = 0, s = 1
       ind <- which( n < 2 )
       d[ind] <- 0
       s[ind] <- 1
     }
     
     return(list(d=d, s=s))
   }
   
   # Multiple groups
   if (length(samples)>2) {
     
     samples.all <- do.call("cbind",samples)
     
     if(!paired) {
       f <- sum(sapply(samples, ncol)) / prod(sapply(samples, ncol))
       r <- vector(mode="numeric", length=nrow(samples.all))
       for(k in 1:length(samples)) {
         r <- r + (rowMeans(samples[[k]], na.rm=TRUE)-rowMeans(samples.all, na.rm=TRUE))^2
       }
       d <- (f*r)^0.5
       
       f <- 1/sum(sapply(samples, ncol)-1) * sum(1/sapply(samples, ncol))
       s <- vector(mode="numeric", length=nrow(samples.all))
       for(k in 1:length(samples)) {
         s <- s + colSums(apply(samples[[k]], 1, function(x) (x-mean(x,na.rm=TRUE))^2), na.rm=TRUE)
       }
       s <- (f*s)^0.5
       
     }
     
     if(paired) {
       stop("Multiple paired groups not supported!")
     }
     
     return(list(d=d, s=s))
   }
   
}

`testStatistic.surv` <- function(samples, time, event){
  samples.all <- do.call("cbind",samples)
  t <- unique(time[event==1])
  
  r <- vector(mode="numeric", length=nrow(samples.all))
  for(k in t) {
    i <- which(time>=k)
    z <- which(time==k)
    d <- z[which(event[which(time==k)]==1)]
    if (length(i)>1) {
      r <- r + (rowSums(as.data.frame(samples.all[,d]), na.rm=TRUE)-length(d)*rowMeans(samples.all[,i], na.rm=TRUE))
    }
  }
  
  s <- vector(mode="numeric", length=nrow(samples.all))
  for(k in t) {
    i <- which(time>=k)
    z <- which(time==k)
    d <- z[which(event[which(time==k)]==1)]
    if (length(i)>1) {
      s <- s + ((length(d)/length(i)) * rowSums((samples.all[,i]-rowMeans(samples.all[,i], na.rm=TRUE))^2, na.rm=TRUE))
    }
  }
  s <- s^0.5
  
  return(list(d=r, s=s))
}



