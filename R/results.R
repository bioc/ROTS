`results` <-
function(object, order=FALSE, pvalue=NULL, FDR=NULL, coef=NULL, logfc=NULL) {
  if(class(object)!="ROTS" & class(object)!="regROTS") {
    stop("Function works only with ROTS or regROTS objects.")
  }
  
  # Gather results from regROTS object
  if (class(object)=="regROTS") {
    out <- lapply(object, function(x) data.frame(feature=names(x$d), d=x$d, coef=x$coef, pvalue=x$pvalue, FDR=x$FDR))
    if(!is.null(coef)) {
      if(length(coef)==1) {
        coef <- rep(coef,length(out))
      }
      if(length(coef)!=length(out)) {
        stop(paste("Coefficient filter needs to be a single value or a vector with a length of ", length(out), sep=""))
      }
    }
    for(i in seq_along(out)) {
      if(order) {
        out[[i]] <- out[[i]][order(abs(out[[i]]$d), abs(out[[i]]$coef), decreasing=TRUE),]
      }
      if(!is.null(pvalue)) {
        out[[i]] <- out[[i]][which(out[[i]]$pvalue <= pvalue),]
      }
      if(!is.null(FDR)) {
        out[[i]] <- out[[i]][which(out[[i]]$FDR <= FDR),]
      }
      if(!is.null(coef)) {
        out[[i]] <- out[[i]][which(abs(out[[i]]$coef) >= coef[i]),]
      }
    }
    return(out)
  }
  
  # Gather results from ROTS object
  if (class(object)=="ROTS") {
    out <- data.frame(feature=names(object$d), d=object$d, logfc=object$logfc, pvalue=object$pvalue, FDR=object$FDR)
    if(order) {
      out<- out[order(abs(out$d), abs(out$logfc), decreasing=TRUE),]
    }
    if(!is.null(pvalue)) {
      out<- out[which(out$pvalue <= pvalue),]
    }
    if(!is.null(FDR)) {
      out <- out[which(out$FDR <= FDR),]
    }
    if(!is.null(logfc)) {
      out<- out[which(abs(out$logfc) >= logfc),]
    }
    return(out)
  }
  
}
