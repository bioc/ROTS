# Function to run the lm model
`runlm` <- function(i, formula, data, metadata) {
  tryCatch({
    fit <- suppressMessages(suppressWarnings(lm(as.formula(paste("datavalue ~",formula)), data=cbind(metadata,datavalue=c(t(data[i,]))))))
    coef <- coefficients(summary(fit))
    co <- coef[-1,1]; names(co) <- paste("coef",names(co),sep=".")
    sd <- coef[-1,2]; names(sd) <- paste("sd",names(sd),sep=".")
    return(c(co,sd))
  }, error = function(e) {
    return(NA)
  })
}

# lmROTS
`lmROTS` <- function(formula, data, metadata, B=100, K=NULL, seed=NULL, BPPARAM=bpparam()) {
  if (is(data, "ExpressionSet"))
    data <- Biobase::exprs(data)
  	
  # Set bootstraps and permutations
  if(!is.null(seed)) {
    set.seed(seed, kind="default")
  }
  boot <- sapply(1:B, function(x) sample(1:nrow(metadata),nrow(metadata),replace=TRUE))
  perm <- sapply(1:B, function(x) sample(1:nrow(metadata),nrow(metadata),replace=FALSE))
  
  # Original
  message("Running initial model")
  lm.original <- bplapply(1:nrow(data), function(i) runlm(i, formula, data, metadata), BPPARAM=BPPARAM)
  names <- names(lm.original[[which.max(sapply(lm.original, length))]])
  for (i in 1:length(lm.original)) {
    lm.original[[i]] <- lm.original[[i]][names]
  }
  lm.original <- do.call("rbind", lm.original)
  colnames(lm.original) <- names
  rownames(lm.original) <- rownames(data)
  
  # Run over bootstraps
  message("Running bootstraps")
  lm.boot <- bplapply(1:B, function(x) {
    out <- lapply(1:nrow(data), function(i) runlm(i, formula, data[,boot[,x]], metadata[boot[,x],]))
    out[is.na(out)] <- list(rep(NA,ncol(lm.original)))
    out <- do.call("rbind",out)
    out <- out[,match(colnames(lm.original),colnames(out))]
    colnames(out) <- colnames(lm.original)
    return(out)
  }, BPPARAM=BPPARAM)
  
  # Run over permutations
  message("Running permutations")
  lm.null <- bplapply(1:B, function(x) {
    out <- lapply(1:nrow(data), function(i) runlm(i, formula, data, metadata[perm[,x],]))
    out[is.na(out)] <- list(rep(NA,ncol(lm.original)))
    out <- do.call("rbind",out)
    out <- out[,match(colnames(lm.original),colnames(out))]
    colnames(out) <- colnames(lm.original)
    return(out)
  }, BPPARAM=BPPARAM)
  
  # Optimize parameters
  ROTS.output <- optimizeModel(data=data, model.original=lm.original, model.boot=lm.boot, model.null=lm.null, B=B, K=K, seed=seed, BPPARAM=BPPARAM)
  class(ROTS.output) <- "regROTS"
  return(ROTS.output)
}
