# Function to run the lm model
runlm <- function(i, formula, data, metadata) {
  tryCatch({
    fit <- suppressMessages(suppressWarnings(lm(as.formula(paste("datavalue ~",formula)), data=cbind(metadata,datavalue=c(t(data[i,]))))))
    coef <- coefficients(summary(fit))
    return(c(coef[-1,1],coef[-1,2]))
  }, error = function(e) {
    return(NULL)
  })
}

# lmROTS
lmROTS <- function(formula, data, metadata, B=100, K=NULL, seed=NULL, BPPARAM=bpparam()) {
  
  # Set bootstraps and permutations
  if(!is.null(seed)) {
    set.seed(seed, kind="default")
  }
  boot <- sapply(1:B, function(x) sample(1:nrow(data.meta),nrow(data.meta),replace=TRUE))
  perm <- sapply(1:B, function(x) sample(1:nrow(data.meta),nrow(data.meta),replace=FALSE))
  
  # Original
  message("Running initial model")
  lm.original <- do.call("rbind", bplapply(1:nrow(data), function(i) runlm(i, formula, data, metadata), BPPARAM=BPPARAM))
  
  # Run over bootstraps
  message("Running bootstraps")
  lm.boot <- bplapply(1:B, function(x) {
    t(sapply(1:nrow(data), function(i) runlm(i, formula, data[,boot[,x]], metadata[boot[,x],])))
  }, BPPARAM=BPPARAM)
  
  # Run over permutations
  message("Running permutations")
  lm.null <- bplapply(1:B, function(x) {
    t(sapply(1:nrow(data), function(i) runlm(i, formula, data, metadata[perm[,x],])))
  }, BPPARAM=BPPARAM)
  
  # Optimize parameters
  ROTS.output <- optimizeModel(model.original=lm.original, model.boot=lm.boot, model.null=lm.null, B=B, K=K, seed=seed, BPPARAM=BPPARAM)
  class(ROTS.output) <- "regROTS"
  return(ROTS.output)
}
