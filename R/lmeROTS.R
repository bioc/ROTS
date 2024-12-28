# Function to run the lme model
`runlme` <- function(i, formula, data, metadata) {
  tryCatch({
    fit <- suppressMessages(suppressWarnings(lme4::lmer(as.formula(paste("datavalue ~",formula)), data=cbind(metadata,datavalue=c(t(data[i,]))))))
    coef <- coefficients(summary(fit))
    return(c(coef[-1,1],coef[-1,2]))
  }, error = function(e) {
    return(NULL)
  })
}

# lmeROTS
`lmeROTS` <- function(formula, data, metadata, B=100, K=NULL, seed=NULL, BPPARAM=bpparam()) {
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
  lme.original <- do.call("rbind", bplapply(1:nrow(data), function(i) runlme(i, formula, data, metadata), BPPARAM=BPPARAM))
  
  # Run over bootstraps
  message("Running bootstraps")
  lme.boot <- bplapply(1:B, function(x) {
    t(sapply(1:nrow(data), function(i) runlme(i, formula, data[,boot[,x]], metadata[boot[,x],])))
  }, BPPARAM=BPPARAM)
  
  # Run over permutations
  message("Running permutations")
  lme.null <- bplapply(1:B, function(x) {
    t(sapply(1:nrow(data), function(i) runlme(i, formula, data, metadata[perm[,x],])))
  }, BPPARAM=BPPARAM)
  
  # Optimize parameters
  ROTS.output <- optimizeModel(data=data, model.original=lme.original, model.boot=lme.boot, model.null=lme.null, B=B, K=K, seed=seed, BPPARAM=BPPARAM)
  class(ROTS.output) <- "regROTS"
  return(ROTS.output)
}
