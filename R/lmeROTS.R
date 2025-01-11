# Function to run the lme model
`runlme` <- function(i, formula, data, metadata) {
  tryCatch({
    fit <- suppressMessages(suppressWarnings(lme4::lmer(as.formula(paste("datavalue ~",formula)), data=cbind(metadata,datavalue=c(t(data[i,]))))))
    coef <- coefficients(summary(fit))
    co <- coef[-1,1]; names(co) <- paste("coef",names(co),sep=".")
    sd <- coef[-1,2]; names(sd) <- paste("sd",names(sd),sep=".")
    return(c(co,sd))
  }, error = function(e) {
    return(NA)
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
  lme.original <- bplapply(1:nrow(data), function(i) runlme(i, formula, data, metadata), BPPARAM=BPPARAM)
  names <- names(lme.original[[which.max(sapply(lme.original, length))]])
  for (i in 1:length(lme.original)) {
    lme.original[[i]] <- lme.original[[i]][names]
  }
  lme.original <- do.call("rbind", lme.original)
  colnames(lme.original) <- names
  rownames(lme.original) <- rownames(data)

  # Run over bootstraps
  message("Running bootstraps")
  lme.boot <- bplapply(1:B, function(x) {
    out <- lapply(1:nrow(data), function(i) runlme(i, formula, data[,boot[,x]], metadata[boot[,x],]))
    out[is.na(out)] <- list(rep(NA,ncol(lme.original)))
    out <- do.call("rbind",out)
    out <- out[,match(colnames(lme.original),colnames(out))]
    colnames(out) <- colnames(lme.original)
    return(out)
  }, BPPARAM=BPPARAM)
  
  # Run over permutations
  message("Running permutations")
  lme.null <- bplapply(1:B, function(x) {
    out <- lapply(1:nrow(data), function(i) runlme(i, formula, data, metadata[perm[,x],]))
    out[is.na(out)] <- list(rep(NA,ncol(lme.original)))
    out <- do.call("rbind",out)
    out <- out[,match(colnames(lme.original),colnames(out))]
    colnames(out) <- colnames(lme.original)
    return(out)
  }, BPPARAM=BPPARAM)
  
  # Optimize parameters
  ROTS.output <- optimizeModel(data=data, model.original=lme.original, model.boot=lme.boot, model.null=lme.null, B=B, K=K, seed=seed, BPPARAM=BPPARAM)
  class(ROTS.output) <- "regROTS"
  return(ROTS.output)
}
