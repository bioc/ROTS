# Function to run the coxph model
`runcoxph` <- function(i, formula, data, metadata) {
  tryCatch({
    fit <- suppressMessages(suppressWarnings(survival::coxph(formula, data=cbind(metadata,exprs=c(t(data[i,]))))))
    coef <- coefficients(summary(fit))
    co <- coef[,1]; names(co) <- paste("coef",rownames(coef),sep=".")
    se <- coef[,3]; names(se) <- paste("se",rownames(coef),sep=".")
    return(c(co,se))
  }, error = function(e) {
    return(NA)
  })
}

# coxphROTS
`coxphROTS` <- function(formula, data, metadata, B=100, K=NULL, seed=NULL, a1=NULL, a2=NULL, BPPARAM=bpparam()) {
  if (is(data, "ExpressionSet"))
    data <- Biobase::exprs(data)
  
  # Check variables in formula
  formula <- as.formula(formula)
  if(!("exprs" %in% all.vars(formula))) {
    stop("Term 'exprs' missing from the formula.")
  }
  
  # Set bootstraps and permutations
  if(!is.null(seed)) {
    set.seed(seed, kind="default")
  }
  boot <- sapply(1:B, function(x) sample(1:nrow(metadata),nrow(metadata),replace=TRUE))
  perm <- sapply(1:B, function(x) sample(1:nrow(metadata),nrow(metadata),replace=FALSE))
  
  # Original
  message("Running initial model")
  coxph.original <- bplapply(1:nrow(data), function(i) runcoxph(i, formula, data, metadata), BPPARAM=BPPARAM)
  names <- names(coxph.original[[which.max(sapply(coxph.original, length))]])
  for (i in 1:length(coxph.original)) {
    coxph.original[[i]] <- coxph.original[[i]][names]
  }
  coxph.original <- do.call("rbind", coxph.original)
  colnames(coxph.original) <- names
  rownames(coxph.original) <- rownames(data)

  # Run over bootstraps
  message("Running bootstraps")
  coxph.boot <- bplapply(1:B, function(x) {
    out <- lapply(1:nrow(data), function(i) runcoxph(i, formula, data[,boot[,x]], metadata[boot[,x],,drop=FALSE]))
    out[is.na(out)] <- list(rep(NA,ncol(coxph.original)))
    out <- do.call("rbind",out)
    out <- out[,match(colnames(coxph.original),colnames(out))]
    colnames(out) <- colnames(coxph.original)
    return(out)
  }, BPPARAM=BPPARAM)
  
  # Run over permutations
  message("Running permutations")
  coxph.null <- bplapply(1:B, function(x) {
    out <- lapply(1:nrow(data), function(i) runcoxph(i, formula, data, metadata[perm[,x],,drop=FALSE]))
    out[is.na(out)] <- list(rep(NA,ncol(coxph.original)))
    out <- do.call("rbind",out)
    out <- out[,match(colnames(coxph.original),colnames(out))]
    colnames(out) <- colnames(coxph.original)
    return(out)
  }, BPPARAM=BPPARAM)
  
  # Optimize parameters
  ROTS.output <- optimizeModel(data=data, model.original=coxph.original, model.boot=coxph.boot, model.null=coxph.null, B=B, K=K, seed=seed, a1=a1, a2=a2, BPPARAM=BPPARAM)
  class(ROTS.output) <- "regROTS"
  return(ROTS.output)
}
