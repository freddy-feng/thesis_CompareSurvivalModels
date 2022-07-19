# ------------------------------------------------------------------------------------------------
# Contains wrapper function, utility function for running experiment using MFPCCox method
# ------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------
# [superseded] Manual cleaning to avoid double entry problem in data
# ------------------------------------------------------------------------------------------------
Remove_invalid_long_data <- function(long) {
  checkcond <- any(training.long[training.long$id %in% 6014, "VISCODE"] == "m0")
  # Equivalent to
  #long <- long %>%
  #  filter(!(id == 6014 & VISCODE == "m0"))
  
  if (checkcond) {
    long <- long %>%
      filter(!(id == 6014 & VISCODE == "m0"))
    print("[Report] Manually removed observation in subject 6014, VISCODE m0.")
  }
  
  return (long)
}
# ------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------
# Transfer longitudinal covariates from long data to 3-dim array
# ------------------------------------------------------------------------------------------------
Convert_long_to_mvarray <- function(long, y.t, obstime, subject_id, n_subject, y.names, n.y) {
  
  # Initialize
  multivar <- array(
    NA, # NA values
    c(n_subject, # dim1 is number of subjects for train + test i.e. all available data
      length(obstime), # dim2 is observation times
      n.y) # dim3 is number of longitudinal covariates
  ) 
  
  # For each subject i
  for (s in 1:n_subject) { 
#    print(s) # Uncomment to debug
    visits <- which(obstime %in% (long[long$id == subject_id[s], y.t])) # get obstime as visit index for each subject i, use as index below
    for(p in 1:n.y){ # for each predictor
      multivar[s, visits, p] <- long[, y.names[p]][long$id == subject_id[s]] # filter obs in covariate p for subject i
    }
  }
  
  return (multivar)
}
# ------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------
# Transfer longitudinal covariates from long data to 3-dim array
# ------------------------------------------------------------------------------------------------
Train_MFPCCox <- function(
    training.surv,
    multivar,
    subject.id,
    y.names,
    baseline.covs,
    argvals,
    pve,
    nbasis
    ) {
  # -----------------------------------------------------------------------------------  
  # MFPCCpx - Get multivar array
  # -----------------------------------------------------------------------------------
  ids.train <- subject.id %in% training.surv$id # Get row indices for train set
  multivar.train <- multivar[ids.train, , ] # Subset for train set

  n.y <- length(y.names)
  # -----------------------------------------------------------------------------------
  # MFPCCox Step 1 - Apply univariate FPCA using PACE algorithm
  # -----------------------------------------------------------------------------------
  print("Step 1 - univariate FPCA")
  t.start <- Sys.time()
  
  # Initialize variables
  Xi.train <- NULL # Xi = FPC scores
  L <- NULL
  phi.train <- NULL # phi = eigenfunctions
  meanFun.train <- NULL # meanFun = mean function
  npc.train <- NULL
  
  for(p in 1:n.y) { # for each long covariate
    print(paste("Apply uFPCA to", y.names[p])) # Uncomment to debug
    
    # wrapper for PACE
    tmp.ufpca <- uPACE(
      testData = multivar.train[, , p], 
      domain = argvals, # scaled obstime in [0,1]
      pve = pve,
      nbasis = nbasis) 
    # nbasis: representing the number of B-spline basis functions used for estimation of the mean function and bivariate smoothing of the covariance surface. Defaults to 10
    
    # Detail about output from uPACE refer to doc of MFPCA library
    Xi.train <- cbind(Xi.train, tmp.ufpca$scores) # estimated FPC scores --> column bind
    L <- c(L, dim(tmp.ufpca$scores)[2]) # vector of Lq
    # ?"@" Extract or replace the contents of a slot in a object with a formal (S4) class structure.
    phi.train[[p]] <- t(tmp.ufpca$functions@X) # estimated functional principal components (eigenfunctions)
    meanFun.train[[p]] <- tmp.ufpca$mu@X # mean functions
    npc.train[[p]] <- tmp.ufpca$npc
  }
  
  t.end <- Sys.time()
  t.step1 <- as.numeric(difftime(t.end, t.start, units = "mins"))
  
  # -----------------------------------------------------------------------------------
  # MFPCCpx Step 2 - Apply multivariate FPCA on uFPCA outputs
  # -----------------------------------------------------------------------------------
  print("Step 2 - MFPCA")
  t.start <- Sys.time()
  
  mFPCA.train <- mFPCA(
    Xi = Xi.train, # Xi = FPC scores
    phi = phi.train, # phi = eigenfunctions
    p = n.y,
    L = L, # vector of Lq
    I = length(training.surv$id)
  )
  
  rho.train <- mFPCA.train$rho # MFPC scores
  pve <- mFPCA.train$pve # percentage of variance explained
  psi <- mFPCA.train$psi # multivariate eigenfunctions
  Cms <- mFPCA.train$Cms # eigenvectors
  
  t.end <- Sys.time()
  t.step2 <- as.numeric(difftime(t.end, t.start, units = "mins"))
  # -----------------------------------------------------------------------------------
  # MFPCCpx Step 3 - Fit Cox model
  # -----------------------------------------------------------------------------------
  print("Step 3 - Cox model")
  t.start <- Sys.time()
  
  # Extract scores
  tmp.rho <- data.frame(rho.train)
  rho.names <- paste0("rho", 1:ncol(rho.train))
  names(tmp.rho) <- rho.names
  training.surv.rho <- data.frame(training.surv, tmp.rho)
  
  lhs <- "Surv(time, event)"
  rhs <- paste(c(baseline.covs, rho.names), collapse = "+")
  formula <- as.formula(paste(lhs, "~", rhs))
  
  mfpccox <- coxph(
    formula, 
    data = training.surv.rho,
    model = TRUE,
    x = TRUE,
    y = TRUE)
  
  t.end <- Sys.time()
  t.step3 <- as.numeric(difftime(t.end, t.start, units = "mins"))
  
  return (list(
    mfpccox = mfpccox,
    mfpca.train = mFPCA.train,
    phi.train = phi.train,
    npc.train = npc.train,
    runtimes = list(
      step1 = t.step1,
      step2 = t.step2,
      step3 = t.step3,
      total = t.step1 + t.step2 + t.step3
    )))
}