# ------------------------------------------------------------------------------------------------
# Contains wrapper function, utility function for running experiment using pencal library
# ------------------------------------------------------------------------------------------------
#install.packages("./package/pencal_1.2.2.tar.gz", repos = NULL, type="source")
require(tidyverse)
require(pencal)
print(paste("pencal version =", packageVersion("pencal")))
# ------------------------------------------------------------------------------------------------
# Wrapper function - Fit PRC model steps 1 to 3
# ------------------------------------------------------------------------------------------------
run_prc_steps <- function(
    long.data, 
    surv.data,
    baseline.covs,
    y.names,
    max.ymissing = 0.2,
    n.boots = 0, 
    n.cores = 1,
    verbose = FALSE,
    penalty = "ridge",
    standardize = TRUE,
    pfac.base.covs = 0 # vector of binary value(s) to indicate regularization for baseline covariates
  ) {
  
  # Fit PRC model steps 1 to 3
  # Input:
  #   long.data: long data
  #   surv.data: survival data
  #   y.names: vector of candidate longitudinal variable names
  # Output:

  # -----------------------------------------------------------------------------------
  # Step 1 of PRC-LMM: estimate the LMMs
  t.start <- Sys.time()
  
  step1 <- fit_lmms(
    y.names = y.names,
    fixefs = ~ age.fup, # Fixed effects
    ranefs = ~ age.fup | id, # Random effects
    long.data = long.data,
    surv.data = surv.data,
    t.from.base = Years.bl, # Name of the variable containing time from baseline in long.data
    max.ymissing = max.ymissing,
    n.boots = n.boots, 
    n.cores = n.cores)
  
  t.end <- Sys.time()
  t.step1 <- as.numeric(difftime(t.end, t.start, units = "mins"))
  
  if(verbose){
    cat("\nStep 1 completed in", t.step1, "min\n")
    cat("--------------------------------------------------------------------------------\n")
  }
  # -----------------------------------------------------------------------------------
  # Step 2 of PRC-LMM: compute the summaries
  t.start <- Sys.time()
  
  step2 <- summarize_lmms(object = step1, n.cores = n.cores)
  
  t.end <- Sys.time()
  t.step2 <- as.numeric(difftime(t.end, t.start, units = "mins"))
  
  if(verbose){
    cat("\nStep 2 completed in", t.step2, "min\n")
    cat("--------------------------------------------------------------------------------\n")
  }
  # -----------------------------------------------------------------------------------
  # Step 3 of PRC-LMM: fit the penalized Cox models
  t.start <- Sys.time()
  
  step3 <- fit_prclmm(
    object = step2, 
    surv.data = surv.data,
    baseline.covs = as.formula(paste("~", paste(baseline.covs, collapse = "+"))),
    penalty = penalty,
    standardize = standardize,
    pfac.base.covs = pfac.base.covs, # Change it to a vector to penalize baseline covariates
    n.cores = n.cores
    )
  
  t.end <- Sys.time()
  t.step3 <- as.numeric(difftime(t.end, t.start, units = "mins"))
  
  if (verbose){
    cat("\nStep 3 completed in", t.step3, "min\n")
    cat("--------------------------------------------------------------------------------\n")
  }
  # -----------------------------------------------------------------------------------
  return (list(
    step1 = step1,
    step2 = step2,
    step3 = step3,
    runtimes = list(
      step1 = t.step1,
      step2 = t.step2,
      step3 = t.step3,
      total = t.step1 + t.step2 + t.step3
    )
    ))
}
# ------------------------------------------------------------------------------------------------

