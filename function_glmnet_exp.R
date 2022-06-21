# ------------------------------------------------------------------------------------------------
# Contains wrapper function, utility function for running experiment using glmnet library
# ------------------------------------------------------------------------------------------------
require(glmnet)
require(survival)
# ------------------------------------------------------------------------------------------------
# Impute missing values in baseline covariates (for numeric only) with mean values
# ------------------------------------------------------------------------------------------------
Imputate.x.mean <- function(x.covs, verbose = FALSE) {
  # Function to fill missing values with mean values
  
  if (verbose) {
    print("----------------------------------------------------------------------")
    print("Check missing before imputation:")
    x.covs %>%
      is.na() %>%
      colSums() %>%
      print()
    print("----------------------------------------------------------------------")
  }
  
  for (i in colnames(x.covs)) {
    #print("----------------------------------------------------------------------")
    if (verbose) {
      print(paste("Check missing values in", i))      
    }
    if (is.numeric(x.covs[, i][[1]])) {
      # Compute median values
      impute <- mean(x.covs[, i][[1]], na.rm = TRUE)
      isNA <- is.na(x.covs[, i][[1]])
      
      if (verbose) {
        print(paste("Median =", median(x.covs[, i][[1]], na.rm = TRUE)))
        print(paste("Mean =", mean(x.covs[, i][[1]], na.rm = TRUE)))
      }
      
      if (any(isNA)) {
        x.covs[, i][[1]][isNA] <- impute
        
        if (verbose) {
          print(paste("[Report] found", sum(isNA), "NA values in", i))
          print(paste("Mean after imputation =", mean(x.covs[, i][[1]])))        
        }
        
      } else {
        if (verbose) {
          print("No imputation required.")  
        }
      }
    } else {
      isNA <- is.na(x.covs[, i][[1]])
      
      if (verbose) {
        print(paste(i, "is not numeric."))
      }

      if (any(isNA)) {
        stop("Not implemented!") # Action for this case to be supplied.
      } else {
        if (verbose) {
          print("No imputation required.")
        }
      }
    }
    #print("----------------------------------------------------------------------")
  }
  
  if (verbose) {
    print("----------------------------------------------------------------------")
    print("Check missing after imputation:")
    x.covs %>%
      is.na() %>%
      colSums() %>%
      print()
    print("----------------------------------------------------------------------")
  }

  return (x.covs)
}
# ------------------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------------------
# Wrapper to fit pCox model with input formatting
# ------------------------------------------------------------------------------------------------
Train_glmnet <- function(
    training.surv,
    training.long,
    scenario,
    baseline.covs,
    baseline.covs.additional
  ) {
  # -----------------------------------------------------------------------------------
  # pCox - prepare data
  # -----------------------------------------------------------------------------------
  
  if (scenario == "scenario1" | scenario == "scenario2") {
    covs.pcox <- c(baseline.covs, baseline.covs.additional)
  } else if (scenario == "scenario0") {
    covs.pcox <- c(baseline.covs)
  } else {
    stop("undefined scenario!")
  }
  
  # Select covariates
  training.x.covs <- training.surv %>%
    select(all_of(covs.pcox))
  # Imputation using train set values
  training.x.covs.imputed <- Imputate.x.mean(training.x.covs)
  
  # Model matrix, of dimension n obs x n vars; each row is an observation vector
  training.x.mat <- model.matrix(~ ., data = training.x.covs.imputed)
  
  # Response
  training.y <- survival::Surv(
    time = training.surv$time,
    event = training.surv$event,
    type = "right"
  )
  
  # -----------------------------------------------------------------------------------  
  # pCox - Fit model
  # -----------------------------------------------------------------------------------
  
  t.start <- Sys.time()
  
  # Model selection i.e. optimize penalty parameter using cross-validation
  cvfit <- cv.glmnet(
    x = training.x.mat,
    y = training.y,
    family = "cox",
    alpha = 1, # Hyperparam
    nfolds = 10, # Hyperparam
    type.measure = "C" # Harrel's concordance, only available for cox models
  )
  
  t.end <- Sys.time()
  t.step1 <- as.numeric(difftime(t.end, t.start, units = "mins"))
  # -----------------------------------------------------------------------------------
  
  # -----------------------------------------------------------------------------------
  # Store results
  # -----------------------------------------------------------------------------------
  res <- list(
    cvfit = cvfit,
    covs.pcox = covs.pcox,
    runtimes = list(
      step1 = t.step1,
      total = t.step1)
  )
}