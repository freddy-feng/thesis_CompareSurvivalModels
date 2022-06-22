library(tidyverse)

source("function_utility.R")
source("function_evaluation.R")

source("function_glmnet_exp.R")

select <- dplyr::select
# -----------------------------------------------------------------------------------
for (seed in 721:730) {
  # -----------------------------------------------------------------------------------
  # Select either method
#  method <- "pCox-bl"
  method <- "pCox-lm"
  # -----------------------------------------------------------------------------------
  # Set experiment parameters here!
  # -----------------------------------------------------------------------------------
  seed <- seed # Vary this for repeated CV
  print(paste("This CV will perform train-test split using seed", seed))
  
  n_fold <- 10 # Cross validation
  
  set_scenario <- "scenario2" # Determine how many longitudinal covariates to use
  
  is_transformed <- "transformed"
  is_scaled <- "scaled"
  
  # Baseline covariates (not time-varying in adnimerge)
  baseline.covs <- c("AGE", "PTGENDER", "PTEDUCAT", "status.bl", "APOE4") # b5 
#  baseline.covs <- c("AGE", "PTGENDER", "PTEDUCAT", "APOE4") # b4
  # -----------------------------------------------------------------------------------
  # -----------------------------------------------------------------------------------
  
  # -----------------------------------------------------------------------------------
  # Load data
  # -----------------------------------------------------------------------------------
  # Load cleaned data
  path.data <- "./data_cleaned/adni_cleaned.RData"
  load(path.data)
  
  # [future] may also reduce the number of columns here to reduce size
  # Note: data.surv and data.long are arranged by id and {id, age.fup} to ensure properly use pencal
  data.surv <- df.surv_preds
  
  if (is_transformed == "transformed") {
    data.long <- df.long_censored_transformed
  } else {
    data.long <- df.long_censored
  }
  # -----------------------------------------------------------------------------------
  # Set up for landmarking and evaluation
  # -----------------------------------------------------------------------------------
  # landmark time
  T.start <- 3
  landmark <- paste0("lm", T.start) # File name description
  
  T.max <- floor(max(data.surv$time)) # Based on last available observation in train set
  # [Warning] foresee a potential bug may happen by chance if the longest observation is in test set, but not train
  
  # Predict 15 years onward from landmark time, but not more than max observed time
  deltaT <- 1:T.max # A vector of prediction times, starting from baseline onward
  deltaT <- deltaT[deltaT > T.start]
  
  # [Future] May need to change deltaT later because the notation is inconsistent with the symbol in paper
  # the true delta T should be prediction time - landmark time! Avoid confusion!
  
  print(paste("[Report] Evaluation on times since baseline (time=0):", paste(deltaT, collapse = " ")))
  # -----------------------------------------------------------------------------------
  
  # Select subjects at risk since landmark time
  data.surv <- data.surv %>%
    filter(time > T.start)
  data.long <- data.long %>%
    filter(time > T.start)
  
  # Remove the repeated observations after the landmark time
  data.long <- data.long %>%
    filter(Years.bl <= T.start)
  
  print(paste("Number of subject at risk after landmark time =", nrow(data.surv)))
  print(paste("Number of visits before landmark time (upper bound of measurements) =", nrow(data.long)))
  
  
  
  # -----------------------------------------------------------------------------------
  # Initialize
  # -----------------------------------------------------------------------------------
  
  
  
  # Manual exclusion
  vars_manual_remove <- c("TAU", "PTAU", "ABETA")
  # Note: Type of TAU, PTAU and ABETA are character
  # need to handle non-numerical values first. currently excluded
  
  # Exclude irrelevant variables
  vars_irrelevant <- c(
    names(data.long)[grepl(".bl", names(data.long))], # Exclude variables with `.bl` suffix including Years.bl and Months.bl
    "id", "RID",
    "time", "event", "status", "DX", # Survival information
    "VISCODE", "EXAMDATE", "Y", "M", "Month", # Time variables
    "AGE", "age.fup", 
    "COLPROT", "ORIGPROT", "PTID", "SITE", # Visit information
    "PTGENDER", "PTEDUCAT", "PTETHCAT", "PTRACCAT", "PTMARRY", "APOE4", # Baseline variables
    "FSVERSION", "IMAGEUID", "FLDSTRENG" # Metadata for image
  )
  
  vars_ignore <- c(vars_manual_remove, vars_irrelevant) # Variables that will not be considered as long covariates
  
  # -----------------------------------------------------------------------------------
  # Update values of time-varying covariates in surv data when landmark time > 0
  # for methods pCox-bl and pCox-lm
  
  # The original values observed at baseline i.e. VISCODE=="bl"
  # are replaced by last observed value on or before landmark
  
  # Step 1: set the covariates to update, should cover the candidate long covariates
  vars_long <- names(data.surv)[!(names(data.surv) %in% vars_ignore)]
  
  # Step 2: update values in surv data
  # For each subject, the latest observed value of time-varying covariate is used
  # The value can be transformed or not, depending on data.long chosen
  use_baseline <- NULL
  if (method == "pCox-bl" | method == "pCox-lm") {
    
    if (method == "pCox-bl") {
      use_baseline = TRUE
    } else {
      use_baseline = FALSE
    }
    
    data.surv <- Update_surv_at_landmark(
      surv = data.surv,
      long = data.long,
      y.names = vars_long,
      use_baseline = use_baseline)
    
    print("[Remind] pCox method is used, the additional covariates will be updated")
  }
  
  # -----------------------------------------------------------------------------------
  # Obtain a list of folds to initialize cross validation
  # Involves:
  # - Select candidate longitudinal covariates based on missingness
  # - Create folds based on stratified train-test split for n-fold CV
  
  # Get scaling table
  folds <- Initialize_exp(
    data.surv = data.surv,
    data.long = data.long,
    baseline.covs = baseline.covs,
    vars_not_long = vars_ignore,
    set_scenario = set_scenario,
    n_fold = n_fold,
    seed = seed
  )
  
  # Note: depending on `is_transformed`, either original or transformed version of data.long will be used
  # the different scaling parameters are different between these cases
  
  # -----------------------------------------------------------------------------------
  
  # -----------------------------------------------------------------------------------
  # Save folds for training and future checking
  subfolder <- "./output/temp/"
  filename <- paste0("output_folds_template_", set_scenario, "_seed", seed, "_", landmark, "_", is_transformed, ".RData")
  path.template <- paste0(subfolder, filename)
  
  save(folds, file = path.template)
  
  print(paste("template of folds saved to path:", path.template))
  
  # -----------------------------------------------------------------------------------
  # You may want to double check data.long and data.surv before proceeding to training.
  # Note that the values in data.surv may be changed in landmarking step.
  # Scaling will be carried out in training step.
  
  #Check_folds(data.surv, folds) # Uncomment to check the stratification / class balance after split
  
  
  # -----------------------------------------------------------------------------------
  # -----------------------------------------------------------------------------------
  # Set model param
  is_scaled <- is_scaled # Set to "scaled" to scale covariates for LMMs
  is_transformed <- is_transformed
  landmark <- paste0("lm", T.start)
  n_basecov <- paste0("b", length(baseline.covs))
  glmnet.alpha <- 0 # 0 for ridge, 1 for lasso
  
  penality.type <- if (glmnet.alpha == 0) {
    "ridge"
  } else if (glmnet.alpha == 1) {
    "lasso"
  } else {
    "elasticnet"
  }
  
  model.hyperparam <- list(
    method = method,
    set_scenario = set_scenario,
    landmark = landmark,
    is_scaled = is_scaled,
    is_transformed = is_transformed,
    glmnet.alpha = glmnet.alpha
  )
  
  hyperparam <- paste(c(set_scenario, n_basecov, is_transformed, is_scaled, penality.type), collapse = "_") # Use hyperparam to describe model
  model.name <- paste(c(method, landmark, hyperparam), collapse = "_")
  
  print(paste("Begin training for model:", model.name))
  # -----------------------------------------------------------------------------------
  # -----------------------------------------------------------------------------------
  rm("folds")
  rm("folds.eval")
  
  print(path.template)
  load(file = path.template) # Load folds template
  # -----------------------------------------------------------------------------------
  # Fit models for each fold
  # -----------------------------------------------------------------------------------
  for (i in 1:n_fold) {
    print("---------------------------------------------------------------------------------------------------")
    print(paste("Start training in fold", i))
    print("---------------------------------------------------------------------------------------------------")
    # -----------------------------------------------------------------------------------
    # General - Subset subjects for fold i
    # -----------------------------------------------------------------------------------
    tmp <- Get_train_test_data(
      data.surv = data.surv, 
      data.long = data.long,  
      ids.test = folds[[i]]$ids.test, 
      is_scaled = is_scaled, 
      scaling_table = folds[[i]]$scaling_table
    )
    
    training.surv <- tmp$training.surv
    training.long <- tmp$training.long
    
    # -----------------------------------------------------------------------------------
    # pCox - train model
    # -----------------------------------------------------------------------------------
    # Baseline covariates used for fitting penalized Cox model in all folds
    baseline.covs <- folds[[i]]$baseline.covs
    baseline.covs.additional <- folds[[i]]$candidate.long.covs
    
    scenario <- folds[[i]]$scenario # Deciding number candidate covariates used
    
    # Note: imputation inside
    res <- Train_glmnet(
      training.surv = training.surv,
      training.long = training.long,
      scenario = scenario,
      baseline.covs = baseline.covs,
      baseline.covs.additional = baseline.covs.additional,
      glmnet.alpha = glmnet.alpha
    )
    
    # -----------------------------------------------------------------------------------
    # Store results
    # -----------------------------------------------------------------------------------
    folds[[i]]$model <- list(
      name = model.name,
      hyperparam = model.hyperparam,
      covariate = list(
        base = res$covs.pcox),
      cvfit = res$cvfit,
      training.time = res$runtimes
    )
  }
  # -----------------------------------------------------------------------------------
  # -----------------------------------------------------------------------------------
  # Test
  folds.eval <- vector(mode = "list", length = n_fold)
  
  for (i in 1:n_fold) {
    # -----------------------------------------------------------------------------------
    # General - Subset subjects for fold i
    # -----------------------------------------------------------------------------------
    tmp <- Get_train_test_data(
      data.surv = data.surv, 
      data.long = data.long,  
      ids.test = folds[[i]]$ids.test, 
      is_scaled = is_scaled, 
      scaling_table = folds[[i]]$scaling_table
    )
    
    surv.new <- tmp$testing.surv
    long.new <- tmp$testing.long
    
    # -----------------------------------------------------------------------------------  
    # pCox - Evaluate model
    # ----------------------------------------------------------------------------------- 
    # Step 1 - prepare new X
    # Subset columns to construct model matrix on new data
    covs.pcox <- folds[[i]]$model$covariate$base
    testing.x.covs <- surv.new %>%
      select(all_of(covs.pcox))
    
    # Mean imputation based on test data
    # No information leakage
    testing.x.covs.imputed <- Imputate.x.mean(testing.x.covs)
    # Convert data to model matrix
    testing.x.mat <- model.matrix(~ ., data = testing.x.covs.imputed)
    # Observed response
    testing.y <- survival::Surv(
      time = surv.new$time,
      event = surv.new$event,
      type = "right"
    )
    # Step 2 - compute linear predictor
    linpred <- predict(
      folds[[i]]$model$cvfit, # Fitted "cv.glmnet" object
      newx = testing.x.mat, # Matrix of new values for x at which predictions are to be made. Must be a matrix
      s = "lambda.min",
      type = "link" # Type "link" (default) returns x^T \beta
    )
    
    # -----------------------------------------------------------------------------------
    # General - Compute tdROC and tdAUC, c-index
    # -----------------------------------------------------------------------------------
    res.tdauc <- Evaluate_tdauc(surv.new, linpred, T.start, deltaT)
    res.c.index <- survcomp::concordance.index(
      x = linpred, # vector of risk predictions
      surv.time = surv.new$time, # vector of event times
      surv.event = surv.new$event, # vector of event occurence indicators
      method = "noether" # conservative, noether or name (see paper Pencina et al. for details)
    )
    # -----------------------------------------------------------------------------------
    
    # -----------------------------------------------------------------------------------
    # Store results
    # -----------------------------------------------------------------------------------
    # Performance
    folds.eval[[i]]$perf <- list(
      landmark = T.start,
      deltaT = deltaT,
      c.index = res.c.index$c.index,
      tdauc = res.tdauc$tdauc,
      tp = res.tdauc$tp,
      fp = res.tdauc$fp)
    # Carry over model information in case training model is not kept
    folds.eval[[i]]$model.info <- list(
      name = folds[[i]]$model$name,
      hyperparam = folds[[i]]$model$hyperparam,
      covariate = folds[[i]]$model$covariate,
      training.time = folds[[i]]$model$training.time
    )
  }
  # -----------------------------------------------------------------------------------
  # Save evaluation result
  folder <- "./output/eval_"
  path.eval <- paste0(folder, model.name, "_seed", seed, ".RData")
  
  save(folds.eval, file = path.eval)
  print(paste("Performance in folds saved:", path.eval))
}