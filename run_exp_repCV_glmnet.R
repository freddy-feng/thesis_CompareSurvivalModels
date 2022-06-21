library(tidyverse)

require(glmnet)
require(survival)

source("function_utility_exp.R")
source("function_evaluation.R")

source("function_glmnet_exp.R")

select <- dplyr::select



for (seed in 721:730) {
  # -----------------------------------------------------------------------------------
  # Set experiment parameters here!
  # -----------------------------------------------------------------------------------
  n_fold <- 10 # Cross validation
  seed <- seed # Vary this for repeated CV
  set_scenario <- "scenario2" # Determine how many longitudinal covariates to use
  is_transformed <- "transformed"
  
  baseline.covs <- c("AGE", "PTGENDER", "PTEDUCAT", "status.bl", "APOE4") # Baseline covariates
  #baseline.covs <- c("AGE", "PTGENDER", "PTEDUCAT", "APOE4") # Baseline covariates
  # -----------------------------------------------------------------------------------
  # -----------------------------------------------------------------------------------
  
  # -----------------------------------------------------------------------------------
  # Load data
  # -----------------------------------------------------------------------------------
  # Contains two dataframes df.surv_preds and df.long_censored
  path.data <- "adni_cleaned.RData"
  load(path.data)
  
  # may also reduce the number of columns here to reduce size
  # Note: make sure data.surv and data.long are arranged by id and {id, age.fup}
  # to ensure proper function with pencal package
  data.surv <- df.surv_preds
  
  if (is_transformed == "transformed") {
    data.long <- df.long_censored_transformed
  } else {
    data.long <- df.long_censored
  }
  
  ## Set up for evaluation
  
  # landmark time
  T.start <- 1
  landmark <- paste0("lm", T.start)
  # We can predict on years from (T.start + 1) up to T.max
  
  T.max <- floor(max(data.surv$time)) # Based on last available observation in train set
  # [Warning] forsee a potential bug may happen by chance if the longest observation is in test set, but not train
  
  # Predict 15 years onward from landmark time, but not more than max observed time
  deltaT <- 1:T.max # A vector of prediction times, starting from baseline onward
  deltaT <- deltaT[deltaT > T.start]
  
  # I need to change deltaT later because the notation is inconsistent with the symbol in paper
  # the true delta T should be prediction time - landmark time!!!!!!!
  
  print(paste("Evaluation on time since baseline =", paste(deltaT, collapse = " ")))
  
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
  # Type of TAU, PTAU and ABETA are character, need to handle non-numerical values first. currently excluded
  vars_manual_remove <- c("TAU", "PTAU", "ABETA")
  
  # Exclude irrelevant variables
  vars_irrelevant <- c(
    names(data.long)[grepl(".bl", names(data.long))], # Exclude baselines, Years.bl, Month.bl
    "RID", "time", "event", "status", "DX",
    "VISCODE", "EXAMDATE", "Y", "M", "Month",
    "AGE", "age.fup", "COLPROT", "ORIGPROT", "PTID", "SITE",
    "PTGENDER", "PTEDUCAT", "PTETHCAT", "PTRACCAT", "PTMARRY", "APOE4",
    "FSVERSION", "IMAGEUID", "FLDSTRENG"
  )
  
  vars_ignore <- c(vars_manual_remove, vars_irrelevant) # Variables that will not be considered as long covariates
  
  # Updating the time-varying covariates in data.surv according to landmark time
  # For each subject, the latest observed value of time-varying covariate is used
  # The value can be transformed or not, depending on data.long assignment
  
  # vars_candidate are the set of candidate long covariates BEFORE screening
  vars_candidate <- names(data.surv)[!(names(data.surv) %in% c("id", vars_ignore))]
  
  # -----------------------------------------------------------------------------------
  # Update landmark value
  # Can it be speed up?
  data.surv <- Update_surv_at_landmark(
    surv = data.surv, 
    long = data.long, 
    y.names = vars_candidate)
  
  # -----------------------------------------------------------------------------------
  folds <- Initialize_exp(
    data.surv = data.surv,
    data.long = data.long,
    baseline.covs = baseline.covs,
    vars_ignore = vars_ignore,
    set_scenario = set_scenario,
    n_fold = n_fold,
    seed = seed
  )
  # Get scaling table
  # Note: depends on is_transformed, different set of data.long will be used, hence different scaling factor will be determined
  
  #Check_folds(data.surv, folds) # Checking the split process, stratification
  
  path.template <- paste0("output_folds_template_", set_scenario, "_seed", seed, "_", landmark, "_", is_transformed, ".RData")
  
  save(folds, file = path.template)
  
  print(paste("template of folds saved to path:", path.template))
  # -----------------------------------------------------------------------------------
  rm("folds")
  rm("folds.eval")
  
  print(path.template)
  load(file = path.template)
  
  # Model param
  is_scaled <- "scaled" # Set to "scaled" to scale covariates for LMMs
  is_transformed <- is_transformed
  landmark <- paste0("lm", T.start)
  
  hyperparam <- paste(c(landmark, is_transformed, is_scaled), collapse = "_") # Use hyperparam to describe model
  
  print(hyperparam)
  
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
      baseline.covs.additional = baseline.covs.additional
    )
    
    # -----------------------------------------------------------------------------------
    # Store results
    # -----------------------------------------------------------------------------------
    folds[[i]]$glmnet <- list(
      cvfit = res$cvfit,
      covs.pcox = res$covs.pcox,
      runtimes = res$runtimes,
      is_scaled = is_scaled
    )
    # -----------------------------------------------------------------------------------
  }
  # -----------------------------------------------------------------------------------
  folds.eval <- vector(mode = "list", length = n_fold)
  
  for (i in 1:n_fold) {
    # -----------------------------------------------------------------------------------
    # General - Subset subjects for fold i
    # -----------------------------------------------------------------------------------
    tmp <- Get_train_test_data(
      data.surv, data.long, 
      folds[[i]]$ids.test,
      is_scaled, 
      folds[[i]]$scaling_table
    )
    
    surv.new <- tmp$testing.surv
    long.new <- tmp$testing.long
    
    # testing.surv <- tmp$testing.surv
    # testing.long <- tmp$testing.long
    # ----------------------------------------------------------------------------------- 
    # Landmarking -only consider subjects at risk at landmark time T.start
    #   ids.at_risk <- testing.surv %>% 
    #     filter(time > T.start) %>%
    #     select(id) %>%
    #     unlist(use.names = FALSE)
    #   
    #   surv.new <- testing.surv %>%
    #     filter(id %in% ids.at_risk)
    #   
    #   # Commented out this part to use all longitudinal information
    # #  long.new <- testing.long %>%
    # #    filter(id %in% ids.at_risk) %>%
    # #    filter(Years.bl <= T.start) # Filter long observations after landmark time
    #   long.new <- testing.long %>%
    #     filter(id %in% ids.at_risk)
    
    #  print(paste("[Count] subject at risk at landmark", length(ids.at_risk)))
    
    # -----------------------------------------------------------------------------------  
    # pCox - Evaluate model
    # ----------------------------------------------------------------------------------- 
    # Step 1 - prepare new X
    # Subset columns to construct model matrix on new data
    covs.pcox <- folds[[i]]$glmnet$covs.pcox
    testing.x.covs <- surv.new %>%
      select(all_of(covs.pcox))
    # Mean imputation
    testing.x.covs.imputed <- Imputate.x.mean(testing.x.covs)
    # Model matrix, of dimension n obs x n vars; each row is an observation vector
    testing.x.mat <- model.matrix(~ ., data = testing.x.covs.imputed)
    # Observed response
    testing.y <- survival::Surv(
      time = surv.new$time,
      event = surv.new$event,
      type = "right"
    )
    # Step 2 - compute linear predictor
    linpred <- predict(
      folds[[i]]$glmnet$cvfit, # Fitted "cv.glmnet" object
      newx = testing.x.mat, # Matrix of new values for x at which predictions are to be made. Must be a matrix
      s = "lambda.min",
      type = "link" # Type "link" (default) returns x^T \beta
    )
    
    # -----------------------------------------------------------------------------------
    # General - Compute tdROC and tdAUC, c-index
    # -----------------------------------------------------------------------------------
    res.tdauc <- Evaluate_tdauc(surv.new, linpred, T.start, deltaT)
    res.c.naive <- survcomp::concordance.index(
      x = linpred, # vector of risk predictions
      surv.time = surv.new$time, # vector of event times
      surv.event = surv.new$event, # vector of event occurence indicators
      method = "noether" # conservative, noether or name (see paper Pencina et al. for details)
    )
    # -----------------------------------------------------------------------------------
    
    # -----------------------------------------------------------------------------------
    # Store results
    # -----------------------------------------------------------------------------------
    #  folds.eval[[i]]$ids.test_no_missing <- ids.valid
    folds.eval[[i]]$perf <- list(
      landmark = T.start,
      deltaT = deltaT,
      c.index = res.c.naive$c.index,
      tdauc = res.tdauc$tdauc,
      tp = res.tdauc$tp,
      fp = res.tdauc$fp
    )
    # -----------------------------------------------------------------------------------
  }
  # Save evaluation result
  tmp <- paste0("./output/eval_", set_scenario, "_pCox_", hyperparam, "_seed", seed, ".RData")
  save(folds.eval, file = tmp)
  
  print(paste("Performance in folds saved to path:", tmp))

}