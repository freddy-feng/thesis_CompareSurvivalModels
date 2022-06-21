library(tidyverse)

source("function_utility_exp.R")
source("function_evaluation.R")

source("function_pencal_exp.R")

select <- dplyr::select

for (seed in 721:730) {
  # -----------------------------------------------------------------------------------
  # Set experiment parameters here!
  # -----------------------------------------------------------------------------------
  n_fold <- 10 # Cross validation
  seed <- seed # Vary this for repeated CV
  set_scenario <- "scenario2" # Determine how many longitudinal covariates to use
  is_transformed <- "transformed"
  
#  baseline.covs <- c("AGE", "PTGENDER", "PTEDUCAT", "status.bl", "APOE4") # Baseline covariates
  baseline.covs <- c("AGE", "PTGENDER", "PTEDUCAT", "APOE4") # Baseline covariates
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
  T.start <- 2
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
  
  # Initialize fold
  print(path.template)
  load(file = path.template)
  
  # Model param
  penalty <- "ridge"
  is_scaled <- "scaled" # Set to "scaled" to scale covariates for LMMs
  is_standardized <- "std" # Refer to step 3 of pencal, standardizing the random effects summary
  is_transformed <- is_transformed
  landmark <- paste0("lm", T.start)
  
  hyperparam <- paste(c(landmark, is_transformed, is_scaled, is_standardized, penalty), collapse = "_") # Use hyperparam to describe model
  
  print(hyperparam)
  
  # -----------------------------------------------------------------------------------
  # Fit models in CV loop
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
    # pencal - Fit model
    # -----------------------------------------------------------------------------------
    res <- run_prc_steps(
      long.data = training.long, 
      surv.data = training.surv,
      baseline.covs = folds[[i]]$baseline.covs,
      y.names = folds[[i]]$candidate.long.covs, # where did you specify baseline covariates
      n.boots = 0,
      n.cores = parallel::detectCores(),
      verbose = TRUE,
      penalty = penalty,
      standardize = is_standardized == "std"
    )
    # -----------------------------------------------------------------------------------
    
    # -----------------------------------------------------------------------------------
    # pencal - Store results
    # -----------------------------------------------------------------------------------
    folds[[i]]$pencal <- list(
      step1 = res$step1,
      step2 = res$step2,
      step3 = res$step3,
      runtimes = res$runtimes,
      hyperparam = list(
        penalty = penalty,
        is_scaled = is_scaled,
        is_standardized = is_standardized
      )
    )
    # -----------------------------------------------------------------------------------
  }
  # -----------------------------------------------------------------------------------
  
  # Initialize result
  folds.eval <- vector(mode = "list", length = n_fold)
  
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
    
    surv.new <- tmp$testing.surv
    long.new <- tmp$testing.long 
    
    
    # testing.surv <- tmp$testing.surv
    # testing.long <- tmp$testing.long
    # ----------------------------------------------------------------------------------- 
    #   # Landmarking -only consider subjects at risk at landmark time T.start
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
    #   
    #  print(paste("[Count] subject at risk at landmark", length(ids.at_risk)))
    
    # -----------------------------------------------------------------------------------  
    # pencal - Evaluate model
    # ----------------------------------------------------------------------------------- 
    # Extract fitted models
    step1 <- folds[[i]]$pencal$step1
    step2 <- folds[[i]]$pencal$step2
    step3 <- folds[[i]]$pencal$step3
    
    # Obtain predicted random effect
    # res object comes from fitted pencal
    preds <- survpred_prclmm(
      step1 = step1, 
      step2 = step2, 
      step3 = step3,
      times = deltaT, # Prediction window(s) for surv prob prediction
      new.longdata = long.new, # Long data in test set
      new.basecovs = surv.new, # Surv data in test set
      keep.ranef = TRUE 
    )
    # -----------------------------------------------------------------------------------  
    # Obtain new X for predict on pcox.orig
    # With baseline covariate, may need to expand into without baseline covariate case later
    # -----------------------------------------------------------------------------------
    X0 <- model.matrix(as.formula(step3$call$baseline.covs),
                       data = surv.new)
    # Join the baseline covariates and predicted random effect summary
    pred_ranefs <- preds$predicted_ranefs
    X.orig <- as.matrix(cbind(X0, as.matrix(pred_ranefs)))
    
    # Drop intercept
    contains.int <- "(Intercept)" %in% colnames(X.orig)
    if (contains.int) {
      X.orig <- X.orig[, -1]
    }
    # -----------------------------------------------------------------------------------
    # Compute the linear predictor
    # -----------------------------------------------------------------------------------
    linpred <- predict(
      object = step3$pcox.orig, # Fitted "cv.glmnet" or "cv.relaxed" object
      newx = X.orig, # Matrix of new values for x at which predictions are to be made. Must be a matrix
      s = "lambda.min",
      type = "link" # Type "link" (default) returns x^T \beta
    )
    # -----------------------------------------------------------------------------------
    
    
    
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
  tmp <- paste0("./output/eval_", set_scenario, "_pencal_", hyperparam, "_seed", seed, ".RData")
  save(folds.eval, file = tmp)
  print(paste("folds.eval saved to path:", tmp))
}