on.alice <- F

# Set folder path on ALICE
folder.alice <- "/data1/s2887592/MFPCCox"
if (on.alice) setwd(folder.alice)

# -----------------------------------------------------------------------------------
# Dependencies
library(tidyverse)

# Brier score
library(survival)
library(pec)

# Common
source("function_utility.R")
source("function_evaluation.R")

# For MFPCCox only
source("function_MFPCCox.R") # From Kan Li github
source("function_MFPCCox_exp.R")

select <- dplyr::select
# -----------------------------------------------------------------------------------
# Set up experiment here!
# -----------------------------------------------------------------------------------
n_fold <- 10 # Cross validation
n_RCV <- 10 # Repeated CV, set to 1 to single CV
T_LMs <- c(2) # Vector of landmark times

seeds <- 721:(721+n_RCV-1) # Seeds for RCV

# -----------------------------------------------------------------------------------
# Set model hyperparam


method <- "MFPCCox"
  
pve <- 0.7 # Hyperparam to choose number of pc
nbasis <- 3 # Mean function

method <- paste(c(
  method,
  paste0("pve", as.character(pve*100)), 
  paste0("nbasis", as.character(nbasis))),
  collapse = "_")


# -----------------------------------------------------------------------------------
# Set data param

set_scenario <- "scenario1" # Determine how many longitudinal covariates to use

# Baseline covariates (not time-varying in adnimerge)
baseline.covs <- c("AGE", "PTGENDER", "PTEDUCAT", "status.bl", "APOE4") # b5 
#  baseline.covs <- c("AGE", "PTGENDER", "PTEDUCAT", "APOE4") # b4 = drop baseline diagnosis

is_transformed <- "transformed" # Transform covariates to reduce skewness
is_scaled <- "scaled" # Set to "scaled" to scale covariates; set to "notScaled" for original

# -----------------------------------------------------------------------------------
# Train test loops
# -----------------------------------------------------------------------------------
# Outer loop - landmark time
for (T.start in T_LMs) {
  #for (T.start in c(1, 2, 3)) {
  print(paste("Start experiment for landmark time T_LM:", T.start))
  # -----------------------------------------------------------------------------------
  # Set identifier
  landmark <- paste0("lm", T.start)
  n_basecov <- paste0("b", length(baseline.covs))
  
  model.hyperparam <- list(
    method = method,
    pve = pve,
    nbasis = nbasis,
    set_scenario = set_scenario,
    landmark = landmark,
    is_scaled = is_scaled,
    is_transformed = is_transformed
  )

  hyperparam <- paste(c(set_scenario, n_basecov, is_transformed, is_scaled), collapse = "_") # Use hyperparam to describe model
  model.name <- paste(c(method, landmark, hyperparam), collapse = "_")
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
  
  # -----------------------------------------------------------------------------------
  # MFPCCox only - Format long data into 3-dim array
  # -----------------------------------------------------------------------------------
  # remove "Fusiform", "Ventricles", "WholeBrain" due to error => relax pve to 0.7 can overcome problem
  # Error in .PACE(X = funDataObject@argvals[[1]], funDataObject@X, Y.pred = Y.pred, : 
  # Measurement error estimated to be zero and there are fewer observed points than PCs; scores cannot be estimated.
  
  # set longitudinal covariates
  #y.names.literature <- c("ADAS13", "MMSE", "RAVLT.learning", "RAVLT.immediate", "FAQ") # longitudinal covariates in literature
  #y.names # Common to pencal
  
  # Use candidate long covariates as basis, then drop covariates that reported error
  # identical over all folds, because the candidate is based on analysis of missing proportions on full data
  
  # time variable for longitudinal variable y
  #y.t <- "Years.bl" # timestamp column name for long covariates
  #y.t <- "M" # use months from baseline instead of years from baseline, less possible values in obstime
  y.t <- "Y" # generic, use years from baseline, rounded
  
  # To prevent a potential bug that the test set has a observed time different from the time domain
  # But it also makes it less robust to new data?
  obstime <- sort(unique(data.long[, y.t])) # get unique obs timestamp in long data
  argvals <- obstime / max(obstime) # scale obstime to [0,1] for uPACE
  
  subject.id <- data.surv$id # subject ids in full data
  nPat <- length(subject.id) # number of subjects in full data
  
  # -----------------------------------------------------------------------------------
  # Middle loop - seed for repeated CV
  for (seed in seeds) { # Vary seed for repeated CV
    #  for (seed in 721:721) { # Only run once
    
    print(paste("This CV will perform train-test split using seed", seed))
    # -----------------------------------------------------------------------------------
    
    # -----------------------------------------------------------------------------------
    # Obtain a list of folds to initialize cross validation
    # Involves:
    # - Select candidate longitudinal covariates based on missingness
    # - Create folds based on stratified train-test split for n-fold CV
    
    # also store scaling table
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
    # [Warning] MFPCCox only - a temporary fix to remove some covariates
    # -----------------------------------------------------------------------------------
    candidate.long.covs <- folds[[1]]$candidate.long.covs
    y.names <- candidate.long.covs[!(candidate.long.covs %in% c("Fusiform", "Ventricles", "WholeBrain", "ICV", "LDELTOTAL"))] 
    #y.names <- vars_long
    n.y <- length(y.names) # number of long covariates
    # -----------------------------------------------------------------------------------
    
    print(paste("Begin training for model:", model.name))
    # -----------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------
    rm("folds")
    rm("folds.eval")
    
    print(path.template)
    load(file = path.template) # Load folds template
    
    folds.eval <- vector(mode = "list", length = n_fold)
    # -----------------------------------------------------------------------------------
    # Fit models in CV loop
    # -----------------------------------------------------------------------------------
    for (i in 1:n_fold) {
      #    for (i in c(1)) { # For debug, run single fold only
      print("---------------------------------------------------------------------------------------------------")
      print(paste("Seed", seed, "- Start training in fold", i))
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
      
      training.surv <- training.surv %>%
        select(all_of(c("id", "time", "event", folds[[i]]$baseline.covs, folds[[i]]$candidate.long.covs)))
      training.long <- training.long %>%
        select(all_of(c("id", "time", "event", "Years.bl", "age.fup", folds[[i]]$baseline.covs, folds[[i]]$candidate.long.covs)))
      
      # -----------------------------------------------------------------------------------  
      # MFPCCox - Initialize multivar array
      # -----------------------------------------------------------------------------------
      # The scaling sould be done within the fold, because the scaling table change on the training data
      if (is_scaled == "scaled") {
        print("[Reminder] scaling is in effect.")
        long.all <- Scale_covariates(data.long, folds[[i]]$scaling_table)
      } else {
        print("[Reminder] no scaling has been done")
        long.all <- data.long
      }
      
      # This multivar array will be shared between train set and test set by indexing
      multivar <- Convert_long_to_mvarray(
        long = long.all, # Full long data, scaled or not scaled
        y.t = y.t, # Time variable used to prepare multivar
        obstime = obstime, # Vector of time to prepare multivar
        subject_id = subject.id, # Provide the order of subject ids in multivar
        n_subject = nPat, 
        y.names = y.names, # Long covariates for multivar
        n.y = n.y)
      
      # -----------------------------------------------------------------------------------
      # MFPCCox - fit model
      # -----------------------------------------------------------------------------------
      res <- Train_MFPCCox(
        training.surv = training.surv,
        multivar = multivar, # Converted long data
        subject.id = subject.id, # Vector of subject ids corresponding to multivar
        y.names = y.names, # Candidate long covariates
        baseline.covs = baseline.covs,
        argvals = argvals, # Scaled time domain
        pve = pve, # Hyperparam to choose number of pc
        nbasis = 3 # Number of basis for mean function
      )
      
      # -----------------------------------------------------------------------------------
      # Store results
      # -----------------------------------------------------------------------------------
      folds[[i]]$model <- list(
        name = model.name, 
        hyperparam = model.hyperparam,
        covariate = list(
          base = folds[[i]]$baseline.covs,
          long = folds[[i]]$candidate.long.covs),
        mfpccox = res$mfpccox,
        mfpca.train = res$mfpca.train,
        phi.train = res$phi.train,
        npc.train = res$npc.train,
        y.names = y.names,
        subject.id = subject.id, # vector of subject id corresponding to multivar array
        obstime = obstime,
        argvals = argvals,
        multivar = multivar,
        training.time = res$runtimes
      )
      print("---------------------------------------------------------------------------------------------------")
      print(paste("Seed", seed, "- Start testing in fold", i))
      # -----------------------------------------------------------------------------------
      # General - Subset subjects for fold i
      # -----------------------------------------------------------------------------------
      # tmp <- Get_train_test_data(
      #   data.surv = data.surv, 
      #   data.long = data.long,  
      #   ids.test = folds[[i]]$ids.test, 
      #   is_scaled = is_scaled, 
      #   scaling_table = folds[[i]]$scaling_table
      # )
      
      surv.new <- tmp$testing.surv
      long.new <- tmp$testing.long
      
      surv.new <- surv.new %>%
        select(all_of(c("id", "time", "event", 
                        folds[[i]]$baseline.covs, 
                        folds[[i]]$candidate.long.covs)))
      long.new <- long.new %>%
        select(all_of(c("id", "time", "event", "Years.bl", "age.fup", 
                        folds[[i]]$baseline.covs, 
                        folds[[i]]$candidate.long.covs)))
      
      # -----------------------------------------------------------------------------------
      # MFPCCox specific landmarking process
      # -----------------------------------------------------------------------------------
      # -----------------------------------------------------------------------------------
      # Subset test data from multivar array
      is_test <- folds[[i]]$model$subject.id %in% surv.new$id # Get subjects at-risk that are event-free at landmark time, t  
      tmp.data <- multivar[is_test, , ] # subset longitudinal outcomes for test set
      
      # Need multivar.train for fold i in UFPCA step
      is_train <- folds[[i]]$model$subject.id %in% training.surv$id # Get (row indices of) subjects in train set
      multivar.train <- multivar[is_train, , ] # Subset for train set
      
      # -----------------------------------------------------------------------------------  
      # MFPCCox - Evaluate model
      # ----------------------------------------------------------------------------------- 
      # -----------------------------------------------------------------------------------  
      # MFPCCox - uFPCA
      # ----------------------------------------------------------------------------------- 
      # univariate FPC 
      Xi.test <- NULL # Xi: FPC scores
      
      for(p in 1:n.y){ # for each longitudinal covariate
        
        print(paste("Computing score for ", folds[[i]]$model$y.names[[p]]))
        
        npc.trained <- folds[[i]]$model$npc.train[[p]]
        
        # estimated trajectories based on a truncated Karhunen-Loeve representation on pred data
        tmp.ufpca <- uPACE(
          testData = multivar.train[, , p], # Specific to train set in each fold
          domain = folds[[i]]$model$argvals, # Should be set of constant
          predData = tmp.data[, , p], 
          nbasis = nbasis,
          #pve = pve,
          npc = npc.trained)
        
        Xi.test <- cbind(Xi.test, tmp.ufpca$scores) # dynamic FPC scores for test subjects 
      }
      # -----------------------------------------------------------------------------------  
      # MFPCCox - MFPCA
      # ----------------------------------------------------------------------------------- 
      # estimate MFPC scores for test subjects
      rho.test <- mfpca.score(Xi.test, folds[[i]]$model$mfpca.train$Cms)
      #tmp.surv.data$rho <- rho.test  
      
      # -----------------------------------------------------------------------------------
      # MFPCCox - Compute the linear predictor
      # ----------------------------------------------------------------------------------- 
      tmp.rho <- data.frame(rho.test)
      rho.names <- paste0("rho", 1:ncol(rho.test))
      names(tmp.rho) <- rho.names
      tmp.surv.data <- data.frame(surv.new, tmp.rho)
      
      X.orig <- tmp.surv.data %>%
        select(all_of(baseline.covs), all_of(rho.names))
      
      linpred <- predict(
        object = folds[[i]]$model$mfpccox, # Fitted "coxph" object
        newdata = X.orig, # new values for x at which predictions are to be made
        type = "lp" # Type "link" (default) returns x^T \beta
      )
      # -----------------------------------------------------------------------------------
      
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
      # MFPCCox Specific - Compute Brier score (work in progress!!!!!!!!!!!!!!!!!!!!!!!)
      # -----------------------------------------------------------------------------------
      
      # pec will return error if the prediction time exceeds the latest survival time in test set
      T.max.test <- floor(max(surv.new$time))
      pred.times <- c(T.start, deltaT)
      pred.times <- pred.times[pred.times <= T.max.test]
      pec.times <- deltaT
      pec.times <- pec.times[pec.times <= T.max.test]
      
      # MFPCCox uses coxph, can feed the trained model to pec, given that the coxph call above has argument X=TRUE to retain X matrix in output
      brier <- tryCatch({
        res.bs <- pec::pec(
          # A matrix with predicted probabilities, dimension of n subjects by m times 
          object = list("model" = folds[[i]]$model$mfpccox),
          formula = Surv(time, event) ~ AGE,
          data = tmp.surv.data, # For computing IPCW
          exact = FALSE, # Do not predict at event times
          times = pec.times, 
          #times = 0:15, 
          cens.model = "cox", # Method for estimating inverse probability of censoring weights:
          splitMethod = "none",
          B = 0,
          verbose = TRUE
        )
        # Return the Brier score evaluated
        res.bs$AppErr$model[-1]
      }, error = function(e) {
        message(e)
        return(NA)
      }, finally = {
      })
      
      # -----------------------------------------------------------------------------------      
      # -----------------------------------------------------------------------------------
      # Store results
      # -----------------------------------------------------------------------------------
      # Performance
      folds.eval[[i]]$perf <- list(
        landmark = T.start,
        deltaT = deltaT,
        brier = brier,
        c.index = res.c.index$c.index,
        tdauc = res.tdauc$tdauc,
        tp = res.tdauc$tp,
        fp = res.tdauc$fp
      )
      # Carry over model information in case training model is not kept
      folds.eval[[i]]$model.info <- list(
        name = folds[[i]]$model$name,
        hyperparam = folds[[i]]$model$hyperparam,
        covariate = folds[[i]]$model$covariate,
        training.time = folds[[i]]$model$training.time
      )
      # -----------------------------------------------------------------------------------
    }
    # -----------------------------------------------------------------------------------
    # Save evaluation result after train test after CV | seed
    folder <- "./output/eval_"
    path.eval <- paste0(folder, model.name, "_seed", seed, ".RData")
    
    save(folds.eval, file = path.eval)
    print(paste("Performance in folds saved:", path.eval))
  }
}