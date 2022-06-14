# ------------------------------------------------------------------------------------------------
# Contains utility function for setting up experiment
# ------------------------------------------------------------------------------------------------
#
require(tidyverse)

Initialize_exp <- function(
    data.surv, data.long, 
    set_scenario,
    n_fold, seed) {
  
  baseline.covs <- c("AGE", "PTGENDER", "PTEDUCAT", "status.bl", "APOE4")
  
  # -----------------------------------------------------------------------------------
  # Prepare candidate longitudinal covariates
  # -----------------------------------------------------------------------------------
  # Set set_scenario to "literature" to use a manually chosen set of long covairates i.e. vars_literature
  # Options available: "scenario1", "scenario2", "literature"

  # Ignore vars_literature if set_scenario is set to TRUE
  vars_literature <- c("ADAS13", "MMSE", "RAVLT.immediate", "RAVLT.learning", "FAQ")
  # This is a subset of longitudinal covariates chosen according to previous literature
  
  # Manual exclusion ----------------------------------------------------------------------------------------
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
  # -----------------------------------------------------------------------------------
  
  ## Filter candidate longitudinal covariates automatically based on missingness criteria
  
  missing_proportions <- Compute_missing_proportions(
    data.long = data.long,
    vars_ignore = c(vars_irrelevant, vars_manual_remove)
  )
  # -----------------------------------------------------------------------------------
  
  # Construct a set of longitudinal candidate covariates
  
  # Set cut-off for maximum proportion of subjects without any values for any covariates
  # Set to 0 to ONLY consider covariates with subjects containing at least 1 measurement
  missing_proportion_limit <- NULL
  if (set_scenario == "scenario1") {
    missing_proportion_limit <- 0 
  } else if (set_scenario == "scenario2") {
    missing_proportion_limit <- 0.1
  }
  
  vars_missing <- Filter_vars_candidate(
    missing_proportions,
    missing_proportion_limit)
  # -----------------------------------------------------------------------------------
  # Compute union of all variables to be excluded in LMMs
  vars_exclude <- Reduce(
    union, list("id", 
                vars_manual_remove, 
                vars_irrelevant,
                vars_missing))
  
  vars_candidate <- names(data.long)[!(names(data.long) %in% vars_exclude)] %>%
    sort()
  # -----------------------------------------------------------------------------------
  # Set y.names
  if (set_scenario == "scenario1" | set_scenario == "scenario2") {
    # Automatic selection based on criteria
    # Exclude longitudinal covariates that exceeded the missing_proportion_limit
    y.names <- vars_candidate
    print("set_scenario is set to TRUE. All available covariates will be considered and screened for missing proportions.")
    print(paste("[Report] Found", length(vars_missing), 
                "covariates with proportion of subjects without any observation exceeded user-defined limit of",
                missing_proportion_limit, "."))
    print("Inpsect `missing_proportions` dataframe to check missing proportion values.")
    print("---------------------------------------------------------------------------------------------------")
    print("The following covariate(s) will NOT be used for LMM:")
    print("---------------------------------------------------------------------------------------------------")
    print(vars_missing)
    print("---------------------------------------------------------------------------------------------------")
  } else if (set_scenario == "literature") {
    print("set_scenario is set to literature")
    y.names <- vars_literature
  } else {
    stop("Undefined sceanrio. Please check set_scenario.")
  }
  
  
  print("The following covariate(s) will be used for LMM:")
  print("---------------------------------------------------------------------------------------------------")
  print(y.names)
  print("---------------------------------------------------------------------------------------------------")
  print(paste("[Count] final longitudinal covariates =", length(y.names)))
  
  # View(missing_proportions)
  
  # -----------------------------------------------------------------------------------
  # Create folds, commonly shared between methods
  # Stratified train-test split for n-fold CV
  # -----------------------------------------------------------------------------------
  # folds is a list containing fitted models and train test ids
  # initialize folds here empty list of n_fold size
  # also contains the ids.test to ensure same subset of subjects are used for training and testing
  folds <- create_folds(data.surv, n_fold, seed)
  #check_folds(data.surv, folds)
  
  # -----------------------------------------------------------------------------------
  # Prepare scaling table and store experiment parameters
  
  vars_scale <- select_if(
    data.long %>% 
      select(c(all_of(baseline.covs), all_of(y.names))) %>%
      select(-c(AGE, APOE4)),
    is.numeric) %>%
    names()
  

  print("---------------------------------------------------------------------------------------------------")
  print("If is_scaled is set to scaled,")
  print("the following covariate(s) will be centered and scaled:")
  print("---------------------------------------------------------------------------------------------------")
  print(vars_scale)
  print("---------------------------------------------------------------------------------------------------")

  
  # -----------------------------------------------------------------------------------
  for (i in 1:n_fold) {
    # -----------------------------------------------------------------------------------
    # Subset data for fold i
    ## Test set
    testing.surv <- data.surv %>% 
      filter(id %in% folds[[i]]$ids.test)
    testing.long <- data.long %>% 
      filter(id %in% folds[[i]]$ids.test)
    
    ## Train set (complement to test set)
    training.surv <- data.surv %>% 
      filter(!(id %in% testing.surv$id))
    training.long <- data.long %>% 
      filter(id %in% training.surv$id)
    # -----------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------
    # Construct scaling table using training data
    scaling_table <- Compute_scaling_table(training.long, vars_scale)
    # -----------------------------------------------------------------------------------
    
    
    # -----------------------------------------------------------------------------------
    # Store experiment data and metadata
    # -----------------------------------------------------------------------------------
    folds[[i]]$scenario <- set_scenario
    folds[[i]]$seed.split <- seed
    
    folds[[i]]$scaling_table <- scaling_table
    
    folds[[i]]$missing_proportions <- missing_proportions
    folds[[i]]$missing_proportion_limit <- missing_proportion_limit
    
    folds[[i]]$baseline.covs <- baseline.covs
    folds[[i]]$candidate.long.covs <- y.names
  }
  
  return(folds)
}
# ------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------
# Generate fold ids for train-test split
# ------------------------------------------------------------------------------------------------

create_folds <- function(surv, n_fold, seed) {
  
  set.seed(seed)
  
  folds <- vector(mode = "list", length = n_fold)
  
  fold_size <- round(1 / n_fold * length(surv$id), 0)
  
  ids.remaining <- as.numeric(surv$id) # Initialize
  for (i in 1:n_fold) {
    j <- n_fold - (i - 1)
    if (j > 1) { # Not last fold
      tmp <- surv %>%
        filter(id %in% ids.remaining) %>% # Excluded subjects from other split
        splitstackshape::stratified(., group = "event", size = 1 / j) 
      # size is set to be proportional to the number of observations per group
      # Store ids
      folds[[i]]$ids.test <- as.numeric(tmp$id)
      # Update remaining subjects for selection
      ids.remaining <- setdiff(ids.remaining, tmp$id)
    } else { # Last fold
      folds[[i]]$ids.test <- ids.remaining
    }
  }
  
  return (folds)
}

# Check fold
check_folds <- function (surv, folds) {
  print("--------------------------------------------------------------------------------")
  print("Check stratification by events in test set for each fold:")
  for (i in 1:length(folds)) {
    print(paste("Fold", i))
    ids.test <- folds[[i]]$ids.test 
    events.test <- surv %>% 
      filter(id %in% ids.test) %>%
      select(event)
    print("Event frequency:")
    print(events.test %>%
            table())
    print("Event proportion:")
    print(events.test %>%
            table() %>%
            prop.table() %>%
            round(digits = 3))
    print("--------------------")
  }
  print("--------------------------------------------------------------------------------")
}
# ------------------------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------------------------
# Summarize tdAUC in different folds from folds.eval and output dataframe
# ------------------------------------------------------------------------------------------------
Summarize.tdAUC <- function(method, path, n_fold, detlaT) {
  # Load results to folds.eval list for corresponding method
  load(path)
  
  # Extract from folds.eval into list of tdAUC
  list.tdauc <- lapply(1:n_fold, function(i) {
    folds.eval[[i]]$perf$tdauc
  })
  
  # Convert list into data frame
  df.tdauc <- data.frame(do.call(rbind, list.tdauc)) # row i = fold i
  colnames(df.tdauc) <- deltaT
  
  # Compute statistics
  df.tdauc <- df.tdauc %>%
    pivot_longer(
      cols = everything(),
      names_to = "prediction_time",
      values_to = "tdAUC"
    ) %>%
    group_by(prediction_time) %>%
    summarise(
      mean = mean(tdAUC),
      ci.upper = quantile(tdAUC, 0.975),
      ci.lower = quantile(tdAUC, 0.025)
    ) 
  
  
  df.tdauc$prediction_time <- as.numeric(df.tdauc$prediction_time)
  df.tdauc <- df.tdauc %>% arrange(prediction_time)
  df.tdauc$method <- method
  
  return(df.tdauc)
}
# ------------------------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------------------------
# Return a table of missing proportions on subjects without any observations per longitudinal covariate
# ------------------------------------------------------------------------------------------------
# Smaller is better, proportion_y_missing = 0 implies all subjects have at least one measurement
Compute_missing_proportions <- function(data.long, vars_ignore) {
  missing_proportions <- data.long %>%
    select(-all_of(vars_ignore)) %>% # Exclude manual exclusion defined above
    group_by(id) %>%
    dplyr::summarise(across(everything(), 
                            function(x) all(is.na(x)))) %>% # Within subject, all y is missing
    dplyr::select(-id) %>%
    dplyr::summarise(across(everything(),
                            function(x) sum(x) / length(unique(data.long$id)))) %>% # Proportion of subjects without any y values
    pivot_longer(cols = everything(), 
                 names_to = "var_name",
                 values_to = "proportion_y_missing") %>%
    arrange(proportion_y_missing)
  return(missing_proportions)
}
# ------------------------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------------------------
# Return a vector of covariate names that exceeded the user defined missing proportion limit
# ------------------------------------------------------------------------------------------------
Filter_vars_candidate <- function(missing_proportions, missing_proportion_limit) {
  ## Filter covariates with missing proportions exceeding pre-defined limit
  vars_missing <- missing_proportions %>%
    filter(proportion_y_missing > missing_proportion_limit) %>% # Higher than limit is undesirable
    select(var_name) %>%
    unlist(use.names = FALSE) %>%
    sort()
  return(vars_missing)
}
# ------------------------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------------------------
# Return a scaled dataframe from data_to_scale using parameters in scaling table
# ------------------------------------------------------------------------------------------------
Scale_covariates <- function(data_to_scale, scaling_table) {
  vars_to_scale <- names(data_to_scale)[names(data_to_scale) %in% scaling_table$vars]
  for (x in vars_to_scale) {
    params <- scaling_table %>%
      filter(vars == x) %>%
      select(mu, sd) %>%
      unlist()
    
    data_to_scale[, x] <- (data_to_scale[, x] - params["mu"]) / params["sd"]
  }
  return(data_to_scale)
}
# ------------------------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------------------------
# Return a dataframe containing the (scalable, numeric) variables names, mu, sd
# ------------------------------------------------------------------------------------------------
Compute_scaling_table <- function(data, vars_scale) {
  scaling_table <- data.frame(do.call(rbind, lapply(vars_scale, function(x) {
    c(mu = mean(data[, x], na.rm = TRUE),
      sd = sd(data[, x], na.rm = TRUE))
  })))
  scaling_table$vars <- vars_scale
  scaling_table <- scaling_table %>%
    select(vars, everything())
  return(scaling_table)
}
# ------------------------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------------------------
# Return a list of train and test data
# ------------------------------------------------------------------------------------------------
Get_train_test_data <- function(
    data.surv, data.long, ids.test,
    is_scaled, scaling_table
    ) {
  ## Test set
  testing.surv <- data.surv %>% 
    filter(id %in% ids.test)
  testing.long <- data.long %>% 
    filter(id %in% ids.test)
  
  ## Train set (complement to test set)
  training.surv <- data.surv %>% 
    filter(!(id %in% ids.test))
  training.long <- data.long %>% 
    filter(!(id %in% ids.test))
  
  # -----------------------------------------------------------------------------------
  # Scaling
  # -----------------------------------------------------------------------------------
  if (is_scaled == "scaled") {
    # Scale using scaling parameters in scaling table
    training.surv <- Scale_covariates(training.surv, scaling_table)
    training.long <- Scale_covariates(training.long, scaling_table)
    # Scaling the testing data using scaling parameters derived from training set | fold i
    testing.surv <- Scale_covariates(testing.surv, scaling_table)
    testing.long <- Scale_covariates(testing.long, scaling_table)
  }
  return(list(
    training.surv = training.surv,
    training.long = training.long,
    testing.surv = testing.surv,
    testing.long = testing.long
  ))
}
# ------------------------------------------------------------------------------------------------


