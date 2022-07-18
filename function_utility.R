# ------------------------------------------------------------------------------------------------
# Contains utility functions for setting up experiment
# ------------------------------------------------------------------------------------------------
require(tidyverse)
require(foreach)
require(ggpubr)
require(moments)
# ------------------------------------------------------------------------------------------------
# Initialize folds list for cross validation
# ------------------------------------------------------------------------------------------------
Initialize_exp <- function(
    data.surv,
    data.long,
    baseline.covs, 
    vars_not_long,
    set_scenario,
    n_fold, 
    seed) {
  
  # -----------------------------------------------------------------------------------
  # Select candidate longitudinal covariates based on missingness
  # -----------------------------------------------------------------------------------
  # Determine a set of longitudinal covariates suitable for model
  # Set set_scenario to "literature" to use a manually chosen set of long covairates i.e. vars_literature
  # Options available: "scenario1", "scenario2", "literature", "scenario0"(not tested)

  # Ignore vars_literature if set_scenario is set to TRUE
  vars_literature <- c("ADAS13", "MMSE", "RAVLT.immediate", "RAVLT.learning", "FAQ")
  # This is a subset of longitudinal covariates chosen according to previous literature
  
  
  # Step 1: compute missingness statistics on full data
  # to filter candidate longitudinal covariates automatically 
  missing_proportions <- Compute_missing_proportions(
    data.long = data.long,
    vars_filter = vars_not_long, 
    is_include = FALSE
  )
  
  # Step2: determine a set of longitudinal covariates suitable for model
  # controlled by a pre-defined limit, hence scenario setting
  
  # Set cut-off for maximum proportion of subjects without any values for any covariates
  # Set to 0 to ONLY consider covariates with subjects containing at least 1 measurement
  missing_proportion_limit <- NULL
  
  if (set_scenario == "scenario1") {
    missing_proportion_limit <- 0 
  } else if (set_scenario == "scenario2") {
    missing_proportion_limit <- 0.1
  }

  # Set y.names
  if (set_scenario == "scenario1" | set_scenario == "scenario2") {
    # Automatic selection based on criteria
    # Exclude longitudinal covariates that exceeded the missing_proportion_limit
    vars_missing <- Filter_vars_candidate(
      missing_proportions,
      missing_proportion_limit)
    
    # Compute union of variables to be excluded in LMMs
    # Ignore too many missing or irrelevant or baseline covariates
    vars_exclude <- Reduce(
      union, list(vars_not_long, vars_missing))
    
    y.names <- names(data.long)[!(names(data.long) %in% vars_exclude)] %>%
      sort()

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
    print("[Remind] set_scenario is set to literature to use pre-defined long covariates")
    y.names <- vars_literature
  } else if (set_scenario == "scenario0") {
    print("[Remind] set_scenario is set to scenario0 to not use long covariates")
    y.names <- NULL
  } else {
    stop("Undefined scenario name. Please check input set_scenario.")
  }
  
  # Print long covariates for model
  print("[Report] The following longitudinal covariate(s) will be used for modeling:")
  print("---------------------------------------------------------------------------------------------------")
  print(y.names)
  print("---------------------------------------------------------------------------------------------------")
  print(paste("[Count] final longitudinal covariates =", length(y.names)))
  
  # View(missing_proportions)
  
  # -----------------------------------------------------------------------------------
  # Create folds based on stratified train-test split for n-fold CV
  # -----------------------------------------------------------------------------------
  # folds is a list containing fitted models and train test ids
  # initialize folds here empty list of n_fold size
  # also contains the ids.test to ensure same subset of subjects are used for training and testing
  folds <- Create_folds(data.surv, n_fold, seed)
  
  # -----------------------------------------------------------------------------------
  # Prepare scaling table and store experiment parameters
  
  # Select variables required scaling (numeric), include baseline and long covariates, exclude categorical
  vars_scale <- select_if(
    data.long %>% 
      select(c(all_of(baseline.covs), all_of(y.names))) %>%
      select(-c(AGE, APOE4)), # Exclude categorical
    is.numeric) %>% 
    names()
  

  print("---------------------------------------------------------------------------------------------------")
  print("If is_scaled is set to scaled,")
  print("the following covariate(s) will be centered and scaled:")
  print("---------------------------------------------------------------------------------------------------")
  print(vars_scale)
  print("---------------------------------------------------------------------------------------------------")

  
  # -----------------------------------------------------------------------------------
  # Compute scaling parameters for each fold
  # Store metadata for reuse in future
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
      filter(!(id %in% folds[[i]]$ids.test))
    training.long <- data.long %>% 
      filter(!(id %in% folds[[i]]$ids.test))
    # -----------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------
    # Compute scaling parameters using TRAINING data
    scaling_table <- Compute_scaling_table(training.long, vars_scale)
    # -----------------------------------------------------------------------------------
    
    # -----------------------------------------------------------------------------------
    # Store experiment data and metadata
    # -----------------------------------------------------------------------------------
    folds[[i]]$scenario <- set_scenario
    folds[[i]]$missing_proportions <- missing_proportions
    folds[[i]]$missing_proportion_limit <- missing_proportion_limit
    
    folds[[i]]$seed.split <- seed
    
    folds[[i]]$scaling_table <- scaling_table
    
    folds[[i]]$baseline.covs <- baseline.covs
    folds[[i]]$candidate.long.covs <- y.names
  }
  
  return(folds)
}
# End of Initialize_exp function
# ------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------
# Generate fold ids for train-test split
# ------------------------------------------------------------------------------------------------
Create_folds <- function(surv, n_fold, seed) {
  
  set.seed(seed)
  
  folds <- vector(mode = "list", length = n_fold)
  
  fold_size <- round(1 / n_fold * length(surv$id), 0)
  
  ids.remaining <- as.numeric(surv$id) # Initialize
  for (i in 1:n_fold) {
    j <- n_fold - (i - 1) # Count down
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
    folds[[i]]$seed <- seed
  }
  
  return (folds)
}
# ------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------
# Check results from train-test split
# ------------------------------------------------------------------------------------------------
Check_folds <- function(surv, folds) {
  print("--------------------------------------------------------------------------------")
  print("Check stratification by events in test set for each fold:")
  print("--------------------------------------------------------------------------------")
  for (i in 1:length(folds)) {
    print("--------------------------------------------------------------------------------")
    print(paste("Checking frequency in fold", i))
    print("--------------------------------------------------------------------------------")
    ids.test <- folds[[i]]$ids.test
    ## -------------------------------------------------------
    events.test <- surv %>% 
      filter(id %in% ids.test) %>%
      select(event)
    
    df.freq.test <- data.frame(events.test %>% table()) %>% 
      mutate(Ratio = Freq / sum(Freq))
    
    colnames(df.freq.test) <- c("Event", "Frequency", "Ratio")
    
    print("Test set")
    print(as.matrix(df.freq.test))
    ## -------------------------------------------------------
    events.train <- surv %>%
      filter(!(id %in% ids.test)) %>% # Complement of test set
      select(event)
    
    df.freq.train <- data.frame(events.train %>% table()) %>% 
      mutate(Ratio = Freq / sum(Freq))
    
    colnames(df.freq.train) <- c("Event", "Frequency", "Ratio")
    
    print("Train set")
    print(as.matrix(df.freq.train))
    print("--------------------------------------------------------------------------------")
  }
  print("--------------------------------------------------------------------------------")
}
# ------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------
# Return a table of missing proportions on subjects without any observations per longitudinal covariate
# ------------------------------------------------------------------------------------------------
# Smaller is better, proportion_y_missing = 0 implies all subjects have at least one measurement
Compute_missing_proportions <- function(data.long, vars_filter, is_include = FALSE) {
  
  if (is_include == TRUE) { # Include vars supplied by vars_filter
    data.long <- data.long %>%
      select(all_of(vars_filter))
  } else if (is_include == FALSE) { # Default is to supply the vars to ignore for filtering
    data.long <- data.long %>%
      select(-all_of(vars_filter), "id")
  }
  
  missing_proportions <- data.long %>%
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
  
  # Count repeated observations per covariate
  obs_counts <- data.long %>%
    group_by(id) %>%
    dplyr::summarise(across(everything(), 
                            function(x) sum(!is.na(x)))) %>% # Within subject, non-NA value
    dplyr::select(-id) %>% # Drop id
    dplyr::summarise(across(everything(),
                            function(x) mean(x))) %>% # Proportion of subjects without any y values
    pivot_longer(cols = everything(), 
                 names_to = "var_name",
                 values_to = "avg_obs_per_subject")
  
  missing_proportions <- missing_proportions %>%
    left_join(obs_counts, by = "var_name")
  
  return(missing_proportions)
}
# ------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------
# Return a vector of covariate names that exceeded the user defined missing proportion limit
# ------------------------------------------------------------------------------------------------
Filter_vars_candidate <- function(missing_proportions, missing_proportion_limit) {
  # Filter covariates with missing proportions exceeding pre-defined limit
  vars_missing <- missing_proportions %>%
    filter(proportion_y_missing > missing_proportion_limit) %>% # Higher than limit is undesirable
    select(var_name) %>%
    unlist(use.names = FALSE) %>%
    sort()
  
  return(vars_missing)
}
# ------------------------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------------------------
# Return a list of train and test data
# ------------------------------------------------------------------------------------------------
Get_train_test_data <- function(
    data.surv, # Full surv data 
    data.long, # Full long data
    ids.test, # Subject ids for test set
    is_scaled, # Set to "scaled" to scale numeric covariates
    scaling_table # Scaling parameters derived from training data
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
  # Scaling - update train test set if scaling is required
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

# ------------------------------------------------------------------------------------------------
# Apply transformation to longitudinal data
# ------------------------------------------------------------------------------------------------
Transform_covariates <- function(
    data.long,
    y.names, # Vars to transform 
    threshold.skew, 
    threshold.sym) {
  
  # Compute statistics for original data
  y.stats <- foreach(var_name = y.names, .combine = rbind) %do% {
    c(skewness(data.long[, var_name], na.rm = TRUE), fivenum(data.long[, var_name], na.rm = TRUE))
  } %>% 
    data.frame(row.names = NULL)
  
  colnames(y.stats) <- c("skewness", "minimum", "lower_hinge", "median", "upper_hinge", "maximum")
  y.stats$var_name <- y.names
  
  y.stats <- y.stats %>% 
    mutate(symmetry = (upper_hinge - median) / (median - lower_hinge)) %>%
    select(var_name, skewness, symmetry, everything())
  # Note: Symmetry ratio constructed according to Sec 4.2, Fox, Applied Regression Analysis and Generalised Linear Models
  
  # If variable contains value <= 0, add a positive constant `start` before power transformation
  # To ensure that the power transformation will be monotone i.e. preserve order after transformation
  # To ensure log transformation has defined input domain
  y.stats <- y.stats %>% 
    mutate(is_neg = minimum < 0,
           start = -floor(minimum),
           is_sym_over = symmetry > 1 + threshold.sym,
           is_sym_under = symmetry < 1 - threshold.sym,
           is_skew_pos = skewness > 0 + threshold.skew,
           is_skew_neg = skewness < 0 - threshold.skew
    ) %>%
    select(var_name, is_neg, start, is_sym_over, is_sym_under, is_skew_pos, is_skew_neg, everything())
  

  # Transformation

  data.long.transformed <- data.long # Copy original data
  
  vec.transformation <- vector(mode = "character", length = length(y.names))
  
  for (i in 1:length(y.names)) {
    
    var_name <- y.names[i]
    #  var_name <- "ADAS11" # Debug
    
    y.raw <- data.long[, var_name]
    y.stats_extract <- y.stats[y.stats$var_name == var_name, ]
    transformation <- var_name # String to describe the transformation
    
    raw.skew <- y.stats_extract$skewness
    raw.min <- y.stats_extract$minimum
    raw.sym <- y.stats_extract$symmetry
    
    
    # Add a positive constant if y <= 0 before transformation
    if (raw.min < 0) {
      start <- y.stats_extract$start
      y.raw <- y.raw + start
      transformation <- paste0(transformation, "+", as.character(start))
    } else if (raw.min == 0) {
      start <- 1
      y.raw <- y.raw + start
      transformation <- paste0(transformation, "+", as.character(start))
    }
    
    # Transform covariate
    if (raw.skew > threshold.skew) { # Positive skew
      y.transformed <- log10(y.raw) # Down the ladder
      transformation <- paste0("log10(", transformation, ")")
    } else if (raw.skew < -threshold.skew) { # Negative skew
      y.transformed <- y.raw^3 # Up the ladder
      transformation <- paste0("(", transformation, ")^3")
    } else { # Within threshold
      y.transformed <- y.raw # Keep original
      transformation <- NA
    }
    # Store the transformed data
    data.long.transformed[, var_name] <- y.transformed
    vec.transformation[i] <- transformation
  }
  
  
  # Recalculate the statistics after transformation
  y.stats.transformed <- foreach(var_name = y.names, .combine = rbind) %do% {
    c(skewness(data.long.transformed[, var_name], na.rm = TRUE), fivenum(data.long.transformed[, var_name], na.rm = TRUE))
  } %>% 
    data.frame(row.names = NULL)
  
  colnames(y.stats.transformed) <- c("skewness_t", "minimum_t", "lower_hinge_t", "median_t", "upper_hinge_t", "maximum_t")
  y.stats.transformed$var_name <- y.names
  
  y.stats.transformed <- y.stats.transformed %>% 
    mutate(symmetry_t = (upper_hinge_t - median_t) / (median_t - lower_hinge_t)) %>%
    select(var_name, skewness_t, symmetry_t, everything())
  
  y.stats.transformed <- y.stats.transformed %>%
    left_join(y.stats, by = "var_name")
  
  y.stats.transformed$transformation <- vec.transformation
  
  # Organize columns
  y.stats.transformed <- y.stats.transformed %>%
    select(var_name, transformation, everything())
  
  return(list(
    data.long.transformed = data.long.transformed,
    summary = y.stats.transformed
  ))
  # End of Transform_covariates
}
# ------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------
# Plot to compare transformation
# Output: Save histogram before and after trasnformation to ./output/transformation/
# ------------------------------------------------------------------------------------------------
Plot_transformation <- function(
    data.long, # Before
    data.long.transformed, # After
    summary, # Transformation summary
    folder = "./output/transformation/") {
  
  y.names <- summary$var_name
  
  for (var_name in y.names) {
    
    stats_extract <- summary[summary$var_name == var_name, ]
    
    if (!is.na(stats_extract$transformation)) { # If transformation took place
      
      # Extract data
      y.raw <- data.long[, var_name]
      y.transformed <- data.long.transformed[, var_name]
      
      # Extract stats
      raw.min <- stats_extract$minimum
      raw.skew <- stats_extract$skewness
      raw.sym <- stats_extract$symmetry
      
      transformed.min <- stats_extract$minimum_t
      transformed.skew <- stats_extract$skewness_t
      transformed.sym <- stats_extract$symmetry_t
      
      # Compare
      g.raw <- ggplot(data = data.frame(y.raw)) + 
        geom_histogram(aes(y.raw), na.rm = TRUE, bins = 50) +
        geom_vline(xintercept = raw.min, color = "blue") + # Minimum
        xlab(var_name) +
        labs(
          title = paste(var_name, "(before transform)"),
          subtitle = paste("skew:", round(raw.skew, 2), "sym:", round(raw.sym, 2))
        )  
      
      g.transformed <- ggplot(data = data.frame(y.transformed)) + 
        geom_histogram(aes(y.transformed), na.rm = TRUE, bins = 50) +
        geom_vline(xintercept = transformed.min, color = "blue") + # Minimum
        xlab(stats_extract$transformation) +
        labs(
          title = paste(var_name, "(after transform)"),
          subtitle = paste("skew:", round(transformed.skew, 2), "sym:", round(transformed.sym, 2))
        )
      
      # Combine plots and save
      g <- ggarrange(plotlist = list(g.raw, g.transformed))
      fn <- paste0("output_transf_hist_", var_name, ".png")
      ggsave(filename = fn, plot = g, path = folder, width = 8, height = 4)
    }
  }
}
# ------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------
# Update the surv data according to landmark time
# Change the time-varying variables in data.surv
# ------------------------------------------------------------------------------------------------
Get_last_nonNA <- function(x) {
  if (all(is.na(x))) {
    NA
  } else {
    tail(x[!is.na(x)], 1)
  }
}

Update_surv_at_landmark <- function(surv, long, y.names, use_baseline) {
  # Note: set use_baseline to TRUE for "glmnet-bl" methods
  
  surv.mod <- surv %>% # Copy original surv but drop the variables to be replaced
    select(-all_of(y.names))
  
  # Use baseline values from data.long instead of latest observed up to landmark
  if (use_baseline) { 
    long <- long %>%
      filter(VISCODE == "bl")
  }
  
  df.last_nonNA <- long %>% 
    group_by(id) %>%
    summarise(across(.cols = y.names, .fns = Get_last_nonNA))

  surv.mod <- surv.mod %>%
    left_join(df.last_nonNA, by = "id") %>%
    select(names(surv))
  
  return(surv.mod)
}
# ------------------------------------------------------------------------------------------------
