# Dependencies
library(tidyverse)
library(pencal)
#library(splitstackshape)

set.seed(721)

# Load data
## Contains two dataframes df.surv_preds and df.long_censored
load("adni_cleaned_prc.RData")

data.surv <- df.surv_preds_prc
data.long <- df.long_censored_prc

# Set up

# Build model using all possible longitudinal variables?
all_long_vars <- TRUE

# Selection of candidate longitudinal covariates

# type of TAU, PTAU and ABETA are character, need to handle non-numerical values first. currently excluded, to be added back later.
exclude_vars <- c("VISCODE", "EXAMDATE", "Years.bl", "Month.bl", "M", "age.fup", "Month", "IMAGEUID", "FSVERSION", "DX", "FLDSTRENG", "COLPROT",
                  "TAU", "PTAU", "ABETA")

# Compute distinct values per variables in dataset and remove covariates with < cut-off obs values per subject

## Cut-off for number of distinct observations to be regarded as candidate longitudinal variables
# [!] Note: there is a bias because I have assumed repetitive values are not meaningful, but that's not always true? e.g. DX
# maybe to count NA instead
avg_obs_lower_lim <- 1.5
temp <- data.long %>%
  group_by(RID) %>%
  dplyr::summarise(across(everything(), n_distinct)) %>%
  dplyr::select(-RID) %>%
  dplyr::summarise(across(everything(), mean)) %>%
  pivot_longer(cols = everything(), names_to = "var_name", values_to = "avg_distinct_values") %>%
  filter(avg_distinct_values > avg_obs_lower_lim) %>%
  arrange(desc(avg_distinct_values))

if (all_long_vars){
  # Exclude irrelevant variables
  y.names <- temp$var_name[!(temp$var_name %in% exclude_vars)]
} else {
  # Subset of longitudinal covariates according to previous literature
  y.names <- c("ADAS13", "MMSE", "RAVLT.immediate", "RAVLT.learning", "FAQ")
}
# Set character vector y.names with the names of the response variables which the LMMs have to be fitted to
# table of candidate long covariates and their average distinct values
candidate_long_vars <- temp %>% filter(var_name %in% y.names)

#View(candidate_long_vars)

# pencal parameters
# n.boots set to large number if CBOCP is required
n.boots = 0
n.cores = 1


# Wrapping function to fit PRC model in steps 1 to 3
## Input:
##   long.data: long data
##   surv.data: survival data
##   y.names: vector of candidate longitudinal variable names
## Output:
##   t.cumulated: total runtime
run_prc_steps <- function(
    long.data, surv.data,
    y.names,
    n.boots=0, n.cores=1,
    verbose = FALSE){
  
  t.start <- Sys.time()
  t.cumulated <- t.start
  
  # Step 1 of PRC-LMM: estimate the LMMs
  step1 <- fit_lmms(
    y.names = y.names,
    fixefs = ~ age.fup, # Fixed effects
    ranefs = ~ age.fup | id, # Random effects
    long.data = training.long,
    surv.data = training.surv, 
    t.from.base = Years.bl, # Name of the variable containing time from baseline in long.data
    n.boots = n.boots, 
    n.cores = n.cores)
  
  t.end <- Sys.time()
  if(verbose){
    cat("\nStep 1 completed in", difftime(t.end, t.start, units = "mins"), "min\n")
  }
  
  t.start <- Sys.time()
  
  # Step 2 of PRC-LMM: compute the summaries
  step2 <- summarize_lmms(object = step1, n.cores = n.cores)
  
  t.end <- Sys.time()
  if(verbose){
    cat("\nStep 2 completed in", difftime(t.end, t.start, units = "mins"), "min\n")
  }
  
  t.start <- Sys.time()
  
  # Step 3 of PRC-LMM: fit the penalized Cox models
  step3 <- fit_prclmm(
    object = step2, 
    surv.data = training.surv, 
    baseline.covs = ~ AGE + PTGENDER + PTEDUCAT + status.bl + APOE4, 
    penalty = 'ridge', 
    n.cores = n.cores)
  
  t.end <- Sys.time()
  if (verbose){
    cat("\nStep 3 completed in", difftime(t.end, t.start, units = "mins"), "min\n")
  }
  
  t.cumulated <- difftime(t.end, t.cumulated, units = "mins")
  
  return(t.cumulated)
}

# Scaling covariates
data.long_scaled <- data.long

for (y in y.names){
  # Scaling
  data.long_scaled[, y] <- scale(data.long[, y], center = T, scale = T)
  # Uncomment to print statistics after transformation
  #cat("variable", y, "after scaling:", "\nmean =", mean(data.long_scaled[, y], na.rm = T), ", sd:", sd(data.long_scaled[, y], na.rm = T), "\n")
}

# Stratified train-test split

# Train set
training.surv <- data.surv %>% splitstackshape::stratified(., group = "event", size = 0.8)
training.long <- data.long %>% filter(id %in% training.surv$id)

# Test set (complement to train set)
testing.surv <- data.surv %>% filter(!(id %in% training.surv$id))
testing.long <- data.long %>% filter(id %in% testing.surv$id)

print("Statistics after split")

# Check size
print(paste("Dataset # subjects =", dim(data.surv)[1])) # Dataset with id column
print(paste("Training set # subjects =", dim(training.surv)[1]))
print(paste("Testing set # subjects =", dim(testing.surv)[1]))

# Check class proportion
print("Dataset event proportion:")
prop.table(table(data.surv$event))
print("Training set event proportion:")
prop.table(table(training.surv$event))
print("Testing set event proportion:")
prop.table(table(testing.surv$event))

#training.long %>% group_by(id) %>% dplyr::summarise(count = n())
#testing.long %>% group_by(id) %>% dplyr::summarise(count = n())

# Fit univariate (longitudinal-wise) PRC model for each candidate long covariate to diagnose problem
## Initiate result column
candidate_long_vars$runtime <- NA
## Run pencal
print("Start fitting models...")
candidate_long_vars$runtime <- sapply(y.names, function(y){
  y.names_subset <- as.vector(y)
  cat("\n-----\nFitting long covariate:", y.names_subset, "\n")
  try(run_prc_steps(
    long.data = training.long, 
    surv.data = training.surv, 
    y.names = y.names_subset, 
    n.boots = n.boots, 
    n.cores = n.cores), silent = FALSE)
})
print("All models fitted. See result in `candidate_long_vars`")
# Result table
## Does it run? (Y)-with total runtime in mins; (N)-with error message
candidate_long_vars
#View(candidate_long_vars)
