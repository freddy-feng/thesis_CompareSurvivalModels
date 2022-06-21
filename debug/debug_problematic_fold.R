# ----------------------------------------------------------------------------------------
# Comparison between two evaluation implementation on pencal model trained on full data
# ----------------------------------------------------------------------------------------
rm(list = ls())
# -----------------------------------------------------------------------------------  
library(tidyverse)
library(pencal)

# -----------------------------------------------------------------------------------  
# Utility functions
# -----------------------------------------------------------------------------------  
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
# -----------------------------------------------------------------------------------  
# Utility functions
# -----------------------------------------------------------------------------------  
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
# -----------------------------------------------------------------------------------  
# Utility functions
# -----------------------------------------------------------------------------------  
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

# -----------------------------------------------------------------------------------  
# Model specification
# -----------------------------------------------------------------------------------  
baseline.covs <- c("AGE", "PTGENDER", "PTEDUCAT", "status.bl", "APOE4")
#y.names <- folds[[1]]$candidate.long.covs # equivalent
y.names <- c("ADAS11", "ADAS13", "ADASQ4", "CDRSB",
             "Entorhinal", "FAQ", "Fusiform", "Hippocampus",          
             "ICV", "LDELTOTAL", "MidTemp", "MMSE", 
             "mPACCdigit", "mPACCtrailsB", "RAVLT.forgetting",
             "RAVLT.immediate", "RAVLT.learning", 
             "RAVLT.perc.forgetting", "TRABSCOR", "Ventricles",
             "WholeBrain")


penalty <- "ridge"

max.ymissing <- 0.2 # Default 
n.boots <- 0
n.cores <- parallel::detectCores()



# -----------------------------------------------------------------------------------  
# Import data
# -----------------------------------------------------------------------------------  
# Get full data
load("adni_cleaned_prc.RData")
data.surv <- df.surv_preds_prc %>% arrange(id) # Arrange by id to use survpred_prclmm properly
data.long <- df.long_censored_prc %>% arrange(id, age.fup) # Arrange by id and time-variable to use survpred_prclmm properly


# Get subjects in problematic fold
load("debug_problematic_fold.RData") # Retrieve subject ids in the problematic fold
print(ids.test.problematic)
ids.train <- data.surv$id[!(data.surv$id %in% ids.test.problematic)] # Get subject ids for training set


# Derive scaling parameters mu, sigma from training set
scaling_table <- Compute_scaling_table(
  data = data.long %>% filter(id %in% ids.train), 
  vars_scale = vars_to_scale)

# -----------------------------------------------------------------------------------  
# Split (and scale) data
# -----------------------------------------------------------------------------------  
# Split and scale data WITH SCALING
tmp.scaled <- Get_train_test_data(
  data.surv = data.surv, # Full surv data 
  data.long = data.long, # Full long data
  ids.test = ids.test.problematic, # Subject ids for test set
  is_scaled = "scaled", # Set to "scaled" to scale numeric covariates
  scaling_table = scaling_table
)

training.surv.scaled <- tmp.scaled$training.surv
training.long.scaled <- tmp.scaled$training.long

testing.surv.scaled <- tmp.scaled$testing.surv
testing.long.scaled <- tmp.scaled$testing.long


# Split and scale data NO SCALING
tmp.notscaled <- Get_train_test_data(
  data.surv = data.surv, # Full surv data 
  data.long = data.long, # Full long data
  ids.test = ids.test.problematic, # Subject ids for test set
  is_scaled = "notscaled", # Set to "scaled" to scale numeric covariates
  scaling_table = scaling_table
)

training.surv.notscaled <- tmp.notscaled$training.surv
training.long.notscaled <- tmp.notscaled$training.long

testing.surv.notscaled <- tmp.notscaled$testing.surv
testing.long.notscaled <- tmp.notscaled$testing.long



# -----------------------------------------------------------------------------------  
# Optional: uncomment to check scaled vs not scaled for all covariates to scale
# -----------------------------------------------------------------------------------  

# for (x in vars_to_scale) {
#   mu.notscaled <- mean(training.long.notscaled[, x], na.rm = TRUE) %>% round(3)
#   sd.notscaled <- sd(training.long.notscaled[, x], na.rm = TRUE) %>% round(3)
#   mu.scaled <- mean(training.long.scaled[, x], na.rm = TRUE) %>% round(3) 
#   sd.scaled <- sd(training.long.scaled[, x], na.rm = TRUE) %>% round(3) 
#   plot(training.long.scaled[, x], training.long.notscaled[, x])
#   title(paste(x, 
#               "\nmean", mu.notscaled,"sd", sd.notscaled,
#               "mean", mu.scaled, "sd", sd.scaled))
# }


# -----------------------------------------------------------------------------------  
# Choose training and testing set here
# -----------------------------------------------------------------------------------  
# not scaled data
#training.surv <- training.surv.notscaled
#training.long <- training.long.notscaled

# scaled data
training.surv <- training.surv.scaled
training.long <- training.long.scaled

# naive evaluation
#testing.surv <- training.surv
#testing.long <- training.long

# evaluate using test set
#testing.surv <- testing.surv.notscaled
#testing.long <- testing.long.notscaled
testing.surv <- testing.surv.scaled
testing.long <- testing.long.scaled

training.surv %>% nrow()
testing.surv %>% nrow()
union(training.surv$id, testing.surv$id) %>% length()
intersect(training.surv$id, testing.surv$id) %>% length

#rm(list = c("preds", "linpred"))

# -----------------------------------------------------------------------------------  
# Fit pencal model on full data
# -----------------------------------------------------------------------------------  
# Step 1 of PRC-LMM: estimate the LMMs
step1 <- fit_lmms(
  y.names = y.names,
  fixefs = ~ age.fup, # Fixed effects
  ranefs = ~ age.fup | id, # Random effects
  long.data = training.long,
  surv.data = training.surv,
  t.from.base = Years.bl, # Name of the variable containing time from baseline in long.data
  max.ymissing = max.ymissing,
  n.boots = n.boots, 
  n.cores = n.cores)
# -----------------------------------------------------------------------------------  
# Step 2 of PRC-LMM: compute the summaries
step2 <- summarize_lmms(object = step1, n.cores = n.cores)
# -----------------------------------------------------------------------------------  
# Step 3 of PRC-LMM: fit the penalized Cox models
step3 <- fit_prclmm(
  object = step2, 
  surv.data = training.surv,
  baseline.covs = as.formula(paste("~", paste(baseline.covs, collapse = "+"))), 
  penalty = penalty,
  n.cores = n.cores)

print("Completed training model on full data")
# ----------------------------------------------------------------------------------- 

# -----------------------------------------------------------------------------------  
# Evaluation set up
# ----------------------------------------------------------------------------------- 
T.start <- 0
deltaT <- 1:15

# -----------------------------------------------------------------------------------  
# Evaluation 1 on full using the built-in function
# -----------------------------------------------------------------------------------  
eval.pencal <- performance_prc(step2, step3, times = deltaT, n.cores = n.cores)
# ----------------------------------------------------------------------------------- 

# -----------------------------------------------------------------------------------  
# Evaluation 2 using custom functions
# -----------------------------------------------------------------------------------  

# Obtain predicted random effect
preds <- survpred_prclmm(
  step1 = step1, 
  step2 = step2, 
  step3 = step3,
  times = deltaT, # Prediction window(s) for surv prob prediction
  new.longdata = testing.long, # Long data in test set
  new.basecovs = testing.surv, # Surv data in test set
  keep.ranef = TRUE # should a data frame with the predicted random effects be included in the output
)
# -----------------------------------------------------------------------------------  
# Obtain new X for predict on pcox.orig
# With baseline covariate, may need to expand into without baseline covariate case later
# -----------------------------------------------------------------------------------
X0 <- model.matrix(as.formula(step3$call$baseline.covs),
                   data = testing.surv)
# Join the baseline covariates and predicted random effect summary
pred_ranefs <- preds$predicted_ranefs
X.orig <- as.matrix(cbind(X0, as.matrix(pred_ranefs)))

# =======================================
# pred_ranefs given full data = step2$ranef.orig
# squared residual is small
#sum((pred_ranefs - step2$ranef.orig) ^ 2) # Uncomment to check
# =======================================

# Drop intercept
contains.int <- "(Intercept)" %in% colnames(X.orig)
if (contains.int) {
  X.orig <- X.orig[, -1]
}
# -----------------------------------------------------------------------------------
# Compute the linear predictor
# -----------------------------------------------------------------------------------
# s - Value(s) of the penalty parameter lambdaat which predictions are required.
# Default is the value s="lambda.1se" stored on the CV object. Alternatively
# s="lambda.min" can be used.

linpred <- predict(
  object = step3$pcox.orig, # Fitted "cv.glmnet" or "cv.relaxed" object
  newx = X.orig, # Matrix of new values for x at which predictions are to be made. Must be a matrix
  s = "lambda.min", # Is it the choice of s? there is difference in using lambda.1se but not close to pencal perf level
  type = "link" # Type "link" (default) returns x^T \beta
)
# -----------------------------------------------------------------------------------
# Compute tdROC and tdAUC
# -----------------------------------------------------------------------------------

# Exploding Evaluate_tdauc()
#res.tdauc <- Evaluate_tdauc(surv.new, linpred, T.start, deltaT)

# Initialize
res.tp <- vector(mode = "list", length = length(deltaT))
res.fp <- vector(mode = "list", length = length(deltaT))
res.tdauc <- vector(mode = "numeric", length = length(deltaT))

for (j in 1:length(deltaT)) {

  predict.time <-  deltaT[j]
  
  #    print(predict.time)
  
  temp <- survivalROC::survivalROC(
    Stime = testing.surv$time, # Event time or censoring time for subjects
    status = testing.surv$event, # Indicator of status, 1 if death or event, 0 otherwise
    marker = linpred, # Predictor or marker value
    entry = NULL, # Entry time for the subjects, default is NULL, why 0?
    predict.time = predict.time, # Time point of the ROC curve
    cut.values = NULL, # marker values to use as a cut-off for calculation of sensitivity and specificity
    method = "NNE", 
    span = 0.25 * nrow(testing.surv)^(-0.2) # small span yield moderate smoothing
  )
  
  res.tdauc[j] <- temp$AUC
  res.tp[[j]] <- temp$TP
  res.fp[[j]] <- temp$FP
}


# -----------------------------------------------------------------------------------
# Compare tdAUC
# -----------------------------------------------------------------------------------
df.compare.tdauc <- data.frame(eval.pencal$tdAUC) %>% dplyr::select(pred.time, tdAUC.naive)
df.compare.tdauc$custom.tdauc <- res.tdauc

df.compare.tdauc <- df.compare.tdauc %>% dplyr::mutate(r = tdAUC.naive / custom.tdauc)
#View(df.compare.tdauc)

# -----------------------------------------------------------------------------------
# Compare C-index
# -----------------------------------------------------------------------------------
res.c.naive <- survcomp::concordance.index(
  x = linpred, # vector of risk predictions
  surv.time = testing.surv$time, # vector of event times
  surv.event = testing.surv$event, # vector of event occurence indicators
  method = "noether" # conservative, noether or name (see paper Pencina et al. for details)
)

c(c.custom = res.c.naive$c.index, c.pencal = eval.pencal$concordance$C.naive)
# -----------------------------------------------------------------------------------

df.compare.tdauc




colnames(step2$ranef.orig)
