# ----------------------------------------------------------------------------------------
# Comparison between two evaluation implementation on pencal model trained on full data
# ----------------------------------------------------------------------------------------
rm(list = ls())
# -----------------------------------------------------------------------------------  
library(tidyverse)
library(pencal)

# Import data
load("adni_cleaned_prc.RData")
data.surv <- df.surv_preds_prc %>% arrange(id) # Arrange by id fix the problem
data.long <- df.long_censored_prc %>% arrange(id, age.fup)

# Stupid but try to reduce unnecessary changes in moving codes
surv.data <- data.surv
long.data <- data.long

# -----------------------------------------------------------------------------------  
# Model specification
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
# Fit pencal model on full data
# -----------------------------------------------------------------------------------  
# Step 1 of PRC-LMM: estimate the LMMs
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
# -----------------------------------------------------------------------------------  
# Step 2 of PRC-LMM: compute the summaries
step2 <- summarize_lmms(object = step1, n.cores = n.cores)
# -----------------------------------------------------------------------------------  
# Step 3 of PRC-LMM: fit the penalized Cox models
step3 <- fit_prclmm(
  object = step2, 
  surv.data = surv.data,
  baseline.covs = as.formula(paste("~", paste(baseline.covs, collapse = "+"))), 
  penalty = penalty,
  n.cores = n.cores)

print("Completed training model on full data")
# ----------------------------------------------------------------------------------- 

# -----------------------------------------------------------------------------------  
# Evaluation 1 on full using the built-in function
# -----------------------------------------------------------------------------------  
eval.pencal <- performance_prc(step2, step3, times = 1:15, n.cores = n.cores)
# ----------------------------------------------------------------------------------- 

# -----------------------------------------------------------------------------------  
# Evaluation 2 using custom functions
# -----------------------------------------------------------------------------------  
T.start <- 0
deltaT <- 1:15

surv.new <- surv.data # dataframe with baseline covariates for the new subjects 
long.new <- long.data # longitudinal data if you want to compute predictions for new subjects

# Obtain predicted random effect
preds <- survpred_prclmm(
  step1 = step1, 
  step2 = step2, 
  step3 = step3,
  times = deltaT, # Prediction window(s) for surv prob prediction
  new.longdata = long.new, # Long data in test set
  new.basecovs = surv.new, # Surv data in test set
  keep.ranef = TRUE # should a data frame with the predicted random effects be included in the output
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

    predict.time <- deltaT[j]
  
  #    print(predict.time)
  
  temp <- survivalROC::survivalROC(
    Stime = surv.new$time, # Event time or censoring time for subjects
    status = surv.new$event, # Indicator of status, 1 if death or event, 0 otherwise
    marker = linpred, # Predictor or marker value
    entry = NULL, # Entry time for the subjects, default is NULL, why 0?
    predict.time = predict.time, # Time point of the ROC curve
    cut.values = NULL, # marker values to use as a cut-off for calculation of sensitivity and specificity
    method = "NNE", 
    span = 0.25 * nrow(surv.new)^(-0.2) # small span yield moderate smoothing
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
  surv.time = surv.new$time, # vector of event times
  surv.event = surv.new$event, # vector of event occurence indicators
  method = "noether" # conservative, noether or name (see paper Pencina et al. for details)
)

c(c.custom = res.c.naive$c.index, c.pencal = eval.pencal$concordance$C.naive)
# -----------------------------------------------------------------------------------
