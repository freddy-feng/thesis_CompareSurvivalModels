
pred_ranefs_scaled <- as.matrix(do.call(cbind,lapply(1:ncol(pred_ranefs), function(x) {
  scale(pred_ranefs[, x])
}))) 

# mean of predicted random effects
pred_ranefs_scaled %>% colSums()

# sd of predicted random effects
sapply(1:ncol(pred_ranefs_scaled), function(x) {
  pred_ranefs_scaled[, x] %>% sd()
})


X.scaled <- as.matrix(cbind(X0, pred_ranefs_scaled))

# Drop intercept
X.scaled <- X.scaled[, -1]
# -----------------------------------------------------------------------------------
# Compute the linear predictor
# -----------------------------------------------------------------------------------
# s - Value(s) of the penalty parameter lambdaat which predictions are required.
# Default is the value s="lambda.1se" stored on the CV object. Alternatively
# s="lambda.min" can be used.

linpred.scaled <- predict(
  object = step3$pcox.orig, # Fitted "cv.glmnet" or "cv.relaxed" object
  newx = X.scaled, # Matrix of new values for x at which predictions are to be made. Must be a matrix
  s = "lambda.min",
  type = "link" # Type "link" (default) returns x^T \beta
)

res.c.naive.scaled <- survcomp::concordance.index(
  x = linpred.scaled, # vector of risk predictions
  surv.time = surv.new$time, # vector of event times
  surv.event = surv.new$event, # vector of event occurence indicators
  method = "noether" # conservative, noether or name (see paper Pencina et al. for details)
)
res.c.naive.scaled$c.index




# mean of predicted random effects
pred_ranefs %>% colSums()

# sd of predicted random effects
sapply(1:ncol(pred_ranefs), function(x) {
  pred_ranefs[, x] %>% sd()
})
# -----------------------------------------------------------------------------------

