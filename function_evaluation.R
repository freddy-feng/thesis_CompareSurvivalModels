# ------------------------------------------------------------------------------------------------
# Return a list of objects containing TP, FP, tdAUC
# ------------------------------------------------------------------------------------------------
require(survivalROC)

Evaluate_tdauc <- function(surv.new, linpred.orig, deltaT) {
  res.tp <- vector(mode = "list", length = length(deltaT))
  res.fp <- vector(mode = "list", length = length(deltaT))
  res.tdauc <- vector(mode = "numeric", length = length(deltaT))
  
  for (j in 1:length(deltaT)) {
    
    predict.time <- deltaT[j]
    
    temp <- survivalROC::survivalROC(
      Stime = surv.new$time, # Event time or censoring time for subjects
      status = surv.new$event, # Indicator of status, 1 if death or event, 0 otherwise
      marker = linpred.orig, # Predictor or marker value
      entry = NULL, # Entry time for the subjects, default is NULL, why 0?
      predict.time = predict.time, # Time point of the ROC curve
      cut.values = NULL, # marker values to use as a cut-off for calculation of sensitivity and specificity
      method = "NNE", 
      span = 0.25 * nrow(surv.new)^(-0.2) # small span yield moderate smoothing, how to select?
    )
    
    res.tdauc[j] <- temp$AUC
    res.tp[[j]] <- temp$TP
    res.fp[[j]] <- temp$FP
  }
  
  return(list(
    tp = res.tp,
    fp = res.fp,
    tdauc = res.tdauc
  ))
}






