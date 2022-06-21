# ------------------------------------------------------------------------------------------------
# Return a list of objects containing TP, FP, tdAUC
# ------------------------------------------------------------------------------------------------
require(survivalROC)

Evaluate_tdauc <- function(
    surv.new, 
    linpred, 
    T.start, 
    deltaT) {
  
  res.tp <- vector(mode = "list", length = length(deltaT))
  res.fp <- vector(mode = "list", length = length(deltaT))
  res.tdauc <- vector(mode = "numeric", length = length(deltaT))
  
  # Remove deltaT before landmark time
  deltaT <- deltaT[deltaT > T.start]
  
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
#
# ------------------------------------------------------------------------------------------------
# Summarize tdAUC in different folds from folds.eval and output dataframe
# ------------------------------------------------------------------------------------------------
Summarize.tdAUC <- function(model, path, n_fold, detlaT) {
  # Load results to folds.eval list for corresponding model
  load(path)
  
  # Extract from folds.eval into list of tdAUC
  list.tdauc <- lapply(1:n_fold, function(i) {
    deltaT.eval <- folds.eval[[i]]$perf$deltaT
    folds.eval[[i]]$perf$tdauc[deltaT.eval %in% detlaT]
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
  df.tdauc$model <- model
  
  return(df.tdauc)
}
# ------------------------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------------------------
# Summarize C-index in different folds from folds.eval and output dataframe
# ------------------------------------------------------------------------------------------------
Summarize.c.index <- function(model, path, n_fold) {
  # Load results to folds.eval list for corresponding model
  load(path)
  
  # Extract from folds.eval into list of tdAUC
  vec.c.index <- sapply(1:n_fold, function(i) {
    folds.eval[[i]]$perf$c.index
  })
  
  # Convert list into data frame
  df.c.index <- data.frame(c.index = vec.c.index) # row i = fold i
  df.c.index$model <- model
  
  return(df.c.index)
}
# ------------------------------------------------------------------------------------------------
#
# 
# ------------------------------------------------------------------------------------------------
# Generate a list of plots containing all tdROC for each fold
# ------------------------------------------------------------------------------------------------
Plot_all_tdROC <- function(folds.eval) {
  out <- list()
  deltaT <- folds.eval[[1]]$perf$deltaT
  
  for (i in 1:n_fold) {
    tdroc.list <- lapply(1:length(deltaT), function(j) {
      tp <- folds.eval[[i]]$perf$tp[[j]]
      fp <- folds.eval[[i]]$perf$fp[[j]]
      ggplot(data = data.frame(TP = tp, FP = fp)) +
        geom_point(aes(x = FP, y = TP), size = 0.2) + 
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
        labs(
          title = paste("tdROC for prediction time =", deltaT[j]),
          subtitle = paste("tdAUC =", round(folds.eval[[i]]$perf$tdauc[[j]], 3))
        )
    })
    out[[i]] <- ggpubr::ggarrange(plotlist = tdroc.list, ncol = 5, nrow = 3)
  }
  return (out)
}
# ------------------------------------------------------------------------------------------------
