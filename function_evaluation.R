# ------------------------------------------------------------------------------------------------
# Return a list of objects containing TP, FP, tdAUC
# ------------------------------------------------------------------------------------------------
require(tidyverse)
require(survivalROC)


Evaluate_tdauc <- function(
    surv.new, 
    linpred, 
    T.start, 
    deltaT) {
  
  res.tp <- vector(mode = "list", length = length(deltaT))
  res.fp <- vector(mode = "list", length = length(deltaT))
  res.tdauc <- vector(mode = "numeric", length = length(deltaT))
  
  # This line is no more required after shifting the time scale
  # # Remove deltaT before landmark time
  # deltaT <- deltaT[deltaT > T.start]
  
  for (j in 1:length(deltaT)) {
    
    predict.time <- deltaT[j]
    
#    print(predict.time)
    event.times <- surv.new$time[surv.new$event == 1]
    if (all(!(event.times <= predict.time))) {
      mess <- paste(
        "No event (surv.new$event == 1) is observed between landmark time", T.start, 
        "and prediction time", predict.time)
      warning(mess)
      # Store no result
      res.tdauc[j] <- NA
      res.tp[[j]] <- NA
      res.fp[[j]] <- NA
    } else {
      temp <- survivalROC::survivalROC(
        Stime = surv.new$time, # Event time or censoring time for subjects
        status = surv.new$event, # Indicator of status, 1 if death or event, 0 otherwise
        marker = linpred, # Predictor or marker value
        entry = NULL, # Entry time for the subjects, default is NULL
        predict.time = predict.time, # Time point of the ROC curve
        cut.values = NULL, # marker values to use as a cut-off for calculation of sensitivity and specificity
        method = "NNE", 
        span = 0.25 * nrow(surv.new)^(-0.2) # small span yield moderate smoothing, how to select?
      )
      # Store result
      res.tdauc[j] <- temp$AUC
      res.tp[[j]] <- temp$TP
      res.fp[[j]] <- temp$FP
    }
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
Summarize.tdauc <- function(name, path, deltaT) {
  # Load folds.eval
  load(path)
  
  n_fold <- length(folds.eval)
  
  # Extract from folds.eval into list of tdAUC
  list.tdauc <- lapply(1:n_fold, function(i) {
    # in case deltaT in fold.eval > deltaT
    deltaT.eval <- folds.eval[[i]]$perf$deltaT + folds.eval[[i]]$perf$landmark # shift back the scale
    folds.eval[[i]]$perf$tdauc[deltaT.eval %in% deltaT]
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
      mean = mean(tdAUC, na.rm = TRUE),
      median = median(tdAUC, na.rm = TRUE),
      sd = sd(tdAUC, na.rm = TRUE),
      ci.upper = quantile(tdAUC, 0.975, na.rm = TRUE),
      ci.lower = quantile(tdAUC, 0.025, na.rm = TRUE)
    ) 
  
  df.tdauc$prediction_time <- as.numeric(df.tdauc$prediction_time)
  df.tdauc <- df.tdauc %>% arrange(prediction_time)
  
  # Get model information, identical over all folds
  df.tdauc$model.id <- name # Different seed, different model name
  df.tdauc$model.name <- folds.eval[[1]]$model.info$name
  df.tdauc$method <- folds.eval[[1]]$model.info$hyperparam$method
  df.tdauc$method.full <- folds.eval[[1]]$model.info$hyperparam$method.full
  df.tdauc$scenario <- folds.eval[[1]]$model.info$hyperparam$set_scenario
  df.tdauc$landmark <- folds.eval[[1]]$model.info$hyperparam$landmark
  df.tdauc$n_bl.covariate <- length(folds.eval[[1]]$model.info$covariate$base)
  df.tdauc$bl.covariate <- paste(folds.eval[[1]]$model.info$covariate$base, collapse = "+")
  
  return(df.tdauc)
}
# ------------------------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------------------------
# Summarize C-index in different folds from folds.eval and output dataframe
# ------------------------------------------------------------------------------------------------
Summarize.c.index <- function(model, path) {
  # Load results to folds.eval list for corresponding model
  load(path)
  
  n_fold <- length(folds.eval)
  
  # Extract from folds.eval into list of tdAUC
  vec.c.index <- sapply(1:n_fold, function(i) {
    folds.eval[[i]]$perf$c.index
  })
  
  # Convert list into data frame
  df.c.index <- data.frame(c.index = vec.c.index) # row i = fold i
  
  # Get model information, identical over all folds
  df.c.index$model.name <- folds.eval[[1]]$model.info$name
  df.c.index$method <- folds.eval[[1]]$model.info$hyperparam$method
  df.c.index$method.full <- folds.eval[[1]]$model.info$hyperparam$method.full
  df.c.index$scenario <- folds.eval[[1]]$model.info$hyperparam$set_scenario
  df.c.index$landmark <- folds.eval[[1]]$model.info$hyperparam$landmark
  df.c.index$n_bl.covariate <- length(folds.eval[[1]]$model.info$covariate$base)
  df.c.index$bl.covariate <- paste(folds.eval[[1]]$model.info$covariate$base, collapse = "+")
  
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
# ------------------------------------------------------------------------------------------------
# Summarize Brier score in different folds from folds.eval and output dataframe
# ------------------------------------------------------------------------------------------------
# Note that the length of Brier score may be shorter in some folds if test set does not contain survival times after the floor(T.max)
Summarize.brier <- function(name, path, deltaT) {
  # Load folds.eval
  load(path)
  
  n_fold <- length(folds.eval)
  
  # Extract from folds.eval into list of Brier
  list.brier <- lapply(1:n_fold, function(i) {
    # in case deltaT in fold.eval > deltaT
    deltaT.eval <- folds.eval[[i]]$perf$deltaT + folds.eval[[i]]$perf$landmark # shift back the scale
    brier <- folds.eval[[i]]$perf$brier
    missing <- length(deltaT.eval) - length(brier)
    if (missing > 0) {
         brier <- c(brier, rep(NA, times = missing))
       }
    brier <- brier[deltaT.eval %in% deltaT]
  })
  
  
  # Convert list into data frame
  df.brier <- data.frame(do.call(rbind, list.brier)) # row i = fold i
  colnames(df.brier) <- deltaT
  
  # Compute statistics
  df.brier <- df.brier %>%
    pivot_longer(
      cols = everything(),
      names_to = "prediction_time",
      values_to = "BS"
    ) %>%
    group_by(prediction_time) %>%
    summarise(
      mean = mean(BS, na.rm = TRUE),
      median = median(BS, na.rm = TRUE),
      sd = sd(BS, na.rm = TRUE),
      ci.upper = quantile(BS, 0.975, na.rm = TRUE),
      ci.lower = quantile(BS, 0.025, na.rm = TRUE)
    ) 
  
  
  df.brier$prediction_time <- as.numeric(df.brier$prediction_time)
  df.brier <- df.brier %>% arrange(prediction_time)
  
  # Get model information, identical over all folds
  df.brier$model.id <- name # Different seed, different model name
  df.brier$model.name <- folds.eval[[1]]$model.info$name
  df.brier$method <- folds.eval[[1]]$model.info$hyperparam$method
  df.brier$method.full <- folds.eval[[1]]$model.info$hyperparam$method.full
  df.brier$scenario <- folds.eval[[1]]$model.info$hyperparam$set_scenario
  df.brier$landmark <- folds.eval[[1]]$model.info$hyperparam$landmark
  df.brier$n_bl.covariate <- length(folds.eval[[1]]$model.info$covariate$base)
  df.brier$bl.covariate <- paste(folds.eval[[1]]$model.info$covariate$base, collapse = "+")
  
  return(df.brier)
}
# -----------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------
# Variable importance (use and interpret with care!)
# -----------------------------------------------------------------------------------
# Not ready for use. The coefficients may not be standardized. Need to double check.
Find.var.importance <- function(object, lambda) {
  beta <- predict(object, s = lambda, type = "coefficients")
  if (is.list(beta)) {
    out <- do.call("cbind", lapply(beta, function(x) x[, 1]))
    out <- as.data.frame(out, stringAsFactors = TRUE)
  } else {
    out <- data.frame(coef = beta[, 1])
  }
  out <- abs(out[rownames(out) != "(Intercept)",, drop = FALSE])
  out 
}
# -----------------------------------------------------------------------------------