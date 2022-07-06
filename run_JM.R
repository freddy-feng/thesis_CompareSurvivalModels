on.alice <- T

# Set folder path on ALICE
folder.alice <- "/data1/s2887592/20220706_JM"

if (on.alice) setwd(folder.alice)
# ------------------------------------------------------------

library(tidyverse)
library(joineRML)

select <- dplyr::select

load("./data_cleaned/adni_cleaned.RData")
# ------------------------------------------------------------
# prepare data

# reduce sample size for debug
n <- 1000
set.seed(721)
ids.train <- sample(df.surv_preds$id, size = n)

training.surv <- df.surv_preds %>%
  filter(id %in% ids.train)
training.long <- df.long_censored %>% # Use transformed
  filter(id %in% ids.train)

ids.all <- training.surv$id

# ------------------------------------------------------------
training.surv <- training.surv %>%
  mutate(age.fup = time + AGE) %>%
  select(id, time, age.fup, everything())

# variables to keep in the df feed to mjoint
baseline.covs <- c("AGE", "PTGENDER", "PTEDUCAT", "status.bl", "APOE4")
keep_ids <- c("id", "time", "event", "age.fup", "Years.bl", baseline.covs)
long.covs <- c("ADAS13", "MMSE")


# here i need to filter out the NAs for each long covariate repectively
list.df.no_NA <- lapply(long.covs, function(y) {
  
  # remove na
  df.tmp <- training.long[!is.na(training.long[, y]), ]
  
  ids.tmp <- unique(df.tmp$id)
  
  # Check if mean imputation at baseline is needed when any subject does not have a measurement for y
  # when removing NAs, some subjects may be removed because they don't have any long measurements at all
  print(length(ids.tmp) - length(ids.all)) # Check, negative => missing subjects => imputation needed
  if (length(ids.tmp) < length(ids.all)) {
    y.mean <- mean(training.long[training.long$VISCODE == "bl", y], na.rm = TRUE)
    # Manually modify
    ids.impute <- ids.all[!ids.all %in% ids.tmp]
    
    df.impute <- training.long %>% # Impute at baseline and maintain same data format
      filter(id %in% ids.impute) %>%
      filter(VISCODE == "bl")
    
    #    df.impute %>% 
    #      select(!!rlang::enquo(y))
    
    #    this cannot copy label and class attributes
    #    df.impute <- df.impute %>%
    #      mutate(across(.cols = all_of(y), .fns = function(x) y.mean))
    
    df.impute[, y] <- y.mean # Impute
    
    # Because ADNI dataframe is labelled. It's probably better to remove all labels at first?
    #    df.impute <- sjlabelled::copy_labels(df_new = df.impute, df_origin = training.long)
    
    df.tmp <- rbind(df.tmp, df.impute) %>%
      arrange(id, Years.bl)
    
    print(paste("Subjects after imputation =", length(unique(df.tmp$id)))) # After imputation
  }
  
  # drop irrelevant data
  df.tmp <- df.tmp[, c(keep_ids, y)]
  
  return(df.tmp)
})

# save memory
rm(df.long_censored, df.long_censored_transformed, df.surv_preds)

training.surv <- training.surv %>%
  select(all_of(keep_ids))
# ------------------------------------------------------------
# training

t.start <- Sys.time()

fit.joineRML = mjoint(
  # a list of formulae for the fixed effects component of each longitudinal outcome
  # The left hand-hand side defines the response, and the right-hand side specifies the fixed effect terms
  formLongFixed = list(
    "ADAS13" = ADAS13 ~ age.fup, 
    "MMSE" = MMSE ~ age.fup),
  # a list of one-sided formulae specifying the model for the random effects effects of each longitudinal outcome.
  formLongRandom = list( 
    "ADAS13" = ~ age.fup | id, 
    "MMSE" = ~ age.fup | id),
  formSurv = Surv(age.fup, event) ~ AGE + PTGENDER + PTEDUCAT + status.bl + APOE4,
  #  list of data.frame objects for each longitudinal outcome in which to interpret the variables named in the formLongFixed and formLongRandom.
  #  If the multiple longitudinal outcomes are measured at the same time points for each patient, then a data.frame object can be given instead of a list
  data = list.df.no_NA,
  survData = training.surv,
  timeVar = "age.fup") # Did not use age.fup as time variable in LMM, due to error

t.end <- Sys.time()

training.time <- difftime(t.end, t.start, units = "mins")

cat("Completed in", round(training.time, 3), "min")
# ------------------------------------------------------------

save(list = c("fit.joineRML", "training.time"), file = "./model_JM.RData")
