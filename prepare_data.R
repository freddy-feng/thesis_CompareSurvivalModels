# -----------------------------------------------------------------------------
# Run script to prepare adnimerge data into two formats
# A survival data and a longitudinal data respectively
# Outputs are saved in .RData file
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------
#rm(list = ls())
require(tidyverse)
require(ADNIMERGE) # Install from tar ball separately
require(moments)

source("function_utility_exp.R")
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Data cleaning and preprocessing
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Filter subjects with only single visit i.e. no repeated measurement
ids_one_visit <- adnimerge %>%
  group_by(RID) %>%
  dplyr::summarize(n_visit = n()) %>%
  filter(n_visit == 1) %>%
  select(RID) %>%
  unlist(use.names = FALSE)

paste("[Count] Subjects with single visit only =", length(ids_one_visit))
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Filter subjects without any diagnosis at baseline
ids_DXbl_NA <- adnimerge %>%
  filter(VISCODE == "bl") %>%
  filter(is.na(DX)) %>%
  select(RID) %>%
  unlist(use.names = FALSE)

paste("[Count] Subjects with DX==NA at baseline =", length(ids_DXbl_NA))
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Filter subjects with Dementia diagnosis at baseline
ids_DXbl_Dem <- adnimerge %>%
  filter(VISCODE == "bl") %>%
  filter(DX == "Dementia") %>%
  select(RID) %>%
  unlist(use.names = FALSE)

paste("[Count] Subjects with DX==Dementia at baseline =", length(ids_DXbl_Dem))
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
vars_baseline <- c("AGE", "PTGENDER", "PTEDUCAT", "PTETHCAT", "PTRACCAT", "PTMARRY", "APOE4")

# Filter subjects with NA or unknown value in baseline covariate of interest
ids_base_NA <- lapply(vars_baseline, function(x) {
  ids_NA <- adnimerge %>%
    filter(VISCODE == "bl") %>%
    filter(is.na(get(x)) | get(x) == "Unknown") %>%
    select(RID) %>%
    unlist(use.names = FALSE)
  print(paste("[Count] Subjects with NA or unknown in baseline covariate", 
              x, " =", length(ids_NA)))
  return (ids_NA)
})
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Construct union of subjects to exclude
#ids_exclude <- union(union(ids_one_visit, ids_DXbl_Dem), ids_DXbl_NA)
ids_exclude <- Reduce(union, 
                      c(ids_base_NA, list(ids_one_visit, ids_DXbl_Dem, ids_DXbl_NA)))
paste("[Count] Subjects to be excluded due to single visit or DX==AD/NA at baseline =", 
      length(ids_exclude))
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Subset adnimerge to exclude subjects identified above
df.clean <- adnimerge %>%
  filter(!(RID %in% ids_exclude))
# -----------------------------------------------------------------------------

# 
# # -----------------------------------------------------------------------------
# # Remove all attributes
# This block does not work yet. same for data.long because some variables are factor and will be affected
# df.clean[] <- lapply(df.clean, function(x) {
#   attributes(x) <- NULL
#   return(x)
# })
# # -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
## Transform categorical variables to factors
vars_factorize <- c("PTGENDER", "PTETHCAT", "PTRACCAT", "PTMARRY", "APOE4")

for (x in vars_factorize){
  df.clean[, x] <- as.factor(df.clean[, x])
  print(paste("[Check] Levels in factor", x, 
              "--------------------------------------------------"))
  print(levels(df.clean[, x]))
}
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# Prepare df.surv: encode survival outcome and compute survival time
# where 0=censored, 1=dementia
# Note that we use `Years.bl` as time variable, which is an exact time since baseline

# Event is DX==Dementia from non-dementia at the EARLIEST follow-up time
df.dementia <- df.clean %>%
  filter(DX == "Dementia", 
         VISCODE != "bl") %>% # Dementia after baseline
  group_by(RID) %>%
  dplyr::summarize(
    time = min(Years.bl), # Earliest follow up time in years from baseline
    event = 1)

# DX==non-Dementia at LATEST follow up time
df.censored <- df.clean %>%
  filter(!(RID %in% df.dementia$RID), 
         VISCODE != "bl") %>% # Exclude subjects with event, exclude baseline
  filter(!(is.na(DX))) %>% # Non-NA, DX should be either CN or MCI
  group_by(RID) %>%
  dplyr::summarize(
    time = max(Years.bl), # Latest follow up in years from baseline
    event = 0) # Censored can be CN or MCI, at latest follow up

# Merge event and censored {RID, time, status}
# Basic layout of surv data layout, without predictors
df.surv <- rbind(df.dementia, df.censored)

# Change status to factor
df.surv$status <- df.surv$event %>% 
  factor(labels = c("censored", "dementia"))

# Add baseline diagnosis e.g. CN/MCI as baseline covariate
df.surv <- df.surv %>% 
  left_join(df.clean %>% 
              filter(VISCODE == "bl") %>%
              select(RID, DX), 
            by = "RID") %>%
  mutate(status.bl = DX) %>%
  select(-DX)

# Drop dementia level in status.bl, which is not used due to previous filtering
df.surv$status.bl <- df.surv$status.bl %>% droplevels()

paste("[Count] Subjects in df.surv =", 
      df.surv %>% nrow())
paste("[Count] Subjects developed to dementia during followup =", 
      df.surv %>% filter(event == 1) %>% nrow())
paste("[Count] Subjects not developed to dementia during followup =", 
      df.surv %>% filter(event == 0) %>% nrow())
paste("[Count] overall censoring rate =", 
      1 - df.surv$event %>% mean() %>% round(3))
paste("[Check] No overlap between subjects in dementia group and censored group =", 
      identical(intersect(df.dementia$RID, df.censored$RID), numeric(0))) # TRUE is ok
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Prepare df.surv_preds: add predictors to survival table

# Join cross-sectional data, measurement at baseline to surv data
# Using RID and baseline VISCODE as unique identifier
df.surv_preds <- df.surv %>%
  left_join(df.clean %>% 
              filter(VISCODE == "bl"), 
            by=c("RID"))
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Prepare df.long: extract long data from adnimerge

# Construct df.long based on subjects in surv_preds

# df.long includes all measurement in adnimerge, without truncation of measurement after event
# Note that in survival analysis we only need the measurement before event, which will be prepared later
df.long <- adnimerge %>%
  filter(RID %in% df.surv$RID) 
# Note: a difference between df.long and df.surv is that the former has factorize some baseline covariates

# # -----------------------------------------------------------------------------
# # Remove all attributes
# df.long[] <- lapply(df.long, function(x) {
#   attributes(x) <- NULL
#   return(x)
# })
# # -----------------------------------------------------------------------------


# Add survival data
df.long <- df.long %>% 
  left_join(
    # Add surv time, outcomes, baseline status (CN/MCI/Dementia)
    df.surv_preds %>% 
      select(RID, time, event, status, status.bl), 
    by = "RID"
  )

# Compute and include age.bl (age at baseline) and age.fup (age at follow up)
df.long <- df.long %>%
  mutate(age.bl = AGE,
         age.fup = AGE + Years.bl)

# Include year based version of M, reminder: M contains integer only
df.long <- df.long %>%
  mutate(Y = M / 12)

# Rearrange column orders: {surv data}, {time vars}, {remaining}
df.long <- df.long %>%
  select(RID, time, event, status, status.bl,
         VISCODE, EXAMDATE, Years.bl, Y, Month.bl, M, Month,
         AGE, age.bl, age.fup,
         everything())

# [Manual cleaning] Remove double entry at time variable M
print("[Reminder] Manual removal of RID==6014, VISCODE==m0
      due to double entry at time variable M==0,
      removed row contains NA measurement")

df.long <- df.long %>% 
  filter(!(RID == 6014 & VISCODE == "m0"))

print("[Report] Manually removed observation in subject 6014, VISCODE `m0`. Because it conflicts with the row of VISCODE `bl`.")
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Prepare df.long_censored: remove observations at or after event
# Prepare df.long_censored by filtering
df.long_censored <- df.long %>% 
  filter(time > Years.bl)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Finalize

# Copy `RID` to `id` to comply with `pencal` naming convention
df.surv_preds$id <- df.surv_preds$RID
df.long_censored$id <- df.long_censored$RID

# Reorganize column order
df.surv_preds <- df.surv_preds %>% 
  select(id, time, event, status, status.bl, everything())
df.long_censored <- df.long_censored %>% 
  select(id, everything())

# Arrange in ascending order by id and id, age.fup
df.surv_preds <- df.surv_preds %>% arrange(id)
df.long_censored <- df.long_censored %>% arrange(id, age.fup)
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Transformation
# -----------------------------------------------------------------------------
data.long <- df.long_censored

# Obtain long covariates to transform
# Manual exclusion
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

vars_transform <- names(data.long)[!(names(data.long) %in% c("id", vars_manual_remove, vars_irrelevant))]

res.transformation <- Transform_covariates(
  data.long = data.long, 
  y.names = vars_transform,
  threshold.skew = 0.5,
  threshold.sym = 0.2) # Not used as criteria

Plot_transformation(
  data.long = data.long,
  data.long.transformed = res.transformation$data.long.transformed,
  y.stats.transformed = res.transformation$summary
)

df.long_censored_transformed <- res.transformation$data.long.transformed
transformation_summary <- res.transformation$summary
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Save prepared datasets
save(df.surv_preds, 
     df.long_censored,
     df.long_censored_transformed,
     transformation_summary,
     file = "adni_cleaned.RData")

print("Finished data preparation for adni_cleaned.RData")
# -----------------------------------------------------------------------------


# Consider adding a rm here to remove the intermediate objects...
