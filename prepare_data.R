# -----------------------------------------------------------------------------
# Run script to prepare adnimerge data into two formats
# A survival data and a longitudinal data respectively
# Outputs are saved in .RData file
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------
#rm(list = ls())
library(tidyverse)

library(moments)
library(consort)

library(ADNIMERGE) # Install from tar ball separately

source("function_utility.R")
# -----------------------------------------------------------------------------

# Get all ids for consort diagram
consort.id.raw <- sort(unique(adnimerge$RID))
n_subjects <- length(consort.id.raw)
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

is_single_visit <- consort.id.raw %in% ids_one_visit
consort.id.singlevisit <- rep(NA, n_subjects) # Initialize
consort.id.singlevisit[is_single_visit] <- "Single visit only" # Reason for exclusion
consort.id.interim1 <- consort.id.raw # Initialize
consort.id.interim1[is_single_visit] <- NA # To next step

#paste("[Count] Subjects removed from previous =", sum(!is.na(consort.id.singlevisit)))
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Filter subjects without any diagnosis at baseline
ids_DXbl_NA <- adnimerge %>%
  filter(VISCODE == "bl") %>%
  filter(is.na(DX)) %>%
  select(RID) %>%
  unlist(use.names = FALSE)

paste("[Count] Subjects with DX==NA at baseline =", length(ids_DXbl_NA))

is_DXbl_NA <- consort.id.interim1 %in% ids_DXbl_NA
consort.id.DXbl_NA <- rep(NA, n_subjects) # Initialize
consort.id.DXbl_NA[is_DXbl_NA] <- "No baseline diagnosis" # Reason for exclusion
consort.id.interim2 <- consort.id.interim1
consort.id.interim2[is_DXbl_NA] <- NA

#paste("[Count] Subjects removed from previous =", sum(!is.na(consort.id.DXbl_NA)))

# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Filter subjects with Dementia diagnosis at baseline
ids_DXbl_Dem <- adnimerge %>%
  filter(VISCODE == "bl") %>%
  filter(DX == "Dementia") %>%
  select(RID) %>%
  unlist(use.names = FALSE)

paste("[Count] Subjects with DX==Dementia at baseline =", length(ids_DXbl_Dem))

is_DXbl_Dem <- consort.id.interim2 %in% ids_DXbl_Dem
consort.id.DXbl_Dem <- rep(NA, n_subjects) # Initialize
consort.id.DXbl_Dem[is_DXbl_Dem] <- "Diagnosed as dementia at baseline" # Reason for exclusion
consort.id.interim3 <- consort.id.interim2
consort.id.interim3[is_DXbl_Dem] <- NA

#paste("[Count] Subjects removed from previous =", sum(!is.na(consort.id.DXbl_Dem)))

# -----------------------------------------------------------------------------
ids_DX_allNA <- adnimerge %>%
  filter(VISCODE != "bl") %>%
  group_by(RID) %>%
  dplyr::summarise(
    all_DX_NA = all(is.na(DX))
  ) %>% 
  filter(all_DX_NA == TRUE) %>%
  select(RID) %>%
  unlist(use.names = FALSE)

paste("[Count] Subjects without any DX after baseline", length(ids_DX_allNA))

is_DX_allNA <- consort.id.interim3 %in% ids_DX_allNA
consort.id.DX_allNA <- rep(NA, n_subjects) # Initialize
consort.id.DX_allNA[is_DX_allNA] <- "No diagnosis after baseline" # Reason for exclusion
consort.id.interim4 <- consort.id.interim3
consort.id.interim4[is_DX_allNA] <- NA

#paste("[Count] Subjects removed from previous =", sum(!is.na(consort.id.DX_allNA)))

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

is_no_basecov <- consort.id.interim4 %in% Reduce(union, ids_base_NA)
consort.id.no_basecov <- rep(NA, n_subjects) # Initialize
consort.id.no_basecov[is_no_basecov] <- "NA(s) in baseline covariate(s)" # Reason for exclusion
consort.id.interim5 <- consort.id.interim4
consort.id.interim5[is_no_basecov] <- NA

#paste("[Count] Subjects removed from previous =", sum(!is.na(consort.id.no_basecov)))
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Construct union of subjects to exclude
#ids_exclude <- union(union(ids_one_visit, ids_DXbl_Dem), ids_DXbl_NA)
ids_exclude <- Reduce(union, 
                      c(ids_base_NA, list(ids_one_visit, ids_DXbl_Dem, ids_DXbl_NA, ids_DX_allNA)))
paste("[Count] Subjects to be excluded due to single visit or DX==AD/NA at baseline or no DX after baseline =", 
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
  filter(DX == "Dementia", VISCODE != "bl") %>% # Redundant step to filter Dementia after baseline, already done in previous step
  group_by(RID) %>%
  dplyr::summarize(
    time = min(Years.bl), # Earliest follow up time in years from baseline
    event = 1)

#paste("[Count] Subjects removed from previous =", sum(!is.na(consort.id.no_basecov)))

# DX==non-Dementia at LATEST follow up time

# Temp part to check the filter DX==NA
df.censored <- df.clean %>%
  filter(!(RID %in% df.dementia$RID), VISCODE != "bl") %>% # Censored outcome after baseline
#  filter(!(is.na(DX))) %>% # Non-NA, DX should be either CN or MCI
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
# Prepare df.long_censored: remove observations at or after outcome
# Should we keep observations at outcome? Make sure information or condition after event is not included in long data
# In other words, to ensure the long data represent the same state
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
  summary = res.transformation$summary
)

df.long_censored_transformed <- res.transformation$data.long.transformed
transformation_summary <- res.transformation$summary
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Consort diagram to describe exclusion of subjects
df.consort <- data.frame(
  consort.id.raw,
  consort.id.singlevisit, 
  consort.id.interim1,
  consort.id.DXbl_NA,
  consort.id.interim2,
  consort.id.DXbl_Dem,
  consort.id.interim3,
  consort.id.DX_allNA,
  consort.id.interim4,
  consort.id.no_basecov,
  consort.id.interim5)

# Unite removal columns
df.consort <- df.consort %>% 
  unite(
    "consort.id.singlevisit", 
    "consort.id.DXbl_NA", 
    "consort.id.DXbl_Dem",
    "consort.id.DX_allNA",
    "consort.id.no_basecov", 
    sep = "", remove = TRUE, 
    na.rm = TRUE) %>%
  rename(consort.id.excl = "consort.id.singlevisit")

df.consort$consort.id.excl[df.consort$consort.id.excl == ""] <- NA

g.consort <- consort_plot(data = df.consort,
             orders = c(consort.id.raw = "ADNI participants",
                        consort.id.excl = "Excluded",
                        consort.id.interim4 = "Subjects for modelling"),
             side_box = c("consort.id.excl"),
             cex = 0.7)

ggsave("consort_diagram.png", plot = g.consort, device = "png",
       path = "./data_cleaned/",
       width = 6, height = 2)


# Note that the subjects can be removed due to one or more reasons, so the diagram depends on the order of exclusion.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Save prepared datasets
save(df.surv_preds, 
     df.long_censored,
     df.long_censored_transformed,
     transformation_summary,
     file = "./data_cleaned/adni_cleaned.RData")

print("Finished data preparation for adni_cleaned.RData")
# -----------------------------------------------------------------------------


# Consider adding a rm here to remove the intermediate objects...
