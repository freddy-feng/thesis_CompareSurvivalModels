# not used scripts, or superseded, or no more required


# ------------------------------------------------------------------------------------------------
# Circumvent data problem in test set created after splitting
# Drawback is bias in performance evaluation
# select subjects with at least one value in each covariate
# ------------------------------------------------------------------------------------------------
exclude_invalid_testing_subjects <- function(testing.long) {
  
  # Count missing data  
  missing_counts <- testing.long %>%
    select(id, all_of(y.names)) %>%
    group_by(id) %>%
    dplyr::summarise(across(everything(), 
                            function(x) all(is.na(x)))) # Within subject, all y is missing
  # Sum missing counts
  all_NA_counts <- missing_counts %>%
    select(-id) %>%
    rowSums()
  
  # Only retain no missing
  ids.valid <- missing_counts$id[all_NA_counts == 0]
  
  print(paste("[Report] removed", 
              length(unique(testing.long$id)) - length(ids.valid), 
              "subjects from original test set of size =", 
              length(unique(testing.long$id))))
  print(paste("[Count] number of subjects in test set with at least one observation in all longitudinal covariates =",
              length(ids.valid)))
  
  return (ids.valid)
}
# ------------------------------------------------------------------------------------------------



# the following script in the fold i test block
# -----------------------------------------------------------------------------------
# pencal - Temp fix for specific problem
# -----------------------------------------------------------------------------------
# Circumvent data problem in test set created after splitting
# by select subjects with at least one value in each covariate
# drawback is smaller test set and larger bias in performance evaluation
# Use modified/filtered test set instead of original test set

# [Detail] Error description and example:
# Problem encountered when using pencal::survpred_prclmm e.g.
# task 2 failed - "Variable ADAS13: all values are NA for at least 1 subject"

  ids.valid <- exclude_invalid_testing_subjects(long.new)

# Update testing set
  long.new <- long.new %>%
    filter(id %in% ids.valid)
  surv.new <- surv.new %>% 
    filter(id %in% ids.valid)