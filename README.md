# Introduction
This folder contains the scripts for reproducing the results in my master theis "An empirical comparison of methods for the prediction of survival from many longitudinal covariates".  
The aim of the thesis is to empirically compare statistical methods for predicting the survival outcome using a large number of longitudinal covariates.  
Four methods were compared:
- pCox-baseline - penalized Cox model using baseline measurements only
- pCox-landmarking - penalized Cox model using last observations carried forward
- PRC-LMM - penalized regression calibration (Signorelli et al., 2019, https://onlinelibrary.wiley.com/doi/full/10.1002/sim.9178)
- MFPCCox - multivariate functional principal component analysis and Cox regression (Li et al., 2019, https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.8334)  
There are various types of scripts in this repo, indicated by their prefixes:
- `function_`: wrapper functions for model training, evaluation; utility functions for data preparation, processing result.
- `prepare_data`: preparation of raw data into cleaned data used for model fitting
- `run_`: running RCV to estimate the predictive performance of each method; result comparisons
- `IDA_`: exploratory data analysis, to produce some figures and tables presented in the thesis

# Acknowledgement
The scripts for MFPCCox was based on the implementation of MFPCCox https://github.com/kan-li/MFPCCox and adapted for application to ADNI data.

# Data availability
Data can be obtained from the Alzheimer’s Disease Neuroimaging Initiative (ADNI) database (https://adni.loni.usc.edu/). 

# How to run

## Step 1 - Set up
- Install the R package `ADNIMERGE` (see Data availability above).
- Check the folder structure:  
Root  
|_ data_cleaned  
|_ figure  
|_ model  
|_ output  


## Step 2 - Prepare data 
- Run `prepare_data.R` to prepare two cleaned dataframes used for model fitting and testing:
    - a surv data `df.surv_preds`: containg survival time, status, and baseline covariates (time-independent)
    - a long data `df.long_censored_transformed` (`df.long_censored` if transformation is not needed)
- The following steps are executed in the script:
    - Data screening
    - Remove data after survival time
    - Reformat `adnimerge` dataframe to `surv` and `long` format respectively
    - Data transformation
- The resulting files are saved to `/data_cleaned/adni_cleaned.RData` for reuse.

## Step 3 - Train and evaluate models
- The following scripts will train and evaluate each method using repeated 10-fold CV, with stratified sampling.
- Run the corresponding script for each method: 
    - pCox-baseline: `run_repCV_glmnet.R`, set `method <- "pCox-bl"` at the start of script
    - pCox-landmarking: `run_repCV_glmnet.R`, set `method <- "pCox-lm"` at the start of script
    - PRC-LMM: `run_repCV_pencal.R`
    - MFPCCox: `run_repCV_MFPCCox.R`
- Trained models are contained within a list `folds` of length 10 (for 10-fold CV).
- Separately, evaluation results for corresponding model are contained within a list `folds.eval`.
- For each CV in RCV, `folds.eval` will be saved to the a `.RData` in subfolder `./output/`
    - file naming is identifiable as "eval_ + model.name + _seed.RData", which model.name is an automatically generated string combining {method, hyperparameters, landmark time, data setting}
    - e.g. `eval_pCox-bl-min-ridge_lm2_scenario2_b5_transformed_scaled_seed721.RData`
- Trained model is not saved by default. Uncomment the code block `# Save folds for training and future checking` to save the trained models in similar fashion.
- Check settings before running:
    - `n_RCV`: number of repetitions for RCV
    - `T_LMs`: landmark times
    - Code block in `# Set model hyperparam` to change method-specific model hyperparameters
        - `is_transformed`: "transformed" (Transform covariates to reduce skewness)
        - `is_scaled`: "scaled" (scale covariates before model fitting)
    - Code block in `# Set data param` to change the data setting i.e. candidate covariates. Options:
        - `set_scenario` accepts the following values: 
            - "scenario2" (default) uses all logitudinal covariates that satistied `missing_proportion_limit` (default set at 0.1)
            - "scenario1": `missing_proportion_limit` set to 0. Every subject will have at least one observations in each longitudinal covariate selected
            - "scenario0": do not use any longitudinal covariates

## Step 4 - Post-processing results
- Open `run_comparison.Rmd`.
- Evaluation files i.e. the `folds.eavl` saved in subfolder `./output/` will be loaded, aggregated, and compared.
- The results are used to reproduce figures in the thesis.

