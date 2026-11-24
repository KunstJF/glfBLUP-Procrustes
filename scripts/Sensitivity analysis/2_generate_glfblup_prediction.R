################################################################################
# generate glfblup predictions using Varimax factor scores, Procrustes rotated factor scores
# or spline characteristic based on Procrustes rotated factor scores
################################################################################

# ----------------------------------------------------
# 1. Load Required Libraries
# ----------------------------------------------------
library(leaps)
library(glfBLUP)
library(magrittr)
library(tibble)
library(rlist)
library(tictoc)
library(tidyverse)

# ----------------------------------------------------
# 2. Source Functions & Load Shared Data
# ----------------------------------------------------
# Dynamically source all custom utility R functions stored in the local 'functions/' folder
for (file.name in list.files("functions")) {
  source(paste0("functions/", file.name))
}

# Load the genomic kinship matrix used across all prediction models
K <- load('./data/processed/kinship.RData')

# ----------------------------------------------------
# 3. Define Experimental Settings Grid
# ----------------------------------------------------
# Create all combinations of Cross-Validation schemes, environments, and plant stages
settings_df <- expand.grid(
  CV = c('CV1', 'CV2'),
  selected_env = c('B5IR', 'B2IR', 'DRIP'),
  plant_stage = c('full', 'heading', 'vegetative'),
  m = 1:3,
  stringsAsFactors = FALSE
)

# ====================================================
# PART A: VARIMAX ROTATED FACTOR SCORES PREDICTION
# ====================================================
for(settings_i in 1:nrow(settings_df)){
  # Extract experimental variables for current configuration
  (CV = settings_df[settings_i,]$CV)
  (env = settings_df[settings_i,]$selected_env)
  (stage = settings_df[settings_i,]$plant_stage)
  (m = settings_df[settings_i,]$m)

  # Check available dataset count for the target environment
  n_datasets <- length(list.files(sprintf('./data/processed/sensitivity_analysis/%s/factors/m_%d', env, m)))

  tic() # Start execution timer
  for(i in 1:n_datasets){
    cat('start run', i, env, CV, stage, m, '\n')

    # Load Varimax factor scores and cross-validation hyper-folds for dataset 'i'
    factor_loadings_file <- sprintf('./data/processed/sensitivity_analysis/%s/factors/m_%d/factors_m_%d_dataset_%d.RData',
                                    env, m, m, i)

    load(factor_loadings_file)
    load(sprintf('./data/processed/hyper_folds/%s/hyper_dataset_%i.RData', env, i))

    # Set temporal threshold index based on environment heading point
    if(env == 'B5IR'){
      heading_i <- 5
    }else{
      heading_i <- 4
    }

    # --- Wrangle Factor Scores ---
    # Filter factor results corresponding to active cross-validation scheme
    cv_loc <- str_detect(names(factor_result), paste0(tolower(CV), '_factor_dfs'))
    factor_scores <- lapply(factor_result[cv_loc], function(setting){bind_rows(setting, .id = 'timepoint')})
    timepoints <- names(factor_result[cv_loc][[1]])

    # Filter timepoints based on target plant stage (heading vs vegetative vs full)
    if(stage == 'heading'){
      timepoints <- timepoints[1:heading_i]
    }else if(stage == 'vegetative'){
      timepoints <- timepoints[1:(heading_i -1)]
    }

    # Reshape factor scores to wide format and combine with response variable (Y)
    fac_df <- factor_scores[[1]] %>%
      filter(timepoint %in% timepoints) |>
      pivot_wider(names_from = timepoint, values_from = contains('F'), names_glue = '{.value}_{timepoint}') |>
      select(G, contains('F'), yield_adjusted_BLUES, Y) %>%
      mutate(Y = yield_adjusted_BLUES, .keep = 'unused')

    # --- Feature Selection via Best Subset Regression ---
    # Compute genotype-level averages across plots for model feature selection
    mean_df <- fac_df %>%
      group_by(G) %>%
      summarise(across(everything(), mean)) %>%
      column_to_rownames('G')

    # Identify optimal subset of factor features minimizing Bayesian Information Criterion (BIC)
    best.fit <- regsubsets(Y~ ., data = na.omit(mean_df), nvmax = ncol(mean_df) - 1)
    sum.fit <- summary(best.fit)
    best.coefs <- coef(best.fit, which.min(sum.fit$bic))

    # --- Fit glfBLUP Model ---
    # Subset kinship matrix to genotypes present in current factor dataset
    k_included <- kinship[unique(fac_df$G), unique(fac_df$G)]

    tryCatch({
      # Fit genomic feature BLUP using features chosen by subset selection
      glfBLUP_result <- custom_glfBLUP(data = as.data.frame(fac_df), selection = names(best.coefs)[-1], K = k_included)

      # Calculate predictive accuracy (Pearson correlation between predictions and target values)
      glfBLUP_result$rho = cor(glfBLUP_result$preds[match(x$pred.target$G, names(glfBLUP_result$preds))],
                              x$pred.target$pred.target)

      # Save output list structure
      list.save(glfBLUP_result,
                file =   sprintf('./data/processed/sensitivity_analysis/%s/glfBLUP/m_%d/%s/%s/glfBLUP_m_%d_dataset_%d.RData',
                                 env, m, CV, stage, m, i))
      rm(glfBLUP_result)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  toc() # Stop execution timer
  }
}








