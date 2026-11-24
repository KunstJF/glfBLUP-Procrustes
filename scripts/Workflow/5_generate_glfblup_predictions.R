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

  # Check available dataset count for the target environment
  n_datasets <- length(list.files(
    sprintf('./data/processed/factors/%s/timepoint_factors/varimax/', env)))

  tic() # Start execution timer
  for(i in 1:100){
    cat('start run', i, env, CV, stage, '\n')

    # Load Varimax factor scores and cross-validation hyper-folds for dataset 'i'
    load(sprintf('./data/processed/factors/%s/timepoint_factors/varimax/factors_dataset_%d.RData', env, i))
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
      filter(timepoint %in% timepoints) %>%
      pivot_wider(names_from = timepoint, values_from = F1:F2) %>%
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
                file = paste0('./data/processed/glfBLUP_timepoints/varimax/',env, '/', CV,
                              '/', stage, '/glfBLUP_', i, '.RData'))
      rm(glfBLUP_result)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  toc() # Stop execution timer
}

# ====================================================
# PART B: PROCRUSTES ROTATED FACTOR SCORES PREDICTION
# ====================================================
for(settings_i in 1:nrow(settings_df)){
  (CV = settings_df[settings_i,]$CV)
  (env = settings_df[settings_i,]$selected_env)
  (stage = settings_df[settings_i,]$plant_stage)

  n_datasets <- length(list.files(
    sprintf('./data/processed/factors/%s/timepoint_factors/varimax/', env)))

  tic()
  for(i in 1:100){
    cat('start run', i, env, CV, stage, '\n')

    # Load Procrustes aligned factor scores and hyper-folds
    load(sprintf('./data/processed/factors/%s/timepoint_factors/procrustes/procrustes_%d.RData', env, i))
    load(sprintf('./data/processed/hyper_folds/%s/hyper_dataset_%i.RData', env, i))

    if(env == 'B5IR'){
      heading_i <- 5
    }else{
      heading_i <- 4
    }

    # --- Wrangle Procrustes Factors ---
    factor_scores <- lapply(procrustes_factor_results[[CV]], function(setting){setting$factor_scores})
    factor_scores_df <- bind_rows(factor_scores, .id = 'timepoint')
    timepoints <- names(factor_scores)

    # Filter timepoints based on target plant stage
    if(stage == 'heading'){
      timepoints <- timepoints[1:heading_i]
    }else if(stage == 'vegetative'){
      timepoints <- timepoints[1:(heading_i -1)]
    }

    # Format factor scores and join with design metadata/yield targets
    fac_df <- factor_scores_df %>%
      filter(timepoint %in% timepoints) %>%
      pivot_wider(names_from = timepoint, values_from = F1:F2)

    full_df <- x$data %>%
      left_join(fac_df, by = c('G', 'row', 'column', 'plot')) %>%
      select(G, contains('F'), yield_adjusted_BLUES.x) %>%
      mutate(Y = yield_adjusted_BLUES.x, .keep = 'unused')

    # --- Feature Selection via Best Subset Regression ---
    mean_df <- full_df %>%
      group_by(G) %>%
      summarise(across(everything(), mean)) %>%
      column_to_rownames('G')

    best.fit <- regsubsets(Y~ ., data = na.omit(mean_df), nvmax = ncol(mean_df) - 1)
    sum.fit <- summary(best.fit)
    best.coefs <- coef(best.fit, which.min(sum.fit$bic))

    # --- Fit glfBLUP Model ---
    k_included <- kinship[unique(full_df$G), unique(full_df$G)]
    tryCatch({
      glfBLUP_result <- custom_glfBLUP(data = as.data.frame(full_df), selection = names(best.coefs)[-1], K = k_included)
      glfBLUP_result$rho = cor(glfBLUP_result$preds[match(x$pred.target$G, names(glfBLUP_result$preds))],
                              x$pred.target$pred.target)

      # Save Procrustes glfBLUP predictions
      list.save(glfBLUP_result,
                file = paste0('./data/processed/glfBLUP_timepoints/procrustes/',env, '/', CV,
                              '/', stage, '/glfBLUP_', i, '.RData'))
      rm(glfBLUP_result)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  toc()
}

# ====================================================
# PART C: SPLINE CHARACTERISTIC PREDICTION
# ====================================================
for(settings_i in 1:nrow(settings_df)){
  (env = settings_df[settings_i,]$selected_env)
  (CV = settings_df[settings_i,]$CV)
  (stage = settings_df[settings_i,]$plant_stage)

  n_datasets <- length(list.files(
    paste0('./data/processed/spline_characteristics/', env, '/', CV,'/', stage, '/')))

  tic()
  for(i in 1:n_datasets){
    cat('start run', i, env, CV, stage, '\n')

    # Read pre-calculated spline characteristics (min, max, AUC for spline and derivatives)
    predictor_df <- read_csv(sprintf('./data/processed/spline_characteristics/%s/%s/%s/characteristic_%i.csv',
                                     env, CV, stage, i), show_col_types = FALSE)

    load(sprintf('./data/processed/hyper_folds/%s/hyper_dataset_%i.RData', env, i))

    # --- Wrangle Spline Predictors with Yield Targets ---
    yield_df <- x$data[,c('G', 'plot', 'yield_adjusted_BLUES')] %>%
      mutate(plot = as.character(plot),
             Y = yield_adjusted_BLUES,
             .keep = 'unused')

    full_df <- predictor_df %>%
      mutate(G = as.character(genotype),
             plot = as.character(plotId),
             .keep = 'unused') %>%
      select(-pop, -contains('time')) %>%
      right_join(yield_df, by = c('G', 'plot')) %>%
      select(G, any_of(colnames(predictor_df[, -c(1:3)])), Y)

    # Scale numeric predictors (standardize mean=0, SD=1 across spline features)
    full_df[, -c(1, 20)] <- scale(full_df[, -c(1, 20)])

    # --- Feature Selection via Best Subset Regression ---
    mean_df <- full_df %>%
      group_by(G) %>%
      summarise(across(everything(), mean)) %>%
      column_to_rownames('G')

    best.fit <- regsubsets(Y~ ., data = na.omit(mean_df), nvmax = ncol(mean_df) - 1)
    sum.fit <- summary(best.fit)
    best.coefs <- coef(best.fit, which.min(sum.fit$bic))

    k_included <- kinship[unique(full_df$G), unique(full_df$G)]

    # --- Fit glfBLUP Model ---
    tryCatch({
      glfBLUP_result <- custom_glfBLUP(data = as.data.frame(full_df), selection = names(best.coefs)[-1], K = kinship)
      glfBLUP_result$rho = cor(glfBLUP_result$preds[match(x$pred.target$G, names(glfBLUP_result$preds))],
                              x$pred.target$pred.target)

      # Save spline-based glfBLUP predictions
      list.save(glfBLUP_result,
                file = paste0('./data/processed/glfBLUP_spline/',env, '/', CV, '/', stage, '/glfBLUP_', i, '.RData'))
      rm(glfBLUP_result)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  toc()
}
