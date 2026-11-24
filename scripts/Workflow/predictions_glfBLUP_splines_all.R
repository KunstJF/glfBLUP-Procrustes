################################################################################
# script name: CIMMYT spline glfBLUP
################################################################################
# ----------------------------------------------------
# 1. Load Required Libraries & Set Seed
# ----------------------------------------------------
library(tidyverse)
library(statgenHTP)
library(rrBLUP)
library(wesanderson)
library(corrplot)
library(leaps)
library(glfBLUP)
library(tidymodels)
library(caret)
library(doParallel)
library(tictoc)

set.seed(1314)

# ----------------------------------------------------
# 2. Data Ingestion & Setup
# ----------------------------------------------------
# --- Load Spline Data ---
# Select spline predictions for target trial (2014, Bed5IR) across all timepoints
fnames_all <- str_subset(list.files('./data/processed/splines/all_data/', full.names = TRUE), 'predsplinedatall')
fnames_all <- str_subset(fnames_all, 'mid_timepoint_5_2014_Bed5IR')
splines_all <- lapply(fnames_all, read_rds)

# --- Load Target & Kinship Data ---
yield_pseudo <- read_csv('./data/processed/2014_yield_pseudo_blue.csv')

# Load pre-computed additive genomic relationship matrix (kinship)
load('./data/processed/kinship.RData')

# Load predefined cross-validation split indices (Melsen folds)
load('./data/processed/hyper_folds.RData')

# ----------------------------------------------------
# 3. Data Wrangling & Parameter Extraction
# ----------------------------------------------------

# change environment and or year here
selected_year <- 2014
selected_env <- 'B5IR'

# Identify and filter out low-frequency test genotypes
check_genos <- table(yield_pseudo$GID)[1:2]

yield_pseudo_select <- yield_pseudo %>%
  mutate(genotype = as.character(GID),
         plotId = as.factor(plot), .keep = 'unused')  %>%
  filter(cycle == selected_year &
           env == selected_env &
           !genotype %in% names(check_genos) &
           genotype %in% rownames(kinship))

# --- Extract Spline Characteristics ---
# Calculate parameter estimates (min, max, AUC) for predictions, 1st derivatives, and 2nd derivatives
extracted_params_all <- lapply(splines_all, function(spline){
  params <- list()
  ### Trajectory Predictions ---
  params[['fgeno']] <- estimateSplineParameters(spline,
                                                estimate = 'predictions',
                                                what = c('min', 'max', 'AUC'),
                                                fitLevel = 'plot')
  ### Growth Rate (1st Derivative) ---
  params[['fderiv1']] <- estimateSplineParameters(spline,
                                                  estimate = 'derivatives',
                                                  what = c('min', 'max', 'AUC'),
                                                  fitLevel = 'plot')
  ### Growth Acceleration (2nd Derivative) ---
  params[['fderiv2']] <- estimateSplineParameters(spline,
                                                  estimate = 'derivatives2',
                                                  what = c('min', 'max', 'AUC'),
                                                  fitLevel = 'plot')
  return(params)
})

# --- Merge & Standardize Trait Datasets ---
extracted_dfs_all <- lapply(extracted_params_all, function(params){
  reduce(params, left_join, by = c('pop','genotype', 'plotId'))
})

names(extracted_dfs_all) <- str_extract(fnames_all, 'factor[0-9]')

# Combine factor datasets to wide format and join adjusted grain yield targets
all_df <- bind_rows(extracted_dfs_all, .id = 'factor') %>%
  pivot_wider(names_from = 'factor', values_from = min_predictions:AUC_derivatives2)

all_yield <- all_df %>%
  left_join(yield_pseudo_select, by = c('plotId', 'genotype')) %>%
  select(!contains(c('min_timePoint', 'max_timePoint', 'Pop', 'timeNumber', 'cycle',
                     'rep', 'subblock', 'row', 'column', 'entry', 'Cross.Name',
                     'trial.name', 'env', 'trial.no'))) %>%
  select(!'gy') %>%
  drop_na()

# Scale secondary predictor variables (mean = 0, SD = 1)
all_yield[, -c(1:2)] <- scale(all_yield[, -c(1:2)])

# Construct plot-level dataset
subset_select_df <- all_yield %>%
  unite(unique_geno, genotype:plotId, sep = '_') %>%
  column_to_rownames('unique_geno')

# Construct genotype-level averaged dataset
subset_mean_df <- all_yield %>%
  select(-plotId) %>%
  group_by(genotype) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames('genotype')

# ----------------------------------------------------
# 4. Phenomic Selection via Best Subsets Regression
# ----------------------------------------------------
# Core evaluation function for phenomic selection across CV folds
phenomic_selection <- function(fold, df){
  train_index <- which(rownames(df) %in% fold$train)
  test_index <- which(!rownames(df) %in% fold$test)

  # Model training and optimal feature selection via BIC minimization
  best.fit <- regsubsets(adj_gy ~ ., data = df[train_index,], nvmax = ncol(df) - 1)
  sum.fit <- summary(best.fit)
  best.coefs <- coef(best.fit, which.min(sum.fit$bic))

  # Model prediction on test set
  test_mat <- model.matrix(adj_gy ~ ., data = df[test_index,], nvmax = ncol(df) - 1)
  predicts <- test_mat[, names(best.coefs)] %*% best.coefs

  # Compute prediction diagnostics (MSE & Pearson's rho)
  mse <- mean((df[test_index,]$adj_gy - predicts)^2)
  rho <- cor(df[test_index,]$adj_gy, predicts)

  return(c(best.coefs,
           'train_bic' = min(sum.fit$bic),
           'train_adjr2' = max(sum.fit$adjr2),
           'mse' = mse,
           'rho' = rho))
}

# --- Parallel Phenomic Selection ---
cl <- makeCluster(detectCores(), type = 'PSOCK')
registerDoParallel(cl)
clusterExport(cl, varlist = c('folds', 'subset_mean_df'))
clusterEvalQ(cl, {library(leaps)})

tic()
phenomic_preds_melsen <- parLapply(cl, folds, df = subset_mean_df, phenomic_selection)
toc()
stopCluster(cl)

# Save phenomic selection outputs
write_rds(phenomic_preds_melsen, './data/processed/glfBLUP_predictions/alltimepoints_subset_prediction_melsenfolds.rds')

# Determine global optimal predictors across all genotype means
best.fit <- regsubsets(adj_gy ~ ., data = subset_mean_df, nvmax = ncol(subset_mean_df) - 1)
sum.fit <- summary(best.fit)
best.coefs <- coef(best.fit, which.min(sum.fit$bic))
best_preds <- names(best.coefs)[-1]

best_preds

# ----------------------------------------------------
# 5. glfBLUP
# ----------------------------------------------------
# Format dataset according to glfBLUP requirements
data_glfblup <- all_yield %>%
  mutate(Y = adj_gy,
         G = genotype,
         .keep = 'unused') %>%
  select(G, all_of(best_preds), Y) %>%
  as.data.frame(.)

yield_mean <- data_glfblup %>%
  group_by(G) %>%
  summarise(mean_y = mean(Y))

# --- Cross-Validation 1 (CV1: Unobserved Genotypes) ---
# Evaluates prediction performance when both secondary traits and target traits are missing in validation set
glfBLUP_cv1 <- function(fold, df, K, yield_mean){
  test_data <- df
  # Mask ALL variables for testing genotypes
  test_data[which(test_data$G %in% fold$test), -1] <- NA

  test_glfblup <- glfBLUP(test_data, K = K)

  pred_y <- data.frame(G = names(test_glfblup$preds),
                       predictions = as.numeric(test_glfblup$preds))

  result_df <- merge(yield_mean, pred_y, by = 'G')
  rho <- cor(result_df$mean_y, result_df$predictions)

  rm(result_df, test_data)

  return(c(rho = rho, glfblup_result = test_glfblup))
}

# Single run benchmark CV1
tic()
test_glfblup_cv1 <- glfBLUP_cv1(folds[[5]], df = data_glfblup, K = kinship, yield_mean = yield_mean)
toc()
test_glfblup_cv1$rho

# Parallel CV1 Execution
cl <- makeCluster(detectCores(), type = 'PSOCK')
registerDoParallel(cl)
clusterExport(cl, varlist = c('folds', 'data_glfblup', 'kinship', 'yield_mean'))
clusterEvalQ(cl, {library(glfBLUP)})

tic()
glfblup_cv1_melsen <- parLapply(cl, folds, df = data_glfblup, K = kinship, yield_mean = yield_mean, fun = glfBLUP_cv1)
toc()
stopCluster(cl)

write_rds(glfblup_cv1_melsen, './data/processed/glfBLUP_predictions/alltimepoints_Bed2IR_glfBLUP_cv1_melsen_folds.rds')

# --- Cross-Validation 2 ---
# Evaluates prediction performance when secondary traits are known, but target trait (Y) is missing
glfBLUP_cv2 <- function(fold, df, K, yield_mean){
  test_data <- df
  # Mask ONLY the primary target variable (Y) for testing genotypes
  test_data[which(test_data$G %in% fold$test), ]$Y <- NA

  test_glfblup <- glfBLUP(test_data, K = K)

  pred_y <- data.frame(G = names(test_glfblup$preds),
                       predictions = as.numeric(test_glfblup$preds))

  result_df <- merge(yield_mean, pred_y, by = 'G')
  rho <- cor(result_df$mean_y, result_df$predictions)

  rm(result_df, test_data)

  return(c(rho = rho, glfblup_result = test_glfblup))
}

# Single run benchmark (Fold 5)
tic()
test_glfBLUP_cv2 <- glfBLUP_cv2(fold = folds[[5]], df = data_glfblup, K = kinship, yield_mean = yield_mean)
toc()
test_glfBLUP_cv2$rho

# Parallel CV2 Execution
cl <- makeCluster(detectCores(), type = 'PSOCK')
registerDoParallel(cl)
clusterExport(cl, varlist = c('folds', 'data_glfblup', 'kinship', 'yield_mean'))
clusterEvalQ(cl, {library(glfBLUP)})

tic()
glfBLUP_cv2_melsen <- parLapply(cl, folds, df = data_glfblup, K = kinship, yield_mean = yield_mean, fun = glfBLUP_cv2)
toc()
stopCluster(cl)

# Save CV2 outputs
write_rds(glfBLUP_cv2_melsen, './data/processed/glfBLUP_predictions/alltimepoints_B5IR_glfBLUP_cv2_melsenfolds.rds')
