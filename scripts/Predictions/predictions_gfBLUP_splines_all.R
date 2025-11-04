# CIMMYT spline gfBLUP
# author: Jonathan Kunst

# libraries ----
library(tidyverse)
library(statgenHTP)
library(rrBLUP)
library(wesanderson)
library(corrplot)
library(leaps)
library(gfBLUP)
library(tidymodels)
library(caret)
library(doParallel)
library(tictoc)

# seed ----
set.seed(1314)

# data ----
## spline data ----
# heading starts 25 february
# grain filling starts 10 march
# load splines with all timepoints
fnames_all <- str_subset(list.files('./data/processed/splines/all_data/', full.names = TRUE), 'predsplinedatall')
fnames_all <- str_subset(fnames_all, 'mid_timepoint_5_2014_Bed5IR')
splines_all <- lapply(fnames_all, read_rds)


## yield data ----
yield_pseudo <- read_csv('./data/processed/2014_yield_pseudo_blue.csv')

## read in kinship matrix ----
# genos <- read_csv('./data/raw/Krause_et_al_2018_Genotypes.csv')
# kinship <- genos %>%
#   column_to_rownames('GID') %>%
#   A.mat(., n.core = 15)
#
# save(kinship, file = './data/processed/kinship.RData')

load('./data/processed/kinship.RData')

## folds ----
load('./data/processed/hyper_folds.RData')

# wrangling ----
selected_year <- 2014
selected_env <- 'B5IR'

check_genos <- table(yield_pseudo$GID)[1:2]

yield_pseudo_select <- yield_pseudo %>%
  mutate(genotype = as.character(GID),
         plotId = as.factor(plot), .keep = 'unused')  %>%
  filter(cycle == selected_year &
           env == selected_env &
           !genotype %in% names(check_genos) &
           genotype %in% rownames(kinship))

# extract parameters for all timepoints ----
# for each factor extract params based on predicted spline, derivative 1 and derivative 2
extracted_params_all <- lapply(splines_all, function(spline){
  params <- list()
  ### fgeno ----
  params[['fgeno']] <- estimateSplineParameters(spline,
                                                estimate = 'predictions',
                                                what = c('min', 'max', 'AUC'),
                                                fitLevel = 'plot')
  ### fgenoderiv1 ----
  params[['fderiv1']] <- estimateSplineParameters(spline,
                                                  estimate = 'derivatives',
                                                  what = c('min', 'max', 'AUC'),
                                                  fitLevel = 'plot')
  ### fgenodriv2 ----
  params[['fderiv2']] <- estimateSplineParameters(spline,
                                                  estimate = 'derivatives2',
                                                  what = c('min', 'max', 'AUC'),
                                                  fitLevel = 'plot')
  return(params)
})

# wrangling the params ----
## using overall ----
# join the predicted, derivative 1 and derivative for each factor
extracted_dfs_all <- lapply(extracted_params_all,function(params){
  reduce(params, left_join, by = c('pop','genotype', 'plotId'))
})

# provide factor names to df list
names(extracted_dfs_all) <- str_extract(fnames_all, 'factor[0-9]')

# combine dataframes from each factor
all_df <- bind_rows(extracted_dfs_all, .id = 'factor') %>%
  pivot_wider(names_from = 'factor', values_from = min_predictions:AUC_derivatives2)

all_yield <- all_df %>%
  left_join(yield_pseudo_select, by = c('plotId', 'genotype')) %>%
  select(!contains(c('min_timePoint', 'max_timePoint', 'Pop', 'timeNumber', 'cycle',
                     'rep', 'subblock', 'row', 'column', 'entry', 'Cross.Name',
                     'trial.name', 'env', 'trial.no'))) %>%
  select(!'gy') %>%
  drop_na()

# first 2 columns are genotype and plotid
all_yield[, -c(1:2)] <- scale(all_yield[, -c(1:2)])

# plot level data
subset_select_df <- all_yield %>%
  unite(unique_geno, genotype:plotId, sep = '_') %>%
  column_to_rownames('unique_geno')

# mean data
subset_mean_df <- all_yield %>%
  select(-plotId) %>%
  group_by(genotype) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames('genotype')

# Melsen folds ----
phenomic_selection <- function(fold, df){
  train_index <- which(rownames(df) %in% fold$train)
  test_index <- which(!rownames(df) %in% fold$test)

  # train
  best.fit <- regsubsets(adj_gy ~ ., data = df[train_index,], nvmax = ncol(df) - 1)
  sum.fit <- summary(best.fit)
  best.coefs <- coef(best.fit, which.min(sum.fit$bic))

  # test
  test_mat <- model.matrix(adj_gy ~ ., data = df[test_index,], nvmax = ncol(df) - 1)
  predicts <- test_mat[, names(best.coefs)] %*% best.coefs

  # mse
  mse <- mean((df[test_index,]$adj_gy - predicts)^2)

  # correlations
  rho <- cor(df[test_index,]$adj_gy, predicts)

  return(c(best.coefs,
           'train_bic' = min(sum.fit$bic),
           'train_adjr2' = max(sum.fit$adjr2),
           'mse'=mse,
           'rho' = rho))
}

## cluster creation ----
cl <- makeCluster(detectCores(), type = 'PSOCK')
registerDoParallel(cl)
clusterExport(cl, varlist = c('folds', 'subset_mean_df'))
clusterEvalQ(cl, {library(leaps)})

## run code paralell
tic()
phenomic_preds_melsen <- parLapply(cl, folds, df = subset_mean_df, phenomic_selection)
toc()
stopCluster(cl)

write_rds(phenomic_preds_melsen, './data/processed/gfBLUP_predictions/alltimepoints__subset_prediction_melsenfolds.rds')

# get best combination of predictors
best.fit <- regsubsets(adj_gy ~ ., data = subset_mean_df, nvmax = ncol(subset_mean_df) - 1)
sum.fit <- summary(best.fit)
best.coefs <- coef(best.fit, which.min(sum.fit$bic))
best_preds <- names(best.coefs)[-1]

best_preds

# gfBLUP ----
# wrangle data so it aligns with gfBLUP restrictions
data_gfblup <- all_yield %>%
  mutate(Y = adj_gy,
         G = genotype,
         .keep = 'unused') %>%
  select(G, all_of(best_preds), Y) %>%
  as.data.frame(.)

# test gfBLUP ----
# cv 1 ----
gfBLUP_cv1 <- function(fold, df, K, yield_mean){
  test_data <- df
  test_data[which(test_data$G %in% fold$test), -1] <- NA

  test_gfblup <- gfBLUP(test_data, K = K)

  pred_y <- data.frame(G = names(test_gfblup$preds),
                       predictions = as.numeric(test_gfblup$preds))

  result_df <- merge(yield_mean, pred_y, by = 'G')
  rho <- cor(result_df$mean_y, result_df$predictions)

  rm(result_df)
  rm(test_data)

  return(
    c(rho = rho,
      gfblup_result = test_gfblup)
  )
}

yield_mean <- data_gfblup %>%
  group_by(G) %>%
  summarise(mean_y = mean(Y))

# how long does 1 fold take?
tic()
test_gfblup_cv1 <- gfBLUP_cv1(folds[[5]], df = data_gfblup, K = kinship, yield_mean = yield_mean)
toc()
test_gfblup_cv1$rho
## cluster creation ----
cl <- makeCluster(detectCores(), type = 'PSOCK')
registerDoParallel(cl)
clusterExport(cl, varlist = c('folds', 'data_gfblup', 'kinship', 'yield_mean'))
clusterEvalQ(cl, {library(gfBLUP)})

## run code paralell
tic()
gfblup_cv1_melsen <- parLapply(cl, folds, df = data_gfblup, K = kinship, yield_mean = yield_mean, fun = gfBLUP_cv1)
toc()
stopCluster(cl)

write_rds(gfblup_cv1_melsen, './data/processed/gfBLUP_predictions/alltimepoints_Bed2IR_gfBLUP_cv1_melsen_folds.rds')

# cv 2----
yield_mean <- data_gfblup %>%
  group_by(G) %>%
  summarise(mean_y = mean(Y))

gfBLUP_cv2 <- function(fold, df, K, yield_mean){
  test_data <- df
  test_data[which(test_data$G %in% fold$test), ]$Y <- NA

  test_gfblup <- gfBLUP(test_data, K = K)

  pred_y <- data.frame(G = names(test_gfblup$preds),
                       predictions = as.numeric(test_gfblup$preds))

  result_df <- merge(yield_mean, pred_y, by = 'G')
  rho <- cor(result_df$mean_y, result_df$predictions)

  rm(result_df)
  rm(test_data)

  return(
    c(rho = rho,
      gfblup_result = test_gfblup)
  )
}

# how long does 1 run take?
tic()
test_gfBLUP_cv2 <- gfBLUP_cv2(fold = folds[[5]], df = data_gfblup, K = kinship, yield_mean = yield_mean)
toc()
test_gfBLUP_cv2$rho
cl <- makeCluster(detectCores(), type = 'PSOCK')
registerDoParallel(cl)
clusterExport(cl, varlist = c('folds', 'data_gfblup', 'kinship', 'yield_mean'))
clusterEvalQ(cl, {library(gfBLUP)})

tic()
gfBLUP_cv2_melsen <- parLapply(cl, folds, df = data_gfblup, K = kinship, yield_mean = yield_mean,  fun = gfBLUP_cv2)
toc()
stopCluster(cl)

write_rds(gfBLUP_cv2_melsen, './data/processed/gfBLUP_predictions/alltimepoints_B5IR_gfBLUP_cv2_melsenfolds.rds')
