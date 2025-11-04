# generate gfblup predictions

# libraries ----
library(leaps)
library(gfBLUP)
library(magrittr)
library(tibble)
library(rlist)
library(tictoc)
library(tidyverse)


# local functions ----
for (file.name in list.files("functions")) {
  source(paste0("functions/", file.name))
}

# kinship used across ----
K <- load('./data/processed/kinship.RData')

# create settings df ----
settings_df <- expand.grid(CV = c('CV1', 'CV2'),
                           selected_env = c('B5IR', 'B2IR', 'DRIP'),
                           plant_stage = c('full', 'heading', 'vegetative'),
                           stringsAsFactors = FALSE)


# naive factor scores ----
for(settings_i in 1:nrow(settings_df)){
  (CV = settings_df[settings_i,]$CV)
  (env = settings_df[settings_i,]$selected_env)
  (stage = settings_df[settings_i,]$plant_stage)

  n_datasets <- length(list.files(
    sprintf('./data/processed/factors/%s/timepoint_factors/varimax/', env)))

  tic()
  for(i in 1:100){
    cat('start run', i, env, CV, stage, '\n')

    load(sprintf('./data/processed/factors/%s/timepoint_factors/varimax/factors_dataset_%d.RData', env, i))

    load(sprintf('./data/processed/hyper_folds/%s/hyper_dataset_%i.RData', env, i))

    if(env == 'B5IR'){
      heading_i <- 5
    }else{
      heading_i <- 4
    }

    ### wrangle factors for setting ----
    cv_loc <- str_detect(names(factor_result), paste0(tolower(CV), '_factor_dfs'))

    factor_scores <- lapply(factor_result[cv_loc], function(setting){bind_rows(setting, .id = 'timepoint')})

    timepoints <- names(factor_result[cv_loc][[1]])

    # determine which timepoints to retain based on plant stages used
    if(stage == 'heading'){
      timepoints <- timepoints[1:heading_i]
    }else if(stage == 'vegetative'){
      timepoints <- timepoints[1:(heading_i -1)]
    }

    # structure secondary factor and focal trait df
    fac_df <- factor_scores[[1]] %>%
      filter(timepoint %in% timepoints) %>%
      pivot_wider(names_from = timepoint, values_from = F1:F2) %>%
      select(G, contains('F'), yield_adjusted_BLUES, Y) %>%
      mutate(Y = yield_adjusted_BLUES, .keep = 'unused')

    ### wrangle yield and secondary trait data together ----
    # scale predictor variables ----
    # create means df for subset selection
    mean_df <- fac_df %>%
      group_by(G) %>%
      summarise(across(everything(), mean)) %>%
      column_to_rownames('G')

    # subset selection ----
    best.fit <- regsubsets(Y~ ., data = na.omit(mean_df), nvmax = ncol(mean_df) - 1)
    sum.fit <- summary(best.fit)
    best.coefs <- coef(best.fit, which.min(sum.fit$bic))

    # gfBLUP ----
    k_included <- kinship[unique(fac_df$G), unique(fac_df$G)]
    tryCatch({
      gfBLUP_result <- custom_gfBLUP(data = as.data.frame(fac_df), selection = names(best.coefs)[-1], K = k_included)
      gfBLUP_result$rho = cor(gfBLUP_result$preds[match(x$pred.target$G, names(gfBLUP_result$preds))],
                              x$pred.target$pred.target)
      # store results ----
      list.save(gfBLUP_result,
                file = paste0('./data/processed/gfBLUP_timepoints/varimax/',env, '/', CV,
                              '/', stage, '/gfBLUP_', i, '.RData'))
      rm(gfBLUP_result)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  toc()
}

# rotated factor scores ----
for(settings_i in 1:nrow(settings_df)){
  (CV = settings_df[settings_i,]$CV)
  (env = settings_df[settings_i,]$selected_env)
  (stage = settings_df[settings_i,]$plant_stage)

  n_datasets <- length(list.files(
    sprintf('./data/processed/factors/%s/timepoint_factors/varimax/', env)))

  tic()
  # just 100 because splines took a long time so stopped at 100 sets
  for(i in 1:100){
    cat('start run', i, env, CV, stage, '\n')

    load(sprintf('./data/processed/factors/%s/timepoint_factors/procrustes/procrustes_%d.RData', env, i))

    load(sprintf('./data/processed/hyper_folds/%s/hyper_dataset_%i.RData', env, i))

    if(env == 'B5IR'){
      heading_i <- 5
    }else{
      heading_i <- 4
    }

    ### wrangle factors for setting ----
    factor_scores <- lapply(procrustes_factor_results[[CV]], function(setting){setting$factor_scores})

    factor_scores_df <- bind_rows(factor_scores, .id = 'timepoint')

    timepoints <- names(factor_scores)

    # determine which timepoints to retain based on plant stages used
    if(stage == 'heading'){
      timepoints <- timepoints[1:heading_i]
    }else if(stage == 'vegetative'){
      timepoints <- timepoints[1:(heading_i -1)]
    }

    # structure secondary factor and focal trait df
    fac_df <- factor_scores_df %>%
      filter(timepoint %in% timepoints) %>%
      pivot_wider(names_from = timepoint, values_from = F1:F2)

    full_df <- x$data %>%
      left_join(fac_df, by = c('G', 'row', 'column', 'plot')) %>%
      select(G, contains('F'), yield_adjusted_BLUES.x) %>%
      mutate(Y = yield_adjusted_BLUES.x, .keep = 'unused')

    ### wrangle yield and secondary trait data together ----
    # scale predictor variables ----
    # create means df for subset selection
    mean_df <- full_df %>%
      group_by(G) %>%
      summarise(across(everything(), mean)) %>%
      column_to_rownames('G')

    # subset selection ----
    best.fit <- regsubsets(Y~ ., data = na.omit(mean_df), nvmax = ncol(mean_df) - 1)
    sum.fit <- summary(best.fit)
    best.coefs <- coef(best.fit, which.min(sum.fit$bic))

    # gfBLUP ----
    k_included <- kinship[unique(full_df$G), unique(full_df$G)]
    tryCatch({
      gfBLUP_result <- custom_gfBLUP(data = as.data.frame(full_df), selection = names(best.coefs)[-1], K = k_included)
      gfBLUP_result$rho = cor(gfBLUP_result$preds[match(x$pred.target$G, names(gfBLUP_result$preds))],
                              x$pred.target$pred.target)
      # store results ----
      list.save(gfBLUP_result,
                file = paste0('./data/processed/gfBLUP_timepoints/procrustes/',env, '/', CV,
                              '/', stage, '/gfBLUP_', i, '.RData'))
      rm(gfBLUP_result)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  toc()
}

# splines characteristics ----
for(settings_i in 1:nrow(settings_df)){
  (env = settings_df[settings_i,]$selected_env)
  (CV = settings_df[settings_i,]$CV)
  (stage = settings_df[settings_i,]$plant_stage)

  n_datasets <- length(list.files(
    paste0('./data/processed/spline_characteristics/', env, '/', CV,'/', stage, '/')))

tic()
for(i in 1:n_datasets){
  cat('start run', i, env, CV, stage, '\n')
  predictor_df <- read_csv(sprintf('./data/processed/spline_characteristics/%s/%s/%s/characteristic_%i.csv',
            env, CV, stage, i), show_col_types = FALSE)

  load(sprintf('./data/processed/hyper_folds/%s/hyper_dataset_%i.RData', env, i))

  # wrangle yield and secondary trait data together ----
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

  # scale predictor variables ----
  full_df[, -c(1, 20)] <- scale(full_df[, -c(1, 20)])

  # create means df for subset selection
  mean_df <- full_df %>%
    group_by(G) %>%
    summarise(across(everything(), mean)) %>%
    column_to_rownames('G')

  # subset selection ----
  best.fit <- regsubsets(Y~ ., data = na.omit(mean_df), nvmax = ncol(mean_df) - 1)
  sum.fit <- summary(best.fit)
  best.coefs <- coef(best.fit, which.min(sum.fit$bic))

  k_included <- kinship[unique(full_df$G), unique(full_df$G)]

  # gfBLUP ----
  tryCatch({
    gfBLUP_result <- custom_gfBLUP(data = as.data.frame(full_df), selection = names(best.coefs)[-1], K = kinship)
    gfBLUP_result$rho = cor(gfBLUP_result$preds[match(x$pred.target$G, names(gfBLUP_result$preds))],
                            x$pred.target$pred.target)
    # store results ----
    list.save(gfBLUP_result,
              file = paste0('./data/processed/gfBLUP_spline/',env, '/', CV, '/', stage, '/gfBLUP_', i, '.RData'))
    rm(gfBLUP_result)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
toc()
}


# single example ----
# i = 1
# env = 'B5IR'
# CV = 'CV1'
# stage = 'full'
# predictor_df <- read_csv(paste0('./data/processed/spline_characteristics/',env, '/',CV,'/',stage, '/characteristic_', i, '.csv'))
# load(paste0('./data/processed/hyper_folds/B5IR/hyper_dataset_', i, '.RData'))
#
# # wrangle yield and secondary trait data together ----
# yield_df <- x$data[,c(1, 4, 625)] %>%
#   mutate(plot = as.character(plot),
#          Y = yield_adjusted_BLUES,
#          .keep = 'unused')
#
# full_df <- predictor_df %>%
#   mutate(G = as.character(genotype),
#          plot = as.character(plotId),
#          .keep = 'unused') %>%
#   select(-pop, -contains('time')) %>%
#   right_join(yield_df, by = c('G', 'plot')) %>%
#   select(G, any_of(colnames(predictor_df[, -c(1:3)])), Y)
#
# # scale predictor variables ----
# full_df[, -c(1, 20)] <- scale(full_df[, -c(1, 20)])
#
# # create means df for subset selection
# mean_df <- full_df %>%
#   group_by(G) %>%
#   summarise(across(everything(), mean)) %>%
#   column_to_rownames('G')
#
# # subset selection ----
# best.fit <- regsubsets(Y~ ., data = na.omit(mean_df), nvmax = ncol(mean_df) - 1)
# sum.fit <- summary(best.fit)
# best.coefs <- coef(best.fit, which.min(sum.fit$bic))
#
# k_included <- kinship[unique(full_df$G), unique(full_df$G)]
#
# # gfBLUP ----
# gfBLUP_result <- gfBLUP(as.data.frame(full_df), selection = names(best.coefs)[-1], K = kinship)
# gfBLUP_result$rho = cor(gfBLUP_result$preds[match(x$pred.target$G, names(gfBLUP_result$preds))],
#                                      x$pred.target$pred.target)
#
#
# # store results ----
# list.save(gfBLUP_result,
#           file = paste0('./data/processed/gfBLUP_spline/',env, '/', CV, '/', stage, '/gfBLUP_', i, '.RData'))


