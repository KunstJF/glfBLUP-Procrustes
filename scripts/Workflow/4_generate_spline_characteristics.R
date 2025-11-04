# generate splines and spline characteristics

# libraries ----
library(statgenHTP)
library(tictoc)
library(tidyverse)
library(magrittr)

# parallelised ----
# create settings df
settings_df <- expand.grid(CV = c('CV1', 'CV2'),
            selected_env = c('B5IR', 'B2IR', 'DRIP'),
            plant_stage = c('full', 'heading', 'vegetative'),
            stringsAsFactors = FALSE)


# no parallelisation because fitting splines takes a lot of processing power and memory,
# so we work sequential
# iterate over each settings n times
tic()
for(setting_i in 1:nrow(settings_df)){
  CV <- settings_df[setting_i,]$CV
  selected_env <- settings_df[setting_i,]$selected_env
  plant_stage <- settings_df[setting_i,]$plant_stage

  for(i in 51:100){
    cat(sprintf('start run %d with settings: %s %s %s',
            i, CV, selected_env, plant_stage))


  # load in files
  factor_file <- sprintf('./data/processed/factors/%s/timepoint_factors/procrustes/procrustes_%d.RData',
                         selected_env,
                         i)
  load(factor_file)

  # wrangle files to statgenhtp format ----
  factor_dfs <- lapply(procrustes_factor_results[[CV]], function(structure){structure$factor_scores})

  factor_proc_score_df <- bind_rows(factor_dfs, .id = 'timePoint')

  # create statgenhtp compatible df
  factor_df <- factor_proc_score_df %>%
    mutate(G = as.factor(G),
           timePoint = as.POSIXct(ymd(timePoint)),
           pop = 'population1',
           colId = column,
           rowId = row,
           plotId = plot,
           .keep = 'unused')

  # timepoint 5 is heading in B5IR while timepoint 4 is heading in DRIP and B2IR
  if(plant_stage == 'vegetative' & selected_env == 'B5IR'){
    factor_df <- factor_df %>%
      filter(timePoint < unique(factor_df$timePoint)[5])
  }
  if(plant_stage == 'heading' & selected_env == 'B5IR'){
    factor_df <- factor_df %>%
      filter(timePoint <= unique(factor_df$timePoint)[5])
  }
  if(plant_stage == 'vegetative' & selected_env != 'B5IR'){
    factor_df <- factor_df %>%
      filter(timePoint < unique(factor_df$timePoint)[4])
  }
  if(plant_stage == 'heading' & selected_env != 'B5IR'){
    factor_df <- factor_df %>%
      filter(timePoint <= unique(factor_df$timePoint)[4])
  }


  # statgenhtp ----
  ## create splines for each factor
  splines <- list(F1 = NA, F2 = NA)

  for(factor_col in names(splines)){
    splines[[factor_col]] <- fitSplineHDM(inDat = factor_df,
                                          trait = factor_col,
                                          genotype = 'G',
                                          pop = 'pop'
    )
  }

  # interpolate across timepoints for errors and spline estimates
  splines <- lapply(splines, function(factor_spline){
    pred_all.psHDM <- predict(object = factor_spline,
                              newtimes = seq(min(factor_spline$time[["timeNumber"]]),
                                             max(factor_spline$time[["timeNumber"]]),
                                             length.out = 100),
                              pred = list(pop = TRUE, geno = TRUE, plot = TRUE),
                              se = list(pop = FALSE, geno = FALSE, plot = FALSE),
                              trace = FALSE)
  })

  # estimate spline characteristics ----
  # for each factor extract params based on predicted spline, derivative 1 and derivative 2
  extracted_params <- lapply(splines, function(spline){
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
  extracted_dfs_all <- lapply(extracted_params,function(params){
    reduce(params, left_join, by = c('pop','genotype', 'plotId'))
  })

  # combine dataframes from each factor
  all_df <- bind_rows(extracted_dfs_all, .id = 'factor') %>%
    pivot_wider(names_from = 'factor', values_from = min_predictions:AUC_derivatives2)

  # store results ----
  ## splines ----
  spline_result_file <- sprintf('./data/processed/splines/%s/%s/%s/spline_%d.RData',
                                selected_env,
                                CV,
                                plant_stage,
                                i)

  save(splines, file = spline_result_file)

  ## characteristics -----
  spline_result_file <- sprintf('./data/processed/spline_characteristics/%s/%s/%s/characteristic_%d.csv',
                                selected_env,
                                CV,
                                plant_stage,
                                i)

  write_csv(all_df, spline_result_file)

  rm(all_df)
  rm(splines)
  rm(extracted_dfs_all)
  rm(procrustes_factor_results)
  rm(extracted_params)
  rm(factor_df)
  rm(factor_proc_score_df)
  rm(factor_dfs)
  }
}
toc()

# # 1 run example ----
# i = 1
# cat(paste('Start run', i))
# # either B5IR, B2IR or DRIP
# env = 'B5IR'
# # either CV1 or CV2
# CV = 'CV1'
# # either full, heading or vegetative
# plant_stage = 'full'
#
# # load in files
# factor_file <- paste0('./data/processed/factors/B5IR/timepoint_factors/procrustes/procrustes_',
#                       i,
#                       '.RData')
# load(factor_file)
#
# # wrangle files to statgenhtp format ----
# factor_dfs <- lapply(procrustes_factor_results[[CV]], function(structure){structure$factor_scores})
#
# factor_proc_score_df <- bind_rows(factor_dfs, .id = 'timePoint')
#
# # create statgenhtp compatible df
# factor_df <- factor_proc_score_df %>%
#   mutate(G = as.factor(G),
#          timePoint = as.POSIXct(ymd(timePoint)),
#          pop = 'population1',
#          colId = column,
#          rowId = row,
#          plotId = plot,
#          .keep = 'unused')
#
# # statgenhtp ----
# ## create splines for each factor
# splines <- list(F1 = NA, F2 = NA)
#
# for(factor_col in names(splines)){
#   splines[[factor_col]] <- fitSplineHDM(inDat = factor_df,
#                                         trait = factor_col,
#                                         genotype = 'G',
#                                         pop = 'pop'
#   )
# }
#
# # interpolate across timepoints for errors and spline estimates
# splines <- lapply(splines, function(factor_spline){
#   pred_all.psHDM <- predict(object = factor_spline,
#                             newtimes = seq(min(factor_spline$time[["timeNumber"]]),
#                                            max(factor_spline$time[["timeNumber"]]),
#                                            length.out = 100),
#                             pred = list(pop = TRUE, geno = TRUE, plot = TRUE),
#                             se = list(pop = FALSE, geno = FALSE, plot = FALSE),
#                             trace = FALSE)
# })
#
# # estimate spline characteristics ----
# # for each factor extract params based on predicted spline, derivative 1 and derivative 2
# extracted_params <- lapply(splines, function(spline){
#   params <- list()
#   ### fgeno ----
#   params[['fgeno']] <- estimateSplineParameters(spline,
#                                                 estimate = 'predictions',
#                                                 what = c('min', 'max', 'AUC'),
#                                                 fitLevel = 'plot')
#   ### fgenoderiv1 ----
#   params[['fderiv1']] <- estimateSplineParameters(spline,
#                                                   estimate = 'derivatives',
#                                                   what = c('min', 'max', 'AUC'),
#                                                   fitLevel = 'plot')
#   ### fgenodriv2 ----
#   params[['fderiv2']] <- estimateSplineParameters(spline,
#                                                   estimate = 'derivatives2',
#                                                   what = c('min', 'max', 'AUC'),
#                                                   fitLevel = 'plot')
#   return(params)
# })
#
# # wrangling the params ----
# ## using overall ----
# # join the predicted, derivative 1 and derivative for each factor
# extracted_dfs_all <- lapply(extracted_params,function(params){
#   reduce(params, left_join, by = c('pop','genotype', 'plotId'))
# })
#
# # combine dataframes from each factor
# all_df <- bind_rows(extracted_dfs_all, .id = 'factor') %>%
#   pivot_wider(names_from = 'factor', values_from = min_predictions:AUC_derivatives2)
#
# # store results ----
# ## splines ----
# spline_result_file <- paste0('./data/processed/splines/', env,'/', CV,'/',plant_stage,'/','spline_', i,'.RData')
#
# save(spline_result_file, file = spline_result_file)
#
# ## characteristics -----
# spline_result_file <- paste0('./data/processed/spline_characteristics/', env,'/', CV,'/',plant_stage,'/','characteristic_', i,'.csv')
#
# write_csv(all_df, spline_result_file)
