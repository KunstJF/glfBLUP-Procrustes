# generate factors per fold

# library ----
library(stringr)
library(tictoc)
library(gfBLUP)
library(doParallel)

# local functions ----
for (file.name in list.files("functions")) {
  source(paste0("functions/", file.name))
}
selected_env = 'DRIP'
# set up cluster
n_datasets <- length(list.files(sprintf('./data/processed/hyper_folds/%s/', selected_env)))
n_cores <- parallel::detectCores() - 2
cl <- parallel::makeCluster(n_cores, outfile = "logs/varimax_factor_loadings_generation.txt")
doParallel::registerDoParallel(cl)

# perform factor generation varimax ----
tic()
foreach::foreach(i = 1:n_datasets, .packages = c("dplyr", "stringr")) %dopar% {
  # load in file
  data_file <- sprintf('./data/processed/hyper_folds/%s/hyper_dataset_%d.RData', selected_env, i)

  load(data_file)
  non_hs_cols <- c('G', 'row', 'column', 'plot', 'yield_adjusted_BLUES', 'Y')
  hs_cols <- !colnames(x$data) %in% non_hs_cols

  # generate factor loadings
  timepoints <- unique(str_split(colnames(x$data)[hs_cols], '_', simplify = TRUE)[, 2])

  generate_factor(df = x$data[, c('G', colnames(x$data)[hs_cols])], timepoint = timepoints[4], m = 2)

  factor_result_cv2 <- lapply(timepoints, generate_factor, df = x$data[, c('G', colnames(x$data)[hs_cols])], m = 2)
  names(factor_result_cv2) <- timepoints

  factor_dfs_cv2 <- lapply(factor_result_cv2, function(result){
    # add field and yield information
    temp_df <- dplyr::bind_cols(x$data[,c(1:4)],result$factor_scores[, -1])

    # join to full data, witholding hyperspectral data
    full_df <- dplyr::left_join(x$data[non_hs_cols], temp_df, by = c('G', 'row', 'column', 'plot'))

    # join to full data, witholding hyperspectral data
    full_df <- dplyr::left_join(x$data[non_hs_cols], temp_df, by = c('G', 'row', 'column', 'plot'))
  })

  # filter genotypes in data for CV1
  df_cv1 <- x$data[which(x$data$G %in% x$train.set), c('G', colnames(x$data)[hs_cols])]

  factor_result_cv1 <- lapply(timepoints, generate_factor, df = df_cv1, m =2)
  names(factor_result_cv1) <- timepoints

  factor_dfs_cv1 <- lapply(factor_result_cv1, function(result){
    # add field and yield information
    temp_df <- dplyr::bind_cols(x$data[which(x$data$G %in% x$train.set), c(1:4)],result$factor_scores[, -1])

    # join to full data, witholding hyperspectral data
    full_df <- dplyr::left_join(x$data[non_hs_cols], temp_df, by = c('G', 'row', 'column', 'plot'))

    return(full_df)
  })


  # bind results
  factor_result <- list(cv1_factor_sctructure = factor_result_cv1,
                        cv1_factor_dfs = factor_dfs_cv1,
                        cv2_factor_structure = factor_result_cv2,
                        cv2_factor_dfs = factor_dfs_cv2)

  # store results
  factor_loadings_file <- sprintf('./data/processed/factors/%s/timepoint_factors/varimax/factors_dataset_%d.RData',
                                  selected_env, i)

  save(factor_result, file = factor_loadings_file)
  }
doParallel::stopImplicitCluster()
parallel::stopCluster(cl)
toc()

# test a single run varimax ----
# i <- 1
# # load in file
# data_file <- paste0('./data/processed/hyper_folds/B5IR/hyper_dataset_',
#                     i,
#                     '.RData')
# load(data_file)
# ncol(x$data)
# colnames(x$data)
# # generate factor loadings
# timepoints <- unique(str_split(colnames(x$data)[-c(1:4, 625:626)], '_', simplify = TRUE)[, 2])
#
# generate_factor(df = x$data[, -c(2:4, 625:626)], timepoint = timepoints[4], m = 2)
#
# factor_result_cv2 <- lapply(timepoints, generate_factor, df = x$data[, -c(2:4, 625:626)], m = 2)
# names(factor_result_cv2) <- timepoints
#
# factor_dfs_cv2 <- lapply(factor_result_cv2, function(result){
#   # add field and yield information
#   temp_df <- dplyr::bind_cols(x$data[,c(1:4)],result$factor_scores[, -1])
#
#   # join to full data, witholding hyperspectral data
#   full_df <- dplyr::left_join(x$data[c(1:4, 625:626)], temp_df, by = c('G', 'row', 'column', 'plot'))
#
#   # join to full data, witholding hyperspectral data
#   full_df <- dplyr::left_join(x$data[c(1:4, 625:626)], temp_df, by = c('G', 'row', 'column', 'plot'))
# })
#
# # filter genotypes in data for CV1
# df_cv1 <- x$data[which(x$data$G %in% x$train.set), -c(2:4, 625:626)]
#
# factor_result_cv1 <- lapply(timepoints, generate_factor, df = df_cv1, m =2)
# names(factor_result_cv1) <- timepoints
#
# factor_dfs_cv1 <- lapply(factor_result_cv1, function(result){
#   # add field and yield information
#   temp_df <- dplyr::bind_cols(x$data[which(x$data$G %in% x$train.set),
#                                      c(1:4)],result$factor_scores[, -1])
#
#   # join to full data, witholding hyperspectral data
#   full_df <- dplyr::left_join(x$data[c(1:4, 625:626)], temp_df, by = c('G', 'row', 'column', 'plot'))
#
#   return(full_df)
# })
#
#
# # bind results
# factor_result <- list(cv1_factor_sctructure = factor_result_cv1,
#                       cv1_factor_dfs = factor_dfs_cv1,
#                       cv2_factor_structure = factor_result_cv2,
#                       cv2_factor_dfs = factor_dfs_cv2)

# store results
# factor_loadings_file <- paste0('./data/processed/factors/timepoint_factors/varimax/factors_dataset_',
#                                i,
#                                '.RData')
#
# save(factor_result, file = factor_loadings_file)



