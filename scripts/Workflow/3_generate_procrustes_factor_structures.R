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

selected_env <- 'B2IR'
# focal times:
# B5IR: 4
# DRIP, B2IR: 3
focal_time <- 3

# set up cluster
n_datasets <- length(list.files(sprintf('./data/processed/factors/%s/timepoint_factors/varimax/', selected_env)))
n_cores <- parallel::detectCores() - 2
cl <- parallel::makeCluster(n_cores, outfile = "logs/varimax_factor_loadings_generation.txt")
doParallel::registerDoParallel(cl)

# perform factor generation varimax ----
tic()
foreach::foreach(i = 1:n_datasets, .packages = c("stringr", 'magrittr', 'dplyr')) %dopar% {

  # load in files
  factor_file <- sprintf('./data/processed/factors/%s/timepoint_factors/varimax/factors_dataset_%d.RData',
                         selected_env,
                         i)
  load(factor_file)

  # load in hyperspectral file
  data_file <- sprintf('./data/processed/hyper_folds/%s/hyper_dataset_%d.RData',
                       selected_env,
                       i)
  load(data_file)

  # select focal timepoint
  timepoints <- names(factor_result$cv2_factor_structure)
  head_timepoint <- timepoints[focal_time]

  target_structure <- factor_result$cv2_factor_structure[[focal_time]]

  # Procrustean smoothing ----
  # create procustes rotated factors for each timepoint
  proc_factors_cv2 <- lapply(factor_result$cv2_factor_structure[-focal_time], function(x){
    procustes(x$factor_structure$loadings,
              target = target_structure$factor_structure$loadings,
              type = 'orthogonal',
              verbose = TRUE)
  })

  smoothed <- lapply(proc_factors_cv2, function(x){
    procustes_smoother(x$rotmatrix)
  })

  # some formatting
  # last timepoint got columns switched
  smoothed_factors <- purrr::map2(factor_result$cv2_factor_structure[-focal_time],
                                  smoothed,
                                  function(timepoint, rotation_matrix){
                                    timepoint$factor_structure$loadings %*% rotation_matrix
                                  })


  origin_loadings <- as.matrix(factor_result$cv2_factor_structure[[focal_time]]$factor_structure$loadings[,])
  colnames(origin_loadings) <- gsub(x = colnames(origin_loadings), pattern = 'F', replacement =  'V')

  smoothed_factors <- append(smoothed_factors, list(origin_loadings), after = focal_time-1)
  names(smoothed_factors) <- timepoints

  # generate factor scores
  proc_factor_scores_cv2 <- vector(mode = 'list', length = length(smoothed_factors))
  for(timepoint_i in 1:length(smoothed_factors)){
    temp_df <-dplyr::select(x$data, G, contains(timepoints[timepoint_i]))

    # factor score
    proc_factor_scores_cv2[[timepoint_i]] <- gfBLUP::factorScores(as.data.frame(temp_df),
                                                                  loadings = smoothed_factors[[timepoint_i]],
                                                                  uniquenesses = factor_result$cv2_factor_structure[[timepoint_i]]$factor_structure$uniquenesses,
                                                                  m = target_structure$factor_structure$m)

    # add plot level data for splines
    proc_factor_scores_cv2[[timepoint_i]] <- x$data %>%
      dplyr::bind_cols(proc_factor_scores_cv2[[timepoint_i]][, -1]) %>%
      dplyr::select(G:plot, F1:F2, yield_adjusted_BLUES:Y)
  }
  names(proc_factor_scores_cv2 ) <- timepoints

  # create result file
  cv2_proc_result <- vector(mode = 'list', length = length(proc_factor_scores_cv2))

  for(timepoint_i in 1:length(proc_factor_scores_cv2)){
    cv2_proc_result[[timepoint_i]]$factor_scores <- proc_factor_scores_cv2[[timepoint_i]]
    cv2_proc_result[[timepoint_i]]$factor_loadings <- smoothed_factors[[timepoint_i]]
    cv2_proc_result[[timepoint_i]]$uniquenesses <- factor_result$cv2_factor_structure[[timepoint_i]]$factor_structure$uniquenesses
  }

  names(cv2_proc_result) <- timepoints

  # cv1 ----
  # generate factor loadings at focal timepoint
  cv1_df <- x$data[which(x$data$G %in% x$train.set), ]
  target_structure <- factor_result$cv1_factor_sctructure[[focal_time]]

  # Procrustean smoothing ----
  # create procustes rotated factors for each timepoint
  proc_factors_cv1 <- lapply(factor_result$cv1_factor_sctructure[-focal_time], function(x){
    procustes(loadings = x$factor_structure$loadings,
              target = target_structure$factor_structure$loadings,
              type = 'orthogonal', verbose = TRUE)
  })

  smoothed <- lapply(proc_factors_cv1, function(x){
    procustes_smoother(x$rotmatrix)
  })

  # some formatting
  # last timepoint got columns switched
  smoothed_factors <- purrr::map2(factor_result$cv1_factor_sctructure[-focal_time],
                                  smoothed, function(factor_object, rotation_matrix){
                                    factor_object$factor_structure$loadings %*% rotation_matrix
                                  })

  origin_loadings <- as.matrix(target_structure$factor_structure$loadings[,])
  colnames(origin_loadings) <- gsub(x = colnames(origin_loadings), pattern = 'F', replacement =  'V')

  smoothed_factors <- append(smoothed_factors, list(origin_loadings), after = 5-1)
  names(smoothed_factors) <- timepoints

  # generate factor scores
  proc_factor_scores_cv1 <- vector(mode = 'list', length = length(smoothed_factors))
  for(timepoint_i in 1:length(smoothed_factors)){
    temp_df <-dplyr::select(cv1_df, G, contains(timepoints[timepoint_i]))

    # factor score
    proc_factor_scores_cv1[[timepoint_i]] <- gfBLUP::factorScores(as.data.frame(temp_df),
                                                                  loadings = smoothed_factors[[timepoint_i]],
                                                                  uniquenesses = factor_result$cv1_factor_sctructure[[timepoint_i]]$factor_structure$uniquenesses,
                                                                  m = target_structure$factor_structure$m)
    # add plot level data for splines
    proc_factor_scores_cv1[[timepoint_i]] <- cv1_df %>%
      dplyr::bind_cols(proc_factor_scores_cv1[[timepoint_i]][, -1]) %>%
      dplyr::select(G:plot, F1:F2, yield_adjusted_BLUES:Y)
  }

  names(proc_factor_scores_cv1) <- timepoints

  # create result file
  cv1_proc_result <- vector(mode = 'list', length = length(proc_factor_scores_cv1))

  for(timepoint_i in 1:length(proc_factor_scores_cv1)){
    cv1_proc_result[[timepoint_i]]$factor_scores <- proc_factor_scores_cv1[[timepoint_i]]
    cv1_proc_result[[timepoint_i]]$factor_loadings <- smoothed_factors[[timepoint_i]]
    cv1_proc_result[[timepoint_i]]$uniquenesses <- factor_result$cv1_factor_sctructure[[timepoint_i]]$factor_structure$uniquenesses
  }

  names(cv1_proc_result) <- timepoints


  # store results
  factor_result_file <- sprintf('./data/processed/factors/%s/timepoint_factors/procrustes/procrustes_%d.RData',
                                selected_env,
                                i
                                )

  procrustes_factor_results <- list(CV1 = cv1_proc_result,
                                    CV2 = cv2_proc_result)

  save(procrustes_factor_results, file = factor_result_file)
}
doParallel::stopImplicitCluster()
parallel::stopCluster(cl)
toc()

# test a single run procrustes rotated ----
# for (file.name in list.files("functions")) {
#   source(paste0("functions/", file.name))
# }
# i <- 1
# # focal times:
# # B5IR: 4
# # DRIP, B2IR: 3
# focal_time <- 4
#
# # load in files
# factor_file <- paste0('./data/processed/factors/B5IR/timepoint_factors/varimax/factors_dataset_',
#                     i,
#                     '.RData')
# load(factor_file)
#
# # load in hyperspectral file
# data_file <- paste0('./data/processed/hyper_folds/B5IR/hyper_dataset_',
#                     i,
#                     '.RData')
# load(data_file)
#
# # select focal timepoint
# # cols 1:4 are G to plot, cols 622:626 are
# timepoints <- names(factor_result$cv2_factor_structure)
# head_timepoint <- timepoints[focal_time]
#
# target_structure <- factor_result$cv2_factor_structure[[focal_time]]
#
# # Procrustean smoothing ----
# # create procustes rotated factors for each timepoint
# proc_factors_cv2 <- lapply(factor_result$cv2_factor_structure[-focal_time], function(x){
#   procustes(x$factor_structure$loadings,
#             target = target_structure$factor_structure$loadings,
#             type = 'orthogonal',
#             verbose = TRUE)
# })
#
# smoothed <- lapply(proc_factors_cv2, function(x){
#   procustes_smoother(x$rotmatrix)
# })
#
# # some formatting
# # last timepoint got columns switched
# smoothed_factors <- purrr::map2(factor_result$cv2_factor_structure[-focal_time],
#                                 smoothed,
#                                 function(timepoint, rotation_matrix){
#                                   timepoint$factor_structure$loadings %*% rotation_matrix
#                                 })
#
#
# origin_loadings <- as.matrix(factor_result$cv2_factor_structure[[focal_time]]$factor_structure$loadings[,])
# colnames(origin_loadings) <- gsub(x = colnames(origin_loadings), pattern = 'F', replacement =  'V')
#
# smoothed_factors <- append(smoothed_factors, list(origin_loadings), after = focal_time-1)
# names(smoothed_factors) <- timepoints
#
# # generate factor scores
# proc_factor_scores_cv2 <- vector(mode = 'list', length = length(smoothed_factors))
# for(timepoint_i in 1:length(smoothed_factors)){
#   temp_df <-dplyr::select(x$data, G, contains(timepoints[timepoint_i]))
#
#   # factor score
#   proc_factor_scores_cv2[[timepoint_i]] <- gfBLUP::factorScores(as.data.frame(temp_df),
#                                                                 loadings = smoothed_factors[[timepoint_i]],
#                                                                 uniquenesses = factor_result$cv2_factor_structure[[timepoint_i]]$factor_structure$uniquenesses,
#                                                                 m = target_structure$factor_structure$m)
#
#   # add plot level data for splines
#   proc_factor_scores_cv2[[timepoint_i]] <- x$data %>%
#     dplyr::bind_cols(proc_factor_scores_cv2[[timepoint_i]][, -1]) %>%
#     dplyr::select(1:4, 627:628, 625:626)
# }
# names(proc_factor_scores_cv2 ) <- timepoints
#
# # create result file
# cv2_proc_result <- vector(mode = 'list', length = length(proc_factor_scores_cv2))
#
# for(timepoint_i in 1:length(proc_factor_scores_cv2)){
#   cv2_proc_result[[timepoint_i]]$factor_scores <- proc_factor_scores_cv2[[timepoint_i]]
#   cv2_proc_result[[timepoint_i]]$factor_loadings <- smoothed_factors[[timepoint_i]]
#   cv2_proc_result[[timepoint_i]]$uniquenesses <- factor_result$cv2_factor_structure[[timepoint_i]]$factor_structure$uniquenesses
# }
#
# names(cv2_proc_result) <- timepoints
#
# # cv1 ----
# # generate factor loadings at focal timepoint
# cv1_df <- x$data[which(x$data$G %in% x$train.set), ]
# target_structure <- factor_result$cv1_factor_sctructure[[focal_time]]
#
# # Procrustean smoothing ----
# # create procustes rotated factors for each timepoint
# proc_factors_cv1 <- lapply(factor_result$cv1_factor_sctructure[-focal_time], function(x){
#   procustes(loadings = x$factor_structure$loadings,
#             target = target_structure$factor_structure$loadings,
#             type = 'orthogonal', verbose = TRUE)
# })
#
# smoothed <- lapply(proc_factors_cv1, function(x){
#   procustes_smoother(x$rotmatrix)
# })
#
# # some formatting
# # last timepoint got columns switched
# smoothed_factors <- purrr::map2(factor_result$cv1_factor_sctructure[-focal_time],
#                                 smoothed, function(factor_object, rotation_matrix){
#                                   factor_object$factor_structure$loadings %*% rotation_matrix
#                                 })
#
# origin_loadings <- as.matrix(target_structure$factor_structure$loadings[,])
# colnames(origin_loadings) <- gsub(x = colnames(origin_loadings), pattern = 'F', replacement =  'V')
#
# smoothed_factors <- append(smoothed_factors, list(origin_loadings), after = 5-1)
# names(smoothed_factors) <- timepoints
#
# # generate factor scores
# proc_factor_scores_cv1 <- vector(mode = 'list', length = length(smoothed_factors))
# for(timepoint_i in 1:length(smoothed_factors)){
#   temp_df <-dplyr::select(cv1_df, G, contains(timepoints[timepoint_i]))
#
#   # factor score
#   proc_factor_scores_cv1[[timepoint_i]] <- gfBLUP::factorScores(as.data.frame(temp_df),
#                                                                 loadings = smoothed_factors[[timepoint_i]],
#                                                                 uniquenesses = factor_result$cv1_factor_sctructure[[timepoint_i]]$factor_structure$uniquenesses,
#                                                                 m = target_structure$factor_structure$m)
#   # add plot level data for splines
#   proc_factor_scores_cv1[[timepoint_i]] <- cv1_df %>%
#     dplyr::bind_cols(proc_factor_scores_cv1[[timepoint_i]][, -1]) %>%
#     dplyr::select(1:4, 627:628, 625:626)
# }
#
# names(proc_factor_scores_cv1) <- timepoints
#
# # create result file
# cv1_proc_result <- vector(mode = 'list', length = length(proc_factor_scores_cv1))
#
# for(timepoint_i in 1:length(proc_factor_scores_cv1)){
#   cv1_proc_result[[timepoint_i]]$factor_scores <- proc_factor_scores_cv1[[timepoint_i]]
#   cv1_proc_result[[timepoint_i]]$factor_loadings <- smoothed_factors[[timepoint_i]]
#   cv1_proc_result[[timepoint_i]]$uniquenesses <- factor_result$cv1_factor_sctructure[[timepoint_i]]$factor_structure$uniquenesses
# }
#
# names(cv1_proc_result) <- timepoints
#
#
# # store results
# factor_result_file <- paste0('./data/processed/factors/B5IR/timepoint_factors/procrustes/procrustes_',
#                              i,
#                              '.RData')
#
# procrustes_factor_results <- list(CV1 = cv1_proc_result,
#                                CV2 = cv2_proc_result)
#
# save(procrustes_factor_results, file = factor_result_file)
#
