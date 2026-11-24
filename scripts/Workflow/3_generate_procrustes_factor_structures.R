################################################################################
# Align latent factor structures across timepoints using Procrustes rotation
################################################################################

# ------------------------------------------------------------------------------
# 1. Package Loading & Global Setup
# ------------------------------------------------------------------------------
library(stringr)
library(tictoc)
library(glfBLUP)
library(doParallel)

# Source custom user functions (e.g., `procustes`, `procustes_smoother`)
for (file.name in list.files("functions")) {
  source(paste0("functions/", file.name))
}

# ------------------------------------------------------------------------------
# 2. Configuration Parameters
# ------------------------------------------------------------------------------
selected_env <- 'B2IR' # Target environment identifier

# Define focal timepoint index used as the reference target matrix for Procrustes rotation
# Environment specific selection: B5IR uses index 4; DRIP and B2IR use index 3
focal_time <- 3

# Count total fold datasets generated for the target environment
n_datasets <- length(list.files(sprintf('./data/processed/factors/%s/timepoint_factors/varimax/', selected_env)))

# Setup parallel processing cluster, leaving 2 CPU cores free for background tasks
n_cores <- parallel::detectCores() - 2
cl <- parallel::makeCluster(n_cores, outfile = "logs/varimax_factor_loadings_generation.txt")
doParallel::registerDoParallel(cl)

# ------------------------------------------------------------------------------
# 3. Parallel Processing Loop (Procrustes Rotation & Factor Scoring)
# ------------------------------------------------------------------------------
tic()
foreach::foreach(i = 1:n_datasets, .packages = c("stringr", 'magrittr', 'dplyr')) %dopar% {

  # ----------------------------------------------------------------------------
  # Data Loading
  # ----------------------------------------------------------------------------
  # Load previously generated Varimax factor loading results
  factor_file <- sprintf('./data/processed/factors/%s/timepoint_factors/varimax/factors_dataset_%d.RData',
                         selected_env,
                         i)
  load(factor_file)

  # Load raw hyperspectral fold partition data
  data_file <- sprintf('./data/processed/hyper_folds/%s/hyper_dataset_%d.RData',
                       selected_env,
                       i)
  load(data_file)

  # Extract timepoint names and isolate target factor loading structure
  timepoints <- names(factor_result$cv2_factor_structure)
  head_timepoint <- timepoints[focal_time]

  # Reference loading matrix to which all other timepoints will be aligned
  target_structure <- factor_result$cv2_factor_structure[[focal_time]]

  # ----------------------------------------------------------------------------
  # CV2 Scenario: Procrustean Rotation Across ALL Genotypes
  # ----------------------------------------------------------------------------
  # Perform orthogonal Procrustes alignment for all non-focal timepoints
  proc_factors_cv2 <- lapply(factor_result$cv2_factor_structure[-focal_time], function(x){
    procustes(x$factor_structure$loadings,
              target = target_structure$factor_structure$loadings,
              type = 'orthogonal',
              verbose = TRUE)
  })

  # Apply smoothing to extracted rotation matrices across time
  smoothed <- lapply(proc_factors_cv2, function(x){
    procustes_smoother(x$rotmatrix)
  })

  # Rotate original factor loadings by multiplying with smoothed rotation matrices
  smoothed_factors <- purrr::map2(factor_result$cv2_factor_structure[-focal_time],
                                  smoothed,
                                  function(timepoint, rotation_matrix){
                                    timepoint$factor_structure$loadings %*% rotation_matrix
                                  })

  # Format target (focal) loadings matrix and rename columns from 'F' to 'V' for appending
  original_loadings <- as.matrix(factor_result$cv2_factor_structure[[focal_time]]$factor_structure$loadings[,])
  colnames(original_loadings) <- gsub(x = colnames(original_loadings), pattern = 'F', replacement =  'V')

  # Re-insert target loading matrix into the timepoint list at its original index position
  smoothed_factors <- append(smoothed_factors, list(original_loadings), after = focal_time-1)
  names(smoothed_factors) <- timepoints

  # Compute factor scores using glfBLUP
  proc_factor_scores_cv2 <- vector(mode = 'list', length = length(smoothed_factors))
  for(timepoint_i in 1:length(smoothed_factors)){
    # Select hyperspectral features corresponding to current timepoint
    temp_df <- dplyr::select(x$data, G, contains(timepoints[timepoint_i]))

    # Calculate factor scores using updated loadings and uniquenesses
    proc_factor_scores_cv2[[timepoint_i]] <- glfBLUP::factorScores(as.data.frame(temp_df),
                                                                  loadings = smoothed_factors[[timepoint_i]],
                                                                  uniquenesses = factor_result$cv2_factor_structure[[timepoint_i]]$factor_structure$uniquenesses,
                                                                  m = target_structure$factor_structure$m)

    # Bind plotdata and yield to calculated factor scores
    proc_factor_scores_cv2[[timepoint_i]] <- x$data %>%
      dplyr::bind_cols(proc_factor_scores_cv2[[timepoint_i]][, -1]) %>%
      dplyr::select(G:plot, F1:F2, yield_adjusted_BLUES:Y)
  }
  names(proc_factor_scores_cv2) <- timepoints

  # Compile CV2 list results
  cv2_proc_result <- vector(mode = 'list', length = length(proc_factor_scores_cv2))

  for(timepoint_i in 1:length(proc_factor_scores_cv2)){
    cv2_proc_result[[timepoint_i]]$factor_scores <- proc_factor_scores_cv2[[timepoint_i]]
    cv2_proc_result[[timepoint_i]]$factor_loadings <- smoothed_factors[[timepoint_i]]
    cv2_proc_result[[timepoint_i]]$uniquenesses <- factor_result$cv2_factor_structure[[timepoint_i]]$factor_structure$uniquenesses
  }

  names(cv2_proc_result) <- timepoints

  # ----------------------------------------------------------------------------
  # CV1 Scenario: Procrustean Rotation Strictly on TRAINING Genotypes
  # ----------------------------------------------------------------------------
  # Filter data to training set entries
  cv1_df <- x$data[which(x$data$G %in% x$train.set), ]

  # Reference loading matrix to which all other timepoints will be aligned
  target_structure <- factor_result$cv1_factor_sctructure[[focal_time]]

  # Align non-focal training factor loadings to target structure
  proc_factors_cv1 <- lapply(factor_result$cv1_factor_sctructure[-focal_time], function(x){
    procustes(loadings = x$factor_structure$loadings,
              target = target_structure$factor_structure$loadings,
              type = 'orthogonal', verbose = TRUE)
  })

  # Apply smoothing to rotation matrices
  smoothed <- lapply(proc_factors_cv1, function(x){
    procustes_smoother(x$rotmatrix)
  })

  # Apply rotated transformation matrices
  smoothed_factors <- purrr::map2(factor_result$cv1_factor_sctructure[-focal_time],
                                  smoothed, function(factor_object, rotation_matrix){
                                    factor_object$factor_structure$loadings %*% rotation_matrix
                                  })

  origin_loadings <- as.matrix(target_structure$factor_structure$loadings[,])
  colnames(origin_loadings) <- gsub(x = colnames(origin_loadings), pattern = 'F', replacement =  'V')

  # Re-insert target loading matrix into the timepoint list at its original index position
  smoothed_factors <- append(smoothed_factors, list(origin_loadings), after = focal_time-1)
  names(smoothed_factors) <- timepoints

  # Compute factor scores for CV1 training subset]
  proc_factor_scores_cv1 <- vector(mode = 'list', length = length(smoothed_factors))
  for(timepoint_i in 1:length(smoothed_factors)){
    temp_df <- dplyr::select(cv1_df, G, contains(timepoints[timepoint_i]))

    proc_factor_scores_cv1[[timepoint_i]] <- glfBLUP::factorScores(as.data.frame(temp_df),
                                                                  loadings = smoothed_factors[[timepoint_i]],
                                                                  uniquenesses = factor_result$cv1_factor_sctructure[[timepoint_i]]$factor_structure$uniquenesses,
                                                                  m = target_structure$factor_structure$m)
    # Bind spatial design metadata and target traits
    proc_factor_scores_cv1[[timepoint_i]] <- cv1_df %>%
      dplyr::bind_cols(proc_factor_scores_cv1[[timepoint_i]][, -1]) %>%
      dplyr::select(G:plot, F1:F2, yield_adjusted_BLUES:Y)
  }

  names(proc_factor_scores_cv1) <- timepoints

  # Compile CV1 list results[cite: 3]
  cv1_proc_result <- vector(mode = 'list', length = length(proc_factor_scores_cv1))

  for(timepoint_i in 1:length(proc_factor_scores_cv1)){
    cv1_proc_result[[timepoint_i]]$factor_scores <- proc_factor_scores_cv1[[timepoint_i]]
    cv1_proc_result[[timepoint_i]]$factor_loadings <- smoothed_factors[[timepoint_i]]
    cv1_proc_result[[timepoint_i]]$uniquenesses <- factor_result$cv1_factor_sctructure[[timepoint_i]]$factor_structure$uniquenesses
  }

  names(cv1_proc_result) <- timepoints

  # ----------------------------------------------------------------------------
  # Export Results
  # ----------------------------------------------------------------------------
  factor_result_file <- sprintf('./data/processed/factors/%s/timepoint_factors/procrustes/procrustes_%d.RData',
                                selected_env,
                                i
  )

  procrustes_factor_results <- list(CV1 = cv1_proc_result,
                                    CV2 = cv2_proc_result)

  save(procrustes_factor_results, file = factor_result_file)
}

# Terminate parallel cluster backend
doParallel::stopImplicitCluster()
parallel::stopCluster(cl)
toc()
