# Script to generate latent factors (via Varimax rotation) per dataset fold

# ------------------------------------------------------------------------------
# 1. Package Loading & Global Setup
# ------------------------------------------------------------------------------
library(stringr)
library(tictoc)
library(gfBLUP)
library(doParallel)

# Load user-defined custom utility functions from local directory
for (file.name in list.files("functions")) {
  source(paste0("functions/", file.name))
}

# Target environment identifier
selected_env <- 'DRIP'

# Count partition dataset files to process
n_datasets <- length(list.files(sprintf('./data/processed/hyper_folds/%s/', selected_env)))

# Allocate CPU cores for parallel execution, leaving 2 cores for OS tasks
n_cores <- parallel::detectCores() - 2

# Spawn parallel socket cluster and direct logs to output file
cl <- parallel::makeCluster(n_cores, outfile = "logs/varimax_factor_loadings_generation.txt")
doParallel::registerDoParallel(cl)

# ------------------------------------------------------------------------------
# 2. Parallel Factor Generation Loop
# ------------------------------------------------------------------------------
tic()
foreach::foreach(i = 1:n_datasets, .packages = c("dplyr", "stringr")) %dopar% {
  # Build path to input fold partition data
  data_file <- sprintf('./data/processed/hyper_folds/%s/hyper_dataset_%d.RData', selected_env, i)

  # Load fold data structure into workspace
  load(data_file)

  # Get genotype, plot identifiers and target yield columns to exclude from factor analysis
  non_hs_cols <- c('G', 'row', 'column', 'plot', 'yield_adjusted_BLUES', 'Y')

  # Boolean vector identifying hyperspectral predictor columns
  hs_cols <- !colnames(x$data) %in% non_hs_cols

  # Extract unique sampling timepoints from hyperspectral column headers (e.g., 'wavelength_timepoint')
  timepoints <- unique(str_split(colnames(x$data)[hs_cols], '_', simplify = TRUE)[, 2])

  # generating 2 factors for timepoint index 4
  generate_factor(df = x$data[, c('G', colnames(x$data)[hs_cols])], timepoint = timepoints[4], m = 2)

  # ----------------------------------------------------------------------------
  # CV2 Scenario: Factor analysis on ALL genotypes (including test entries)
  # ----------------------------------------------------------------------------
  # Extract m = 2 latent factors across every timepoint for the full dataset
  factor_result_cv2 <- lapply(timepoints, generate_factor, df = x$data[, c('G', colnames(x$data)[hs_cols])], m = 2)
  names(factor_result_cv2) <- timepoints

  # Re-attach genotype, plot identifiers and yield to generated factor scores
  factor_dfs_cv2 <- lapply(factor_result_cv2, function(result){
    # Combine plot identifiers metadata (columns 1 to 4) with calculated factor score columns
    temp_df <- dplyr::bind_cols(x$data[,c(1:4)],result$factor_scores[, -1])
    full_df <- dplyr::left_join(x$data[non_hs_cols], temp_df, by = c('G', 'row', 'column', 'plot'))
  })

  # ----------------------------------------------------------------------------
  # CV1 Scenario: Factor analysis strictly on TRAINING genotypes (prevents leakage)[cite: 2]
  # ----------------------------------------------------------------------------
  # Subset hyperspectral observations to training set entries only
  df_cv1 <- x$data[which(x$data$G %in% x$train.set), c('G', colnames(x$data)[hs_cols])]

  # Extract m = 2 latent factors strictly on training subset
  factor_result_cv1 <- lapply(timepoints, generate_factor, df = df_cv1, m =2)
  names(factor_result_cv1) <- timepoints

  # Re-attach spatial design and yield metrics for CV1 partitions
  factor_dfs_cv1 <- lapply(factor_result_cv1, function(result){
    # Combine training spatial metadata with calculated training factor scores
    temp_df <- dplyr::bind_cols(x$data[which(x$data$G %in% x$train.set), c(1:4)],result$factor_scores[, -1])

    # Join metadata back to training factor scores
    full_df <- dplyr::left_join(x$data[non_hs_cols], temp_df, by = c('G', 'row', 'column', 'plot'))

    return(full_df)
  })

  # Compile structure and dataframes for CV1 and CV2 cross-validation schemes
  factor_result <- list(cv1_factor_sctructure = factor_result_cv1,
                        cv1_factor_dfs = factor_dfs_cv1,
                        cv2_factor_structure = factor_result_cv2,
                        cv2_factor_dfs = factor_dfs_cv2)

  # Output target file path
  factor_loadings_file <- sprintf('./data/processed/factors/%s/timepoint_factors/varimax/factors_dataset_%d.RData',
                                  selected_env, i)

  # Export latent factor results to disk
  save(factor_result, file = factor_loadings_file)
}

# Terminate parallel cluster execution
doParallel::stopImplicitCluster()
parallel::stopCluster(cl)
toc()

