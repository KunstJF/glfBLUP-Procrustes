################################################################################
# Perform univariate gBLUP for yield across environments
################################################################################

# ----------------------------------------------------
# 1. Load Required Libraries
# ----------------------------------------------------
library(doParallel)
library(tictoc)
library(dplyr)
library(magrittr)
library(rrBLUP)

# ----------------------------------------------------
# 2. Source Helper Functions & Cluster Setup
# ----------------------------------------------------
# Load all custom utility functions from local directory
for (file.name in list.files("functions")) {
  source(paste0("functions/", file.name))
}

# Configure multi-core parallel processing cluster (reserving 2 cores for OS)
n.cores <- parallel::detectCores() - 2
cl <- parallel::makeCluster(n.cores, outfile = "logs/hyper_data_generation.txt")
doParallel::registerDoParallel(cl)

# ----------------------------------------------------
# 3. Environment & Kinship Initialization
# ----------------------------------------------------
environments <- c('B5IR', 'B2IR', 'DRIP')
load('./data/processed/kinship.RData') # Load global genomic relationship matrix

n_datasets <- 100
selected_env = environments[1] # Select target environment (e.g., 'B5IR')

# ----------------------------------------------------
# 4. Parallel gBLUP Execution Loop
# ----------------------------------------------------
tic() # Start performance timer

invisible(
  foreach::foreach(
    i = 1:n_datasets,
    .packages = c("magrittr", 'dplyr', 'rrBLUP', 'rlist')
  ) %dopar% {

    # --- Load Cross-Validation Fold Data ---
    load(sprintf('./data/processed/hyper_folds/%s/hyper_dataset_%d.Rdata', selected_env, i))

    # Extract distinct genotype IDs (G) and observed phenotype values (Y)
    df <- x$data %>%
      select(G, Y) %>%
      distinct()

    # Subset kinship matrix to match genotypes present in training dataset
    kinship <- kinship[rownames(kinship) %in% df$G, colnames(kinship) %in% df$G]

    # --- Fit Genomic BLUP Model ---
    # Estimate Genomic Estimated Breeding Values (GEBVs) via mixed model solver
    gblup_result <- kin.blup(df, geno = 'G', pheno = 'Y', K = kinship)

    # Format predictions into a clean data frame
    gblup_df <- data.frame(
      G = names(gblup_result$pred),
      Y_pred = as.numeric(gblup_result$pred)
    )

    # Optional visual diagnostic of predicted vs actual validation targets
    x$pred.target %>%
      left_join(gblup_df, by = 'G') %$%
      plot(pred.target, Y_pred)

    # --- Calculate Predictive Accuracy ---
    # Compute Pearson correlation (rho) between predictions and true target values
    gblup_result$rho <- x$pred.target %>%
      left_join(gblup_df, by = 'G') %$%
      cor(pred.target, Y_pred)

    # --- Save Iteration Results ---
    list.save(gblup_result, file = sprintf('./data/processed/gblup/%s/gblup_%d.RData', selected_env, i))
  }
)

toc() # Stop timer and report execution time
