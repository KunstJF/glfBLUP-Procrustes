################################################################################
# Generate splines and spline characteristics of Procrustes rotated
# factor scores
################################################################################

# ----------------------------------------------------
# 1. Load Required Libraries
# ----------------------------------------------------
library(statgenHTP)
library(tictoc)
library(tidyverse)
library(magrittr)

# ----------------------------------------------------
# 2. Experimental Setup & Settings Grid
# ----------------------------------------------------
# Generate all combinations of Cross-Validation folds (CV), environments, and plant stages
settings_df <- expand.grid(
  CV = c('CV1', 'CV2'),
  selected_env = c('B5IR', 'B2IR', 'DRIP'),
  plant_stage = c('full', 'heading', 'vegetative'),
  stringsAsFactors = FALSE
)

# ----------------------------------------------------
# 3. Execution Loop
# ----------------------------------------------------
# Sequential loop over each experimental setting combination
tic() # Start timer

for(setting_i in 1:nrow(settings_df)){
  # Extract specific scenario variables for the current iteration
  CV <- settings_df[setting_i,]$CV
  selected_env <- settings_df[setting_i,]$selected_env
  plant_stage <- settings_df[setting_i,]$plant_stage

  # Run iterations
  for(i in 1:100){
    cat(sprintf('start run %d with settings: %s %s %s',
                i, CV, selected_env, plant_stage))

    # --- Data Loading ---
    # Construct input file path for the current Procrustes factors
    factor_file <- sprintf('./data/processed/factors/%s/timepoint_factors/procrustes/procrustes_%d.RData',
                           selected_env,
                           i)
    load(factor_file)

    # --- Data Wrangling to statgenHTP Format ---
    # Extract factor scores from loaded data structure across time points
    factor_dfs <- lapply(procrustes_factor_results[[CV]], function(structure){structure$factor_scores})
    factor_proc_score_df <- bind_rows(factor_dfs, .id = 'timePoint')

    # Reformat columns to match 'statgenHTP' standard requirements
    factor_df <- factor_proc_score_df %>%
      mutate(
        G = as.factor(G),
        timePoint = as.POSIXct(ymd(timePoint)),
        pop = 'population1',
        colId = column,
        rowId = row,
        plotId = plot,
        .keep = 'unused'
      )

    # --- Plant Stage Filtering Rules ---
    # Timepoint 5 corresponds to heading in B5IR; timepoint 4 in DRIP and B2IR
    if(plant_stage == 'vegetative' & selected_env == 'B5IR'){
      factor_df <- factor_df %>% filter(timePoint < unique(factor_df$timePoint)[5])
    }
    if(plant_stage == 'heading' & selected_env == 'B5IR'){
      factor_df <- factor_df %>% filter(timePoint <= unique(factor_df$timePoint)[5])
    }
    if(plant_stage == 'vegetative' & selected_env != 'B5IR'){
      factor_df <- factor_df %>% filter(timePoint < unique(factor_df$timePoint)[4])
    }
    if(plant_stage == 'heading' & selected_env != 'B5IR'){
      factor_df <- factor_df %>% filter(timePoint <= unique(factor_df$timePoint)[4])
    }

    # --- Spline Fitting ---
    # Fit hierarchical Generalized Additive Models (PS-splines) for factors F1 and F2
    splines <- list(F1 = NA, F2 = NA)

    for(factor_col in names(splines)){
      splines[[factor_col]] <- fitSplineHDM(
        inDat = factor_df,
        trait = factor_col,
        genotype = 'G',
        pop = 'pop'
      )
    }

    # --- Interpolation across Time ---
    # Predict spline trajectories over evenly spaced timepoints
    splines <- lapply(splines, function(factor_spline){
      pred_all.psHDM <- predict(
        object = factor_spline,
        newtimes = seq(min(factor_spline$time[["timeNumber"]]),
                       max(factor_spline$time[["timeNumber"]]),
                       length.out = 100),
        pred = list(pop = TRUE, geno = TRUE, plot = TRUE),
        se = list(pop = FALSE, geno = FALSE, plot = FALSE),
        trace = FALSE
      )
    })

    # --- Extract Spline Characteristics ---
    # Derive physiological parameters (min, max, area under curve) from spline curves and derivatives
    extracted_params <- lapply(splines, function(spline){
      params <- list()
      # Fitted predictions (plot level)
      params[['fgeno']] <- estimateSplineParameters(
        spline,
        estimate = 'predictions',
        what = c('min', 'max', 'AUC'),
        fitLevel = 'plot'
      )
      # First derivatives (growth rates)
      params[['fderiv1']] <- estimateSplineParameters(
        spline,
        estimate = 'derivatives',
        what = c('min', 'max', 'AUC'),
        fitLevel = 'plot'
      )
      # Second derivatives (growth acceleration)
      params[['fderiv2']] <- estimateSplineParameters(
        spline,
        estimate = 'derivatives2',
        what = c('min', 'max', 'AUC'),
        fitLevel = 'plot'
      )
      return(params)
    })

    # --- Parameter Dataset Formatting ---
    # Merge predictions, 1st derivatives, and 2nd derivatives per factor
    extracted_dfs_all <- lapply(extracted_params, function(params){
      reduce(params, left_join, by = c('pop', 'genotype', 'plotId'))
    })

    # Pivot into a single wide-format dataset combining F1 and F2 parameters
    all_df <- bind_rows(extracted_dfs_all, .id = 'factor') %>%
      pivot_wider(names_from = 'factor', values_from = min_predictions:AUC_derivatives2)

    # --- Save Outputs ---
    # Save fitted spline models as RData
    spline_result_file <- sprintf('./data/processed/splines/%s/%s/%s/spline_%d.RData',
                                  selected_env, CV, plant_stage, i)
    save(splines, file = spline_result_file)

    # Save extracted characteristics as CSV
    spline_result_file <- sprintf('./data/processed/spline_characteristics/%s/%s/%s/characteristic_%d.csv',
                                  selected_env, CV, plant_stage, i)
    write_csv(all_df, spline_result_file)

    # --- Clean Memory ---
    # Remove large objects to ensure memory efficiency for subsequent iterations
    rm(all_df, splines, extracted_dfs_all, procrustes_factor_results,
       extracted_params, factor_df, factor_proc_score_df, factor_dfs)
  }
}
toc() # Print execution time
