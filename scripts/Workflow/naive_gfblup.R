################################################################################
# perform naive gfBLUP as done in Melsen et al. (2025)
# but with Bayesian Information Criterion (BIC)
################################################################################

# ----------------------------------------------------
# 1. Load Required Libraries
# ----------------------------------------------------
library(doParallel)
library(tictoc)
library(dplyr)
library(magrittr)
library(stringr)
library(gfBLUP)
library(rlist)

# ----------------------------------------------------
# 2. Source Functions & Define Hyperparameters
# ----------------------------------------------------
# Load all custom helper functions from the local directory
for (file.name in list.files("functions")) {
  source(paste0("functions/", file.name))
}

# Experimental configurations
environments <- c('B5IR', 'B2IR', 'DRIP')
plant_stages <- c('vegetative', 'heading', 'full') # Growth stage subsets
selection_methods <- c('bic', 'adjr2')            # Criteria for regsubsets factor selection
n_datasets <- 100                                 # Total cross-validation folds/iterations

# ----------------------------------------------------
# 3. Processing & Prediction Pipeline
# ----------------------------------------------------
for(selection_method in selection_methods){
  for(plant_stage in plant_stages){
    print(plant_stage)
    for(selected_env in environments){
      # Load genomic relationship matrix (kinship)
      load('./data/processed/kinship.RData')

      # Define environment-specific heading timepoint thresholds
      # B5IR reaches heading at timepoint 5; B2IR and DRIP reach heading at timepoint 4
      if(selected_env == 'B5IR'){
        heading <- 5
      } else {
        heading <- 4
      }

      for(i in 1:n_datasets){
        # --- Data Ingestion & Timepoint Filtering ---
        # Load hyperspectral fold dataset for current iteration
        datalist <- list.load(file = sprintf('./data/processed/hyper_folds/%s/hyper_dataset_%d.Rdata',
                                             selected_env, i))

        d <- datalist$data
        # Extract hyperspectral band feature column names containing numbers
        hs_cols <- colnames(d)[str_detect(colnames(d), '[0-9]+')]

        # Parse unique timepoints encoded within column names (e.g., "wavelength_timepoint")
        timepoints <- unique(str_split(hs_cols, '_', simplify = TRUE)[, 2])

        # Filter hyperspectral features according to selected growth stage
        if(plant_stage == 'heading'){
          d <- d %>%
            select(G, any_of(contains(timepoints[1:heading])), yield_adjusted_BLUES, Y)
        }
        if(plant_stage == 'vegetative'){
          vegetative <- heading - 1
          d <- d %>%
            select(G, any_of(contains(timepoints[1:vegetative])), yield_adjusted_BLUES, Y)
        }
        pred.target <- datalist$pred.target

        # --- Model Construction & Execution ---

        # 1. Prepare training data and extract feature / focal trait identifiers
        d.train <- droplevels(na.omit(d))
        sec <- names(d[5:(ncol(d) - 2)]) # Secondary traits (hyperspectral features)
        foc <- names(d)[ncol(d)-1]       # Focal target trait (yield_adjusted_BLUES)

        # 2. Apply redundancy filter to secondary features on training set only (threshold tau = 0.95)
        temp <- gfBLUP::redundancyFilter(data = d.train[c("G", sec)], tau = 0.95, verbose = FALSE)
        d.train.RF <- cbind(temp$data.RF, d.train[foc])
        d.RF <- d[names(d.train.RF)]
        sec.RF <- names(d.RF[2:(ncol(d.RF) - 1)])

        # 3. Estimate regularized genetic and residual correlation matrices using cross-validation folds
        folds <- gfBLUP::createFolds(genos = unique(as.character(d.train.RF$G)))
        tempG <- gfBLUP::regularizedCorrelation(data = d.train.RF[c("G", sec.RF)], folds = folds, what = "genetic", dopar = FALSE, verbose = FALSE)
        tempE <- gfBLUP::regularizedCorrelation(data = d.train.RF[c("G", sec.RF)], folds = folds, what = "residual", dopar = FALSE, verbose = FALSE)
        Rg.RF.reg <- tempG$optCor # Optimal regularized genetic correlation matrix

        # 4. Fit latent factor model to regularized genetic correlations
        FM.fit <- gfBLUP::factorModel(data = d.train.RF[c("G", sec.RF)], cormat = Rg.RF.reg, what = "genetic", verbose = FALSE)

        # 5. Calculate Thomson genetic factor scores (rescaling loadings and uniquenesses to covariance scale)
        D <- sqrt(diag(tempG$Sg)) # Standard deviations from genetic covariance matrix
        L.cov <- diag(D) %*% FM.fit$loadings
        PSI.cov <- outer(D, D) * FM.fit$uniquenesses

        # Calculate CV1 Factor Scores (Mask test set yield values with NA)
        CV1.d.RF <- d.RF
        CV1.d.RF[which(is.na(CV1.d.RF$yield_adjusted_BLUES)), 2:ncol(CV1.d.RF)] <- NA
        CV1.F.scores <- gfBLUP::factorScores(data = CV1.d.RF[c("G", sec.RF)],
                                             loadings = L.cov,
                                             uniquenesses = PSI.cov,
                                             m = FM.fit$m,
                                             type = "genetic-thomson-repdiv",
                                             Se = outer(sqrt(diag(tempE$Se)), sqrt(diag(tempE$Se))) * tempE$optCor)

        CV1.d.final <- cbind(CV1.F.scores, CV1.d.RF$yield_adjusted_BLUES)
        names(CV1.d.final)[ncol(CV1.d.final)] <- "Y"
        names(CV1.d.final)[1] <- "G"

        # Calculate CV2 Factor Scores (Standardize secondary features across both train and test sets)
        d.RF[sec.RF] <- sapply(d.RF[sec.RF], scale)
        CV2.F.scores <- gfBLUP::factorScores(data = d.RF[c("G", sec.RF)],
                                             loadings = L.cov,
                                             uniquenesses = PSI.cov,
                                             m = FM.fit$m,
                                             type = "genetic-thomson-repdiv",
                                             Se = outer(sqrt(diag(tempE$Se)), sqrt(diag(tempE$Se))) * tempE$optCor)

        CV2.d.final <- cbind(CV2.F.scores, d.RF$yield_adjusted_BLUES)
        names(CV2.d.final)[ncol(CV2.d.final)] <- "Y"
        names(CV2.d.final)[1] <- "G"

        # 6. Select informative factors using specified criterion (BIC or Adjusted R^2)
        selection <- regsubset_factorSelect(CV1.d.final, procedure = "regsubsets", method = selection_method, verbose = FALSE)

        # 7. Fit multi-trait genomic prediction models (gfBLUP) for CV1 and CV2 schemes
        K <- kinship[rownames(kinship) %in% CV1.d.final$G, rownames(kinship) %in% CV1.d.final$G]
        CV1.temp <- gfBLUP::gfBLUP(data = CV1.d.final, selection = selection, K = K, sepExp = FALSE, verbose = F)
        CV2.temp <- gfBLUP::gfBLUP(data = CV2.d.final, selection = selection, K = K, sepExp = FALSE, verbose = F)

        # --- Predictive Accuracy Calculation & Result Archiving ---
        # Calculate predictive accuracy (Pearson's rho) against ground-truth validation targets
        rho_cv1 <- cor(pred.target$pred.target, CV1.temp$preds[match(pred.target$G, names(CV1.temp$preds))])
        rho_cv2 <- cor(pred.target$pred.target, CV2.temp$preds[match(pred.target$G, names(CV2.temp$preds))])

        # Compile model diagnostics, variance components, and prediction correlations
        result_list <- list(loadings = FM.fit$loadings,
                            uniquenesses = FM.fit$uniquenesses,
                            m = FM.fit$m,
                            m.selected = selection,
                            Sg = CV1.temp$Sg,
                            Se = CV1.temp$Se,
                            h2s = CV1.temp$h2s,
                            rho_cv1 = rho_cv1,
                            rho_cv2 = rho_cv2)

        # Save iteration results to output file structure
        list.save(result_list, sprintf('./data/processed/gfBLUP_naive_regsubset_%s/%s/%s/result_%d.RData',
                                       selection_method, selected_env, plant_stage, i))
      }
    }
  }
}
