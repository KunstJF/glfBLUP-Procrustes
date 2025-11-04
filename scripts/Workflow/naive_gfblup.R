# naive gfBLUP

# libraries ----
library(doParallel)
library(tictoc)
library(dplyr)
library(magrittr)
library(stringr)
library(gfBLUP)
library(rlist)

# Loading functions:
for (file.name in list.files("functions")) {
  source(paste0("functions/", file.name))
}

environments <- c('B5IR', 'B2IR', 'DRIP')
plant_stage <-  'vegetative'
n_datasets <- 100
for(selected_env in environments){
  load('./data/processed/kinship.RData')

  # plant stage
  # B5IR heading happens at timepoint 5
  # B2IR and DRIP heading happens at timepoint 4
  if(selected_env == 'B5IR'){
    heading <- 5
  } else {
    heading <- 4
  }
  for(i in 1:n_datasets){
    # Loading hyperspectral dataset:
    datalist <- list.load(file = sprintf('./data/processed/hyper_folds/%s/hyper_dataset_%d.Rdata',
                                         selected_env, i))

    # Storing data and prediction target:
    d <- datalist$data
    hs_cols <- colnames(d)[str_detect(colnames(d), '[0-9]+')]

    timepoints <- unique(str_split(hs_cols, '_', simplify = TRUE)[, 2])

    if(plant_stage == 'heading'){
      d <- d %>%
        select(G,any_of(contains(timepoints[1:heading])), yield_adjusted_BLUES, Y)
    }
    if(plant_stage == 'vegetative'){
      vegetative <- heading - 1
      d <- d %>%
        select(G,any_of(contains(timepoints[1:vegetative])), yield_adjusted_BLUES, Y)
    }
    pred.target <- datalist$pred.target

    ### Model ##############################################################

    ### 1. Make training data and store feature/focal trait names ------------------------------------------------------------------------
    d.train <- droplevels(na.omit(d))
    sec <- names(d[5:(ncol(d) - 2)])
    foc <- names(d)[ncol(d)-1]

    ### 2. Redundancy filter the secondary features using training data only -------------------------------------------------------------
    temp <- gfBLUP::redundancyFilter(data = d.train[c("G", sec)], tau = 0.95, verbose = FALSE)
    d.train.RF <- cbind(temp$data.RF, d.train[foc])
    d.RF <- d[names(d.train.RF)]
    sec.RF <- names(d.RF[2:(ncol(d.RF) - 1)])

    ### 3. Regularization ----------------------------------------------------------------------------------------------------------------
    folds <- gfBLUP::createFolds(genos = unique(as.character(d.train.RF$G)))
    tempG <- gfBLUP::regularizedCorrelation(data = d.train.RF[c("G", sec.RF)], folds = folds, what = "genetic", dopar = FALSE, verbose = FALSE)
    tempE <- gfBLUP::regularizedCorrelation(data = d.train.RF[c("G", sec.RF)], folds = folds, what = "residual", dopar = FALSE, verbose = FALSE)
    Rg.RF.reg <- tempG$optCor

    ### 4. Fitting factor model ----------------------------------------------------------------------------------------------------------
    # data is only used to determine the sample size for the MP-bound. what specifies that it's a genetic correlation matrix, so the
    # number of training genotypes should be used, and not the number of training individuals (= genotypes * replicates).
    FM.fit <- gfBLUP::factorModel(data = d.train.RF[c("G", sec.RF)], cormat = Rg.RF.reg, what = "genetic", verbose = FALSE)

    #### 5. Getting factor scores (also for the test set in CV2!) ------------------------------------------------------------------------
    # Loadings and uniquenesses were estimated on the correlation scale, but should be on the covariance scale for genetic-thomson scores:
    D <- sqrt(diag(tempG$Sg)) # Getting standard deviations
    L.cov <- diag(D) %*% FM.fit$loadings
    PSI.cov <- outer(D, D) * FM.fit$uniquenesses

    # CV1 Factor scores:
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

    # CV2 Factor scores:
    # First recenter/rescale the training and test secondary data together:
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

    #### 6. Selecting the relevant factors -----------------------------------------------------------------------------------------------
    selection <- gfBLUP::factorSelect(CV1.d.final, procedure = "leaps", verbose = FALSE)

    #### 7. Multi-trait genomic prediction -----------------------------------------------------------------------------------------------
    K <- kinship[rownames(kinship) %in% CV1.d.final$G, rownames(kinship) %in% CV1.d.final$G]
    CV1.temp <- gfBLUP::gfBLUP(data = CV1.d.final, selection = selection, K = K, sepExp = FALSE, verbose = F)
    CV2.temp <- gfBLUP::gfBLUP(data = CV2.d.final, selection = selection, K = K, sepExp = FALSE, verbose = F)

    ########################################################################
    rho_cv1 <- cor(pred.target$pred.target, CV1.temp$preds[match(pred.target$G, names(CV1.temp$preds))])
    rho_cv2 <- cor(pred.target$pred.target, CV2.temp$preds[match(pred.target$G, names(CV2.temp$preds))])

    result_list <- list(loadings = FM.fit$loadings,
                        uniquenesses = FM.fit$uniquenesses,
                        m = FM.fit$m,
                        m.selected = selection,
                        Sg = CV1.temp$Sg,
                        Se = CV1.temp$Se,
                        h2s = CV1.temp$h2s,
                        rho_cv1 = rho_cv1,
                        rho_cv2 = rho_cv2)

    list.save(result_list, sprintf('./data/processed/gfBLUP_naive/%s/%s/result_%d.RData',
                                   selected_env, plant_stage, i))
  }
}

