# visualisation AUC correlation optimal flat, moderate drought & heat

# libraries ----
library(corrplot)
library(tidyverse)
library(statgenHTP)

# load itidyverse# load in data ----
moderate_factor1 <- read_rds('./data/processed/splines/all_data/predsplinedatall_factor1_mid_timepoint_4_2014_Bed2IR.RDS')
moderate_factor2 <- read_rds('./data/processed/splines/all_data/predsplinedatall_factor2_mid_timepoint_4_2014_Bed2IR.RDS')
severe_factor1 <- read_rds('./data/processed/splines/all_data/predsplinedatall_factor1_mid_timepoint_4_2014_Drip.RDS')
severe_factor2 <- read_rds('./data/processed/splines/all_data/predsplinedatall_factor2_mid_timepoint_4_2014_Drip.RDS')
optimal_factor1 <- read_rds('./data/processed/splines/all_data/predsplinedatall_factor1_mid_timepoint_4_2014_Bed5IR.RDS')
optimal_factor2 <- read_rds('./data/processed/splines/all_data/predsplinedatall_factor2_mid_timepoint_4_2014_Bed5IR.RDS')

AUC_moderate_factor1 <- estimateSplineParameters(moderate_factor1,
                                                 estimate = 'predictions',
                                                 what = 'AUC',
                                                 fitLevel = 'geno')
AUC_moderate_factor2 <- estimateSplineParameters(moderate_factor2,
                                                 estimate = 'predictions',
                                                 what = 'AUC',
                                                 fitLevel = 'geno')
AUC_severe_factor1 <- estimateSplineParameters(severe_factor1,
                                                 estimate = 'predictions',
                                                 what = 'AUC',
                                                 fitLevel = 'geno')
AUC_severe_factor2 <- estimateSplineParameters(severe_factor2,
                                                 estimate = 'predictions',
                                                 what = 'AUC',
                                                 fitLevel = 'geno')
AUC_optimal_factor1 <- estimateSplineParameters(optimal_factor1,
                                               estimate = 'predictions',
                                               what = 'AUC',
                                               fitLevel = 'geno')
AUC_optimal_factor2 <- estimateSplineParameters(optimal_factor2,
                                               estimate = 'predictions',
                                               what = 'AUC',
                                               fitLevel = 'geno')

auc_df <- data.frame(`Moderate Drought Factor 1` = AUC_moderate_factor1$AUC_predictions,
                     `Moderate Drought Factor 2` = AUC_moderate_factor2$AUC_predictions,
                     `Severe Drought Factor 1` = AUC_severe_factor1$AUC_predictions,
                     `Severe Drought Factor 2` = AUC_severe_factor2$AUC_predictions,
                     `Optimal Flat Factor 1` = AUC_optimal_factor1$AUC_predictions,
                     `Optimal Flat Factor 2` = AUC_optimal_factor2$AUC_predictions)

(plot_cor <- corrplot(cor(auc_df),
         method = 'ellipse',
         type = 'lower',
         diag = FALSE,
         order= 'AOE',
         addCoef.col = 'black'))



