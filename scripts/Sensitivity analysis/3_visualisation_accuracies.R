################################################################################
# Visualisation of prediction accuracies of all methods across all
# CV scenario, plant stages and environments
################################################################################

# libraries ----
library(emmeans)
library(multcomp)
library(tidyverse)
library(magrittr)
library(correctR)

# settings ----
# create settings df
settings_df <- expand.grid(
  CV = c('CV1', 'CV2'),
  selected_env = c('B5IR', 'B2IR', 'DRIP'),
  plant_stage = c('full', 'heading', 'vegetative'),
  m = paste('m',1:3, sep = '_'),
  stringsAsFactors = FALSE
)

# palette
type_colours_2 <- c('#009E73',
                    '#CC79A7',
                    '#E69F00'
                    )


# timepoint factors ----
## varimax ----
accuracy_list <- by(settings_df, seq_len(nrow(settings_df)), function(row){
  selected_env <-   row$selected_env
  plant_stage <- row$plant_stage
  m <- row$m
  CV <- row$CV

  # retrieve filenames for specific setting
  file_names <- list.files(sprintf('./data/processed/sensitivity_analysis/%s/glfBLUP/%s/%s/%s/',
                                   selected_env, m, CV, plant_stage), full.names = TRUE)

  # retrieve accuracy from datafile
  accuracies <- lapply(file_names, function(file_name){
    load(file_name)
    return(data.frame(rho = x$rho))
  })

  #wrangle accuracies to a data frame with dataset as identifier
  names(accuracies) <- str_extract(str_extract(file_names, '[0-9]+.RData'), '[0-9]+')
  accuracy_df <- bind_rows(accuracies, .id = 'dataset')
  return(accuracy_df)
}, simplify = FALSE)

# provide setting names to accuracy list
names(accuracy_list) <- paste(settings_df$selected_env,
                                      # because _ will be separator
                                      str_remove(settings_df$m, '_'),
                                      settings_df$CV, settings_df$plant_stage, sep = '_')

# by function added attributes we want to remove
accuracy_list <- lapply(accuracy_list, function(x){x})

# create accuracy df with setting columns
accuracy_df <- bind_rows(accuracy_list, .id = 'setting') |>
  separate_wider_delim(setting, names = c('environment', 'm', 'CV', 'plant_stage'), delim = '_') |>
  mutate(m = str_remove(m, 'm')) |>
  filter(as.double(dataset) <= 100)

## perform significance testing using Nadeau and Bengio (2003) procedure ----
# train-test folds are related over 100 subsamples, so we use Nadeau & Bengio (2003) procedure
# this is implemented in correctR::resampled_ttest
# n = number of subsamples (train-test fold splits).
# n1 = number of train genotypes 1092 - 366 = 722
# n2 = number of test genotypes 311
split_dfs <- split(accuracy_df, ~environment + CV + plant_stage)

m_pairs <- combn(unique(split_dfs[[1]]$m), m = 2)

# generate compact letter display (cld) based on p value groupings
pval_list <- lapply(split_dfs, function(x){
  wide_df <- x |>
    pivot_wider(names_from = m, values_from = rho)

  pair_df <- as.data.frame(t(m_pairs))
  pair_df$pval <- NULL

  for (i in 1:nrow(pair_df)) {
    tmp_out <- resampled_ttest(x = wide_df[[pair_df$V1[i]]],
                               y = wide_df[[pair_df$V2[i]]], n = 100, n1 = 722, n2 = 311, tailed = 'two')
    pair_df$pval[i] <- tmp_out[[2]]
  }

  pvals <- pair_df$pval
  names(pvals) <- paste(pair_df$V1, pair_df$V2, sep = '-')

  # change denominator to relevant number methods/m (here it is 3 for factor dimension 1 to 3)
  cld <- multcompView::multcompLetters(pvals, threshold = 0.05/3)$Letters
  return(cld)
})

# create cld dataframe for plotting
cld_df <- bind_rows(pval_list, .id = 'environment') |>
  pivot_longer(cols = `1`:`3`, names_to = 'm', values_to = 'cld') |>
  separate_wider_delim('environment', delim = '.', names = c('environment', 'CV', 'plant_stage')) |>
  group_by(CV, environment, plant_stage) |>
  # get counts per grouping of comparisons
  add_count(cld) |>
  # if it only contains one letter (a) it will sum to number of factors (3)
  filter(n != 3) |>
  # remove white space around letters and set position for plot
  mutate(cld = str_trim(cld),
         y = .9,
         .keep = 'unused')

## create boxplot with cld all data ----
complete_signif_boxplot <- accuracy_df |>
  left_join(cld_df) |>
  mutate(plant_stage = fct_recode(factor(plant_stage, levels = c('vegetative', 'heading', 'full')),
                                  Vegetative = 'vegetative',
                                  `Vegetative & Heading` = 'heading',
                                  Full = 'full'),
         environment = fct_recode(factor(environment, levels = c('B5IR', 'B2IR', 'DRIP')),
                                  'Optimal flat' = 'B5IR',
                                  'Moderate drought' = 'B2IR',
                                  'Severe drought' = 'DRIP')) |>
  ggplot(aes(x = plant_stage, y = rho, fill = m)) +
  geom_boxplot() +
  geom_text(aes(y = y, label = cld, colour = m), , size = 3, position = position_dodge(width = .85), stat = 'unique') +
  facet_wrap(~CV + environment) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(fill=guide_legend(title= "Number of factors")) +
  scale_colour_manual(values = type_colours_2, guide = 'none') +
  scale_fill_manual(values = type_colours_2) +
  labs(
    x = 'Plant stages included',
    y = 'Predictive ability')

complete_signif_boxplot

ggsave(filename = './figures/boxplot_sensitivity_analysis.pdf', plot = complete_signif_boxplot, width = 12, height = 10, device = 'pdf')



