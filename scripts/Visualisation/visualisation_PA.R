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
settings_df <- expand.grid(CV = c('CV1', 'CV2'),
                           selected_env = c('B5IR', 'B2IR', 'DRIP'),
                           plant_stage = c('full', 'heading', 'vegetative'),
                           stringsAsFactors = FALSE)

# palette
type_colours <- c('#E18727',
                    '#BC3C29',
                    '#0072B5',
                    '#20854E')

# glfBLUP with splines ----
spline_pa_list <- by(settings_df, seq_len(nrow(settings_df)), function(row){
  selected_env <-   row$selected_env
  CV <- row$CV
  plant_stage <- row$plant_stage

  # retrieve filenames for specific setting
  file_names <- list.files(sprintf('./data/processed/glfBLUP_spline/%s/%s/%s/',
                                   selected_env, CV, plant_stage))

  # retrieve pa from datafile
  accuracies <- lapply(file_names, function(file_name){
    file_name_full <- sprintf('./data/processed/glfBLUP_spline/%s/%s/%s/%s',
                              selected_env, CV, plant_stage, file_name)

    load(file_name_full)
    return(data.frame(rho = x$rho))
  })

  #wrangle accuracies to a data frame with dataset as identifier
  names(accuracies) <- str_extract(file_names, '[0-9]+')
  pa_df <- bind_rows(accuracies, .id = 'dataset')
  return(pa_df)
}, simplify = FALSE)

# provide setting names to pa list
names(spline_pa_list) <- paste(settings_df$selected_env, settings_df$CV, settings_df$plant_stage, sep = '_')

# by function added attributes we want to remove
spline_pa_list <- lapply(spline_pa_list, function(x){x})

# create PA df with setting columns
spline_pa_df <- bind_rows(spline_pa_list, .id = 'setting') %>%
  separate_wider_delim(setting, names = c('environment', 'CV', 'plant_stage'), delim = '_')

# naive glfBLUP ----
settings_naive <- settings_df %>%
  select(selected_env, plant_stage) %>%
  distinct()

selected_env <- settings_naive$selected_env[1]
plant_stage <- settings_naive$plant_stage[1]
naive_pa_list <- by(settings_naive, seq_len(nrow(settings_naive)), function(row){
  selected_env <-   row$selected_env
  plant_stage <- row$plant_stage

  # retrieve filenames for specific setting
  file_names <- list.files(sprintf('./data/processed/glfBLUP_naive_regsubset_bic/%s/%s/',
                                   selected_env, plant_stage))

  # retrieve pa from datafile
  accuracies <- lapply(file_names, function(file_name){
    file_name_full <- sprintf('./data/processed/glfBLUP_naive_regsubset_bic/%s/%s/%s',
                              selected_env, plant_stage, file_name)

    load(file_name_full)
    return(data.frame(CV = c('CV1', 'CV2'),
                      rho = c(x$rho_cv1, x$rho_cv2)))
  })

  #wrangle accuracies to a data frame with dataset as identifier
  names(accuracies) <- str_extract(file_names, '[0-9]+')
  pa_df <- bind_rows(accuracies, .id = 'dataset')
  return(pa_df)
}, simplify = FALSE)

# provide setting names to pa list
names(naive_pa_list) <- paste(settings_naive$selected_env,
                                    settings_naive$plant_stage, sep = '_')

# by function added attributes we want to remove
naive_pa_list <- lapply(naive_pa_list, function(x){x})

# create pa df with setting columns
naive_pa_df <- bind_rows(naive_pa_list, .id = 'setting') %>%
  separate_wider_delim(setting, names = c('environment', 'plant_stage'), delim = '_')

# gBLUP ----
settings_gblup <- settings_naive %>%
  select(selected_env) %>%
  distinct()

# gblup accuracies collections
gblup_pa_list <- by(settings_gblup, seq_len(nrow(settings_gblup)), function(row){
  selected_env <-   row$selected_env

  # retrieve filenames for specific setting
  file_names <- list.files(sprintf('./data/processed/gBLUP/%s/',
                                   selected_env))

  # retrieve pa from datafile
  accuracies <- lapply(file_names, function(file_name){
    file_name_full <- sprintf('./data/processed/gBLUP/%s/%s',
                              selected_env,  file_name)

    load(file_name_full)
    return(data.frame(CV = c('univariate'),
                      rho = c(x$rho)))
  })

  #wrangle accuracies to a data frame with dataset as identifier
  names(accuracies) <- str_extract(file_names, '[0-9]+')
  pa_df <- bind_rows(accuracies, .id = 'dataset')
  return(pa_df)
}, simplify = FALSE)

# provide setting names to pa list
names(gblup_pa_list) <- paste(settings_gblup$selected_env)

# by function added attributes we want to remove
gblup_pa_list <- lapply(gblup_pa_list, function(x){x})

# create pa df with setting columns
gblup_pa_df <- bind_rows(gblup_pa_list, .id = 'environment') |>
  dplyr::select(-CV)

# timepoint factors ----
## varimax ----
varimax_pa_list <- by(settings_df, seq_len(nrow(settings_df)), function(row){
  selected_env <-   row$selected_env
  plant_stage <- row$plant_stage
  CV <- row$CV

  # retrieve filenames for specific setting
  file_names <- list.files(sprintf('./data/processed/glfBLUP_timepoints/varimax/%s/%s/%s/',
                                   selected_env, CV, plant_stage))

  # retrieve pa from datafile
  accuracies <- lapply(file_names, function(file_name){
    file_name_full <- sprintf('./data/processed/glfBLUP_timepoints/varimax/%s/%s/%s/%s',
                              selected_env, CV, plant_stage, file_name)

    load(file_name_full)
    return(data.frame(rho = x$rho))
  })

  #wrangle accuracies to a data frame with dataset as identifier
  names(accuracies) <- str_extract(file_names, '[0-9]+')
  pa_df <- bind_rows(accuracies, .id = 'dataset')
  return(pa_df)
}, simplify = FALSE)

# provide setting names to pa list
names(varimax_pa_list) <- paste(settings_df$selected_env, settings_df$CV, settings_df$plant_stage, sep = '_')

# by function added attributes we want to remove
varimax_pa_list <- lapply(varimax_pa_list, function(x){x})

# create pa df with setting columns
varimax_pa_df <- bind_rows(varimax_pa_list, .id = 'setting') %>%
  separate_wider_delim(setting, names = c('environment', 'CV', 'plant_stage'), delim = '_')

## procrustes ----
procrustes_pa_list <- by(settings_df, seq_len(nrow(settings_df)), function(row){
  selected_env <-   row$selected_env
  plant_stage <- row$plant_stage
  CV <- row$CV

  # retrieve filenames for specific setting
  file_names <- list.files(sprintf('./data/processed/glfBLUP_timepoints/procrustes/%s/%s/%s/',
                                   selected_env, CV, plant_stage))

  # retrieve pa from datafile
  accuracies <- lapply(file_names, function(file_name){
    file_name_full <- sprintf('./data/processed/glfBLUP_timepoints/procrustes/%s/%s/%s/%s',
                              selected_env, CV, plant_stage, file_name)

    load(file_name_full)
    return(data.frame(rho = x$rho))
  })

  #wrangle accuracies to a data frame with dataset as identifier
  names(accuracies) <- str_extract(file_names, '[0-9]+')
  pa_df <- bind_rows(accuracies, .id = 'dataset')
  return(pa_df)
}, simplify = FALSE)

# provide setting names to pa list
names(procrustes_pa_list) <- paste(settings_df$selected_env, settings_df$CV, settings_df$plant_stage, sep = '_')

# by function added attributes we want to remove
procrustes_pa_list <- lapply(procrustes_pa_list, function(x){x})

# create pa df with setting columns
procrustes_pa_df <- bind_rows(procrustes_pa_list, .id = 'setting') %>%
  separate_wider_delim(setting, names = c('environment', 'CV', 'plant_stage'), delim = '_')

# combined dataframe ----
# because glfblup with splines had some issues, we need to only use those datasets included
# gblup is the same across plant stages within an environment, so we give it a seperate "stage"
complete_df <- spline_pa_df %>%
  left_join(naive_pa_df, by = c('environment', 'plant_stage', 'dataset', 'CV'), suffix = c('_spline', '_naive')) %>%
  left_join(varimax_pa_df, by = c('environment', 'plant_stage', 'dataset', 'CV')) %>%
  left_join(procrustes_pa_df, by = c('environment', 'plant_stage', 'dataset', 'CV'), suffix = c('_varimax', '_procrustes')) %>%
  left_join(gblup_pa_df, by = c('environment', 'dataset')) %>%
  mutate(rho_gblup = rho, .keep = 'unused') %>%
  pivot_longer(cols = contains('rho'), names_to = 'type', values_to = 'pa') %>%
  separate_wider_delim(type, delim = '_', names = c(NA, 'type'))

## perform significance testing using Nadeau and Bengio (2003) procedure ----
# train-test folds are related over 100 subsamples, so we use Nadeau & Bengio (2003) procedure
# this is implemented in correctR::resampled_ttest
# n = number of subsamples (train-test fold splits).
# n1 = number of train genotypes 1092 - 311 = 722
# n2 = number of test genotypes 311
split_dfs <- split(filter(complete_df, type != 'varimax'), ~CV + environment + plant_stage)
model_pairs <- combn(unique(split_dfs$CV2.B2IR.full$type), m = 2)
pval_list <- lapply(split_dfs, function(x){

  wide_df <- x |>
    pivot_wider(names_from = type, values_from = pa)

  pair_df <- as.data.frame(t(model_pairs))
  pair_df$pval <- NULL

  for (i in 1:nrow(pair_df)) {
    tmp_out <- resampled_ttest(x = wide_df[[pair_df$V1[i]]],
                               y = wide_df[[pair_df$V2[i]]], n = 100, n1 = 722, n2 = 311, tailed = 'two')
    pair_df$pval[i] <- tmp_out[[2]]
  }

  pvals <- pair_df$pval
  names(pvals) <- paste(pair_df$V1, pair_df$V2, sep = '-')

  # change denominator to relevant number methods/m (here it is 4 for methods univariate, concatenated, procrustes and, spline)
  cld <- multcompView::multcompLetters(pvals, threshold = 0.05/4)$Letters
  return(cld)
})

# make df with cld letter combinations
cld_df <- bind_rows(pval_list, .id = 'environment') |>
  pivot_longer(cols = spline:gblup, names_to = 'type', values_to = 'cld') |>
  separate_wider_delim('environment', delim = '.', names = c('CV', 'environment', 'plant_stage')) |>
  group_by(CV, environment, plant_stage) |>
  # get counts per grouping of comparisons
  add_count(cld) |>
  # if it only contains one letter (a) it will sum to number of methods (4)
  filter(n != 4) |>
  # remove white space around letters and set position for plot
  mutate(cld = str_trim(cld),
         y = .9,
         .keep = 'unused')
cld_df

## create boxplot with cld all data ----
complete_signif_boxplot <- complete_df |>
  left_join(cld_df) |>
  filter(type != 'varimax') |>
  mutate(plant_stage = fct_recode(factor(plant_stage, levels = c('vegetative', 'heading', 'full')),
                                               Vegetative = 'vegetative',
                                               `Vegetative & Heading` = 'heading',
                                               All = 'full'),
                      environment = fct_recode(factor(environment, levels = c('B5IR', 'B2IR', 'DRIP')),
                                               'Optimal flat' = 'B5IR',
                                               'Moderate drought' = 'B2IR',
                                               'Severe drought' = 'DRIP'),
                      type = fct_recode(factor(type, levels = c('gblup', 'naive', 'varimax', 'procrustes', 'spline')),
                                        'Univariate gBLUP' = 'gblup',
                                        'Concatenated glfBLUP' = 'naive',
                                        'Varimax factors glfBLUP' = 'varimax',
                                        'Procrustes factors glfBLUP' = 'procrustes',
                                        'Spline characteristics glfBLUP' = 'spline')) |>
  # for plot without spline
  ggplot(aes(x = plant_stage, y = pa, fill = type)) +
  geom_boxplot() +
  geom_text(aes(y = y, label = cld, colour = type), , size = 3, position = position_dodge(width = .85), stat = 'unique') +
  facet_wrap(~CV + environment) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(fill=guide_legend(title= "Method")) +
  scale_colour_manual(values = type_colours, guide = 'none') +
  scale_fill_manual(values = type_colours) +
  labs(
    x = 'Plant stages included',
    y = 'Predictive ability')

complete_signif_boxplot

ggsave(filename = './figures/complete_boxplot_pa_signif_spline.pdf', plot = complete_signif_boxplot, width = 12, height = 10, device = 'pdf')

# boxplot with cld without splines ----
## perform significance testing using Nadeau and Bengio (2003) procedure ----
# train-test folds are related over 100 subsamples, so we use Nadeau & Bengio (2003) procedure
# this is implemented in correctR::resampled_ttest
# n = number of subsamples (train-test fold splits).
# n1 = number of train genotypes 1033 - 311 = 722
# n2 = number of test genotypes 311
nospline_split_dfs <- split(filter(complete_df, !type %in% c('spline', 'varimax')), ~CV + environment + plant_stage)
model_pairs <- combn(unique(nospline_split_dfs$CV2.B2IR.full$type), m = 2)
nospline_pval_list <- lapply(nospline_split_dfs, function(x){
  wide_df <- x |>
    pivot_wider(names_from = type, values_from = pa)

  pair_df <- as.data.frame(t(model_pairs))
  pair_df$pval <- NULL

  for (i in 1:nrow(pair_df)) {
    tmp_out <- resampled_ttest(x = wide_df[[pair_df$V1[i]]],
                               y = wide_df[[pair_df$V2[i]]], n = 100, n1 = 722, n2 = 311, tailed = 'two')
    pair_df$pval[i] <- tmp_out[[2]]
  }

  pvals <- pair_df$pval
  names(pvals) <- paste(pair_df$V1, pair_df$V2, sep = '-')

  # change denominator to relevant number methods/m (here it is 3 for methods univariate, concatenated and, procrustes)
  cld <- multcompView::multcompLetters(pvals, threshold = 0.05/3)$Letters
  return(cld)
})

nospline_cld_df <- bind_rows(nospline_pval_list, .id = 'environment') |>
  pivot_longer(cols = naive:gblup, names_to = 'type', values_to = 'cld') |>
  separate_wider_delim('environment', delim = '.', names = c('CV', 'environment', 'plant_stage')) |>
  group_by(CV, environment, plant_stage) |>
  # get counts per grouping of comparisons
  add_count(cld) |>
  # if it only contains one letter (a) it will sum to number of methods (3)
  filter(n != 3) |>
  # remove white space around letters and set position for plot
  mutate(cld = str_trim(cld),
         y = .9,
         .keep = 'unused')

nospline_cld_df

# join cld display to boxplot df ----

# create boxplot with cld ----
nospline_signif_boxplot <- complete_df |>
  filter(!type %in% c('varimax', 'spline')) |>
  left_join(nospline_cld_df) |>
  mutate(plant_stage = fct_recode(factor(plant_stage, levels = c('vegetative', 'heading', 'full')),
                                  Vegetative = 'vegetative',
                                  `Vegetative & Heading` = 'heading',
                                  All = 'full'),
         environment = fct_recode(factor(environment, levels = c('B5IR', 'B2IR', 'DRIP')),
                                  'Optimal flat' = 'B5IR',
                                  'Moderate drought' = 'B2IR',
                                  'Severe drought' = 'DRIP'),
         type = fct_recode(factor(type, levels = c('gblup', 'naive', 'varimax', 'procrustes', 'spline')),
                           'Univariate gBLUP' = 'gblup',
                           'Concatenated glfBLUP' = 'naive',
                           'Varimax factors glfBLUP' = 'varimax',
                           'Procrustes factors glfBLUP' = 'procrustes',
                           'Spline characteristics glfBLUP' = 'spline')) %>%
  # for plot without spline
  ggplot(aes(x = plant_stage, y = pa, fill = type)) +
  geom_boxplot() +
  geom_text(aes(y = y, label = cld, colour = type), , size = 3, position = position_dodge(width = .85), stat = 'unique') +
  facet_wrap(~CV + environment) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(fill=guide_legend(title= "Method")) +
  scale_colour_manual(values = type_colours, guide = 'none') +
  scale_fill_manual(values = type_colours) +
  labs(
    x = 'Plant stages included',
    y = 'Predictive ability')

nospline_signif_boxplot

ggsave(filename = './figures/complete_boxplot_pa_signif_nospline.pdf', plot = nospline_signif_boxplot, width = 12, height = 10, device = 'pdf')


# get means and CIs ----
complete_df |>
  group_by(environment, CV, plant_stage, type) |>
  mutate(mean_pa = mean(pa),
         lower_pa = mean_pa + qt(0.025, 99, lower.tail = TRUE)* sd(pa)/sqrt(100),
         upper_pa = mean_pa + qt(0.025, 99, lower.tail = FALSE) * sd(pa)/sqrt(100),
         .keep = 'unused') |>
  select(-dataset) |>
  distinct()


