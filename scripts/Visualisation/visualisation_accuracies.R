# collect accuracy scores

# libraries ----
library(tidyverse)

# settings ----
# create settings df
settings_df <- expand.grid(CV = c('CV1', 'CV2'),
                           selected_env = c('B5IR', 'B2IR', 'DRIP'),
                           plant_stage = c('full', 'heading', 'vegetative'),
                           stringsAsFactors = FALSE)

# gfBLUP with splines ----
spline_accuracy_list <- by(settings_df, seq_len(nrow(settings_df)), function(row){
  selected_env <-   row$selected_env
  CV <- row$CV
  plant_stage <- row$plant_stage

  # retrieve filenames for specific setting
  file_names <- list.files(sprintf('./data/processed/gfBLUP_spline/%s/%s/%s/',
                                   selected_env, CV, plant_stage))

  # retrieve accuracy from datafile
  accuracies <- lapply(file_names, function(file_name){
    file_name_full <- sprintf('./data/processed/gfBLUP_spline/%s/%s/%s/%s',
                              selected_env, CV, plant_stage, file_name)

    load(file_name_full)
    return(data.frame(rho = x$rho))
  })

  #wrangle accuracies to a data frame with dataset as identifier
  names(accuracies) <- str_extract(file_names, '[0-9]+')
  accuracy_df <- bind_rows(accuracies, .id = 'dataset')
  return(accuracy_df)
}, simplify = FALSE)

# provide setting names to accuracy list
names(spline_accuracy_list) <- paste(settings_df$selected_env, settings_df$CV, settings_df$plant_stage, sep = '_')

# by function added attributes we want to remove
spline_accuracy_list <- lapply(spline_accuracy_list, function(x){x})

# create accuracy df with setting columns
spline_accuracy_df <- bind_rows(spline_accuracy_list, .id = 'setting') %>%
  separate_wider_delim(setting, names = c('environment', 'CV', 'plant_stage'), delim = '_')

# get number of datasets and median value to position label in boxplot
spline_accuracy_df <- spline_accuracy_df %>%
  filter(!is.na(rho))

plot_spline_accuracy_df <- spline_accuracy_df %>%
  group_by(environment, plant_stage, CV) %>%
  mutate(median_rho = median(as.numeric(rho)),
         n_datasets = n())

## boxplot visualisation ----
spline_accuracy_boxplot <- plot_spline_accuracy_df %>%
  mutate(plant_stage = fct_recode(factor(plant_stage, levels = c('vegetative', 'heading', 'full')),
                                  Vegetative = 'vegetative',
                                  `Vegetative & Heading` = 'heading',
                                  Full = 'full'),
         environment = fct_recode(factor(environment, levels = c('B5IR', 'B2IR', 'DRIP')),
                                  'Optimal flat' = 'B5IR',
                                  'Moderate drought' = 'B2IR',
                                  'Severe drought' = 'DRIP')) %>%
  ggplot(aes(x = plant_stage, y = rho, fill = CV)) +
  geom_boxplot() +
  geom_text(aes(label = paste('n =', n_datasets), y = median_rho + 0.01),
            position = position_dodge(width = .75),
            size = 2.5) +
  facet_wrap(~environment) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(
    x = 'Plant stages included',
    y = 'Prediction accuracy')

spline_accuracy_boxplot
ggsave(filename = './figures/boxplot_accuracy.pdf', plot = spline_accuracy_boxplot, width = 5, height = 4)

# naive gfBLUP ----
settings_naive <- settings_df %>%
  select(selected_env, plant_stage) %>%
  distinct()

naive_accuracy_list <- by(settings_naive, seq_len(nrow(settings_naive)), function(row){
  selected_env <-   row$selected_env
  plant_stage <- row$plant_stage

  # retrieve filenames for specific setting
  file_names <- list.files(sprintf('./data/processed/gfBLUP_naive/%s/%s/',
                                   selected_env, plant_stage))

  # retrieve accuracy from datafile
  accuracies <- lapply(file_names, function(file_name){
    file_name_full <- sprintf('./data/processed/gfBLUP_naive/%s/%s/%s',
                              selected_env, plant_stage, file_name)

    load(file_name_full)
    return(data.frame(CV = c('CV1', 'CV2'),
                      rho = c(x$rho_cv1, x$rho_cv2)))
  })

  #wrangle accuracies to a data frame with dataset as identifier
  names(accuracies) <- str_extract(file_names, '[0-9]+')
  accuracy_df <- bind_rows(accuracies, .id = 'dataset')
  return(accuracy_df)
}, simplify = FALSE)

# provide setting names to accuracy list
names(naive_accuracy_list) <- paste(settings_naive$selected_env,
                                    settings_naive$plant_stage, sep = '_')

# by function added attributes we want to remove
naive_accuracy_list <- lapply(naive_accuracy_list, function(x){x})

# create accuracy df with setting columns
naive_accuracy_df <- bind_rows(naive_accuracy_list, .id = 'setting') %>%
  separate_wider_delim(setting, names = c('environment', 'plant_stage'), delim = '_')

## boxplot visualisation ----
naive_accuracy_boxplot <- naive_accuracy_df %>%
  mutate(plant_stage = fct_recode(factor(plant_stage, levels = c('vegetative', 'heading', 'full')),
                                  Vegetative = 'vegetative',
                                  Heading = 'heading',
                                  Full = 'full'),
         environment = fct_recode(factor(environment, levels = c('B5IR', 'B2IR', 'DRIP')),
                                  'Optimal flat' = 'B5IR',
                                  'Moderate drought' = 'B2IR',
                                  'Severe drought' = 'DRIP')) %>%
  ggplot(aes(x = plant_stage, y = rho, fill = CV)) +
  geom_boxplot() +
  # geom_text(aes(label = paste('n =', n_datasets), y = median_rho + 0.01),
  #           position = position_dodge(width = .75),
  #           size = 2.5) +
  facet_wrap(~environment) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(
    x = 'Spline length',
    y = 'Accuracy')

naive_accuracy_boxplot

# gBLUP ----
settings_gblup <- settings_naive %>%
  select(selected_env) %>%
  distinct()

# gblup accuracies collections
gblup_accuracy_list <- by(settings_gblup, seq_len(nrow(settings_gblup)), function(row){
  selected_env <-   row$selected_env

  # retrieve filenames for specific setting
  file_names <- list.files(sprintf('./data/processed/gBLUP/%s/',
                                   selected_env))

  # retrieve accuracy from datafile
  accuracies <- lapply(file_names, function(file_name){
    file_name_full <- sprintf('./data/processed/gBLUP/%s/%s',
                              selected_env,  file_name)

    load(file_name_full)
    return(data.frame(CV = c('univariate'),
                      rho = c(x$rho)))
  })

  #wrangle accuracies to a data frame with dataset as identifier
  names(accuracies) <- str_extract(file_names, '[0-9]+')
  accuracy_df <- bind_rows(accuracies, .id = 'dataset')
  return(accuracy_df)
}, simplify = FALSE)

# provide setting names to accuracy list
names(gblup_accuracy_list) <- paste(settings_gblup$selected_env)

# by function added attributes we want to remove
gblup_accuracy_list <- lapply(gblup_accuracy_list, function(x){x})

# create accuracy df with setting columns
gblup_accuracy_df <- bind_rows(gblup_accuracy_list, .id = 'environment')

## create boxplot -----
gblup_accuracy_df %>%
  mutate(environment = fct_recode(factor(environment, levels = c('B5IR', 'B2IR', 'DRIP')),
                                  'Optimal flat' = 'B5IR',
                                  'Moderate drought' = 'B2IR',
                                  'Severe drought' = 'DRIP')) %>%
  ggplot(aes(x = environment, y = rho)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(
    x = 'Spline length',
    y = 'Accuracy')

# timepoint factors ----
## varimax ----
varimax_accuracy_list <- by(settings_df, seq_len(nrow(settings_df)), function(row){
  selected_env <-   row$selected_env
  plant_stage <- row$plant_stage
  CV <- row$CV

  # retrieve filenames for specific setting
  file_names <- list.files(sprintf('./data/processed/gfBLUP_timepoints/varimax/%s/%s/%s/',
                                   selected_env, CV, plant_stage))

  # retrieve accuracy from datafile
  accuracies <- lapply(file_names, function(file_name){
    file_name_full <- sprintf('./data/processed/gfBLUP_timepoints/varimax/%s/%s/%s/%s',
                              selected_env, CV, plant_stage, file_name)

    load(file_name_full)
    return(data.frame(rho = x$rho))
  })

  #wrangle accuracies to a data frame with dataset as identifier
  names(accuracies) <- str_extract(file_names, '[0-9]+')
  accuracy_df <- bind_rows(accuracies, .id = 'dataset')
  return(accuracy_df)
}, simplify = FALSE)

# provide setting names to accuracy list
names(varimax_accuracy_list) <- paste(settings_df$selected_env, settings_df$CV, settings_df$plant_stage, sep = '_')

# by function added attributes we want to remove
varimax_accuracy_list <- lapply(varimax_accuracy_list, function(x){x})

# create accuracy df with setting columns
varimax_accuracy_df <- bind_rows(varimax_accuracy_list, .id = 'setting') %>%
  separate_wider_delim(setting, names = c('environment', 'CV', 'plant_stage'), delim = '_')

# get number of datasets and median value to position label in boxplot
varimax_accuracy_df <- varimax_accuracy_df %>%
  filter(!is.na(rho))

plot_varimax_accuracy_df <- varimax_accuracy_df %>%
  group_by(environment, plant_stage, CV) %>%
  mutate(median_rho = median(as.numeric(rho)),
         n_datasets = n())

### boxplot visualisation ----
varimax_accuracy_boxplot <- plot_varimax_accuracy_df %>%
  mutate(plant_stage = fct_recode(factor(plant_stage, levels = c('vegetative', 'heading', 'full')),
                                  Vegetative = 'vegetative',
                                  `Vegetative & Heading` = 'heading',
                                  Full = 'full'),
         environment = fct_recode(factor(environment, levels = c('B5IR', 'B2IR', 'DRIP')),
                                  'Optimal flat' = 'B5IR',
                                  'Moderate drought' = 'B2IR',
                                  'Severe drought' = 'DRIP')) %>%
  ggplot(aes(x = plant_stage, y = rho, fill = CV)) +
  geom_boxplot() +
  facet_wrap(~environment) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(
    x = 'Plant stages included',
    y = 'Prediction accuracy')

varimax_accuracy_boxplot
ggsave(filename = './figures/boxplot_varimax_accuracy.pdf', plot = spline_accuracy_boxplot, width = 5, height = 4)

## procrustes ----
procrustes_accuracy_list <- by(settings_df, seq_len(nrow(settings_df)), function(row){
  selected_env <-   row$selected_env
  plant_stage <- row$plant_stage
  CV <- row$CV

  # retrieve filenames for specific setting
  file_names <- list.files(sprintf('./data/processed/gfBLUP_timepoints/procrustes/%s/%s/%s/',
                                   selected_env, CV, plant_stage))

  # retrieve accuracy from datafile
  accuracies <- lapply(file_names, function(file_name){
    file_name_full <- sprintf('./data/processed/gfBLUP_timepoints/procrustes/%s/%s/%s/%s',
                              selected_env, CV, plant_stage, file_name)

    load(file_name_full)
    return(data.frame(rho = x$rho))
  })

  #wrangle accuracies to a data frame with dataset as identifier
  names(accuracies) <- str_extract(file_names, '[0-9]+')
  accuracy_df <- bind_rows(accuracies, .id = 'dataset')
  return(accuracy_df)
}, simplify = FALSE)

# provide setting names to accuracy list
names(procrustes_accuracy_list) <- paste(settings_df$selected_env, settings_df$CV, settings_df$plant_stage, sep = '_')

# by function added attributes we want to remove
procrustes_accuracy_list <- lapply(procrustes_accuracy_list, function(x){x})

# create accuracy df with setting columns
procrustes_accuracy_df <- bind_rows(procrustes_accuracy_list, .id = 'setting') %>%
  separate_wider_delim(setting, names = c('environment', 'CV', 'plant_stage'), delim = '_')

# get number of datasets and median value to position label in boxplot
procrustes_accuracy_df <- procrustes_accuracy_df %>%
  filter(!is.na(rho))

plot_procrustes_accuracy_df <- procrustes_accuracy_df %>%
  group_by(environment, plant_stage, CV) %>%
  mutate(median_rho = median(as.numeric(rho)),
         n_datasets = n())

### boxplot visualisation ----
procrustes_accuracy_boxplot <- plot_procrustes_accuracy_df %>%
  mutate(plant_stage = fct_recode(factor(plant_stage, levels = c('vegetative', 'heading', 'full')),
                                  Vegetative = 'vegetative',
                                  `Vegetative & Heading` = 'heading',
                                  Full = 'full'),
         environment = fct_recode(factor(environment, levels = c('B5IR', 'B2IR', 'DRIP')),
                                  'Optimal flat' = 'B5IR',
                                  'Moderate drought' = 'B2IR',
                                  'Severe drought' = 'DRIP')) %>%
  ggplot(aes(x = plant_stage, y = rho, fill = CV)) +
  geom_boxplot() +
  facet_wrap(~environment) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(
    x = 'Plant stages included',
    y = 'Prediction accuracy')

procrustes_accuracy_boxplot
ggsave(filename = './figures/boxplot_procrustes_accuracy.pdf', plot = spline_accuracy_boxplot, width = 5, height = 4)

# combined dataframe ----
# because gfblup with splines had some issues, we need to only use those datasets included
# gblup is the same across plant stages within an environment, so we give it a seperate "stage"
complete_df <- spline_accuracy_df %>%
  left_join(naive_accuracy_df, by = c('environment', 'plant_stage', 'dataset', 'CV'), suffix = c('_spline', '_naive')) %>%
  left_join(varimax_accuracy_df, by = c('environment', 'plant_stage', 'dataset', 'CV')) %>%
  left_join(procrustes_accuracy_df, by = c('environment', 'plant_stage', 'dataset', 'CV'), suffix = c('_varimax', '_procrustes')) %>%
  left_join(gblup_accuracy_df, by = c('environment', 'dataset')) %>%
  mutate(rho_gblup = rho, .keep = 'unused') %>%
  pivot_longer(cols = contains('rho'), names_to = 'type', values_to = 'accuracy') %>%
  separate_wider_delim(type, delim = '_', names = c(NA, 'type'))

head(complete_df)

# complete boxplot ----
# vector for colours ----
type_colours <- c('#E18727'= 'gblup',
                  '#BC3C29' = 'concatenated glfBLUP',
                  '#0072B5' = 'Varimax factors glfBLUP',
                  '#20854E' = 'Procrustes factors glfBLUP',
                  '#EE4C97' = 'spline characteristics glfBLUP')

type_colours_2 <- c('#E18727',
                 '#BC3C29',
                 '#0072B5',
                 '#20854E',
                 '#EE4C97')
complete_boxplot <- complete_df %>%
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
  # filter(type != 'Spline characteristics glfBLUP') %>%
  ggplot(aes(x = plant_stage, y = accuracy, fill = type)) +
  geom_boxplot() +
  facet_wrap(~CV.x + environment) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(fill=guide_legend(title= "Method")) +
  scale_fill_manual(values = type_colours_2) +
  labs(
    x = 'Plant stages included',
    y = 'Prediction accuracy')

complete_boxplot

ggsave(filename = './figures/complete_boxplot_accuracy_spline.pdf', plot = complete_boxplot, width = 10, height = 10, device = 'pdf')

