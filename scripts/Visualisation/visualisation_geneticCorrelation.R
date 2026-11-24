# libraries ----
library(emmeans)
library(multcomp)
library(tidyverse)
library(magrittr)
library(correctR)

# settings ----
# create settings df
settings_df <- expand.grid(CV = c('CV1','CV2'),
                           selected_env = c('B5IR', 'B2IR', 'DRIP'),
                           plant_stage = c('full', 'heading', 'vegetative'),
                           stringsAsFactors = FALSE)

# procrustes ----
procrustes_sg_list <- by(settings_df, seq_len(nrow(settings_df)), function(row){
  selected_env <-   row$selected_env
  plant_stage <- row$plant_stage
  CV <- row$CV

  # retrieve filenames for specific setting
  file_names <- list.files(sprintf('./data/processed/glfBLUP_timepoints/procrustes/%s/%s/%s/',
                                   selected_env, CV, plant_stage), full.names = FALSE)

  # retrieve genetic correlations from datafile
  Sg <- lapply(file_names, function(file_name){
    file_name_full <- sprintf('./data/processed/glfBLUP_timepoints/procrustes/%s/%s/%s/%s',
                              selected_env, CV, plant_stage, file_name)

    load(file_name_full)
    return((as.data.frame(t(x$Sg[,'Y']))))
  })

  # wrangle correlations to a data frame with dataset as identifier
  names(Sg) <- str_extract(file_names, '[0-9]+')
  sg_df <- bind_rows(Sg, .id = 'dataset')
  return(sg_df)
}, simplify = FALSE)


# provide setting names to genetic correlations list
names(procrustes_sg_list) <- paste(settings_df$selected_env, settings_df$CV, settings_df$plant_stage, sep = '_')

# by function added attributes we want to remove
procrustes_sg_list <- lapply(procrustes_sg_list, function(x){x})

# create accuracy df with setting columns
procrustes_sg_df <- bind_rows(procrustes_sg_list, .id = 'setting') %>%
  separate_wider_delim(setting, names = c('environment', 'CV', 'plant_stage'), delim = '_')

# plotting ----
plot_df <- procrustes_sg_df |>
  pivot_longer(F1_150110:F2_150315, names_to = 'variable', values_to = 'correlation') |>
  filter(variable != 'Y') |>
  separate_wider_delim(variable, names = c('Factor', 'timepoint'), delim = '_') |>
  mutate(timepoint = paste0('20',str_replace_all(timepoint, "(\\d{2})(?=\\d)", "\\1-")),
         environment = fct_recode(factor(environment, levels = c('B5IR', 'B2IR', 'DRIP')),
                                  'Optimal flat' = 'B5IR',
                                  'Moderate drought' = 'B2IR',
                                  'Severe drought' = 'DRIP'),
         Factor = str_remove(Factor, 'F'))

## boxplot ----
boxplot_correlation <- plot_df |>
  ggplot(aes(x = timepoint, y = correlation)) +
  geom_boxplot(aes(fill = Factor)) +
  facet_wrap(~CV + environment, scales = 'free_x') +
  xlab('Date') +
  ylab('Genetic correlation with grain yield') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


boxplot_correlation

## save plot ----
ggsave(filename = './figures/boxplot_geneticCorrelation.pdf',
       plot = boxplot_correlation, width = 12, height = 10, device = 'pdf')

## summary stats ----
plot_df |>
  group_by(environment, CV, plant_stage, Factor, timepoint) |>
  mutate(mean_cor = mean(correlation)) |>
  select(-dataset, -correlation) |>
  distinct() |>
  filter(CV == 'CV2' & environment == 'Severe drought')




