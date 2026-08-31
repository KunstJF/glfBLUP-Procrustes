################################################################################
# Visualisation of subset selection of factors across folds per environment
################################################################################

# libaries ----
library(tidyverse)

# setting grid ----
settings_df <- expand.grid(CV = c('CV1', 'CV2'),
                           selected_env = c('B5IR', 'B2IR', 'DRIP'),
                           plant_stage = c('full', 'heading', 'vegetative'),
                           stringsAsFactors = FALSE)
# result list
result_list  <- vector(mode = 'list', length = nrow(settings_df))
names(result_list) <- unite(settings_df, col = setting, sep = '_')$setting

# collect selected factor per setting ----
for(row_i in 1:nrow(settings_df)){
  env <- settings_df[row_i,]$selected_env
  cv <- settings_df[row_i,]$CV
  plant_stage <- settings_df[row_i,]$plant_stage

  files <- list.files(sprintf('./data/processed/gfBLUP_timepoints/procrustes/%s/%s/%s/', env, cv, plant_stage), full.names = TRUE)

  selected_facs <- lapply(files, function(f){
    load(f)
    fac_df <- as.data.frame(str_split(x$factorSelection, '_', simplify = TRUE))
    colnames(fac_df) <- c('factor', 'date')
    fac_df
  })
  result_list[[row_i]] <- bind_rows(selected_facs, .id = 'dataset')
}

# wrangle result df -----
result_df <- bind_rows(result_list, .id = 'setting') %>%
  separate_wider_delim(setting, delim = '_', names = c('CV', 'environment', 'plant_stage')) %>%
  mutate(plant_stage = fct_recode(factor(plant_stage, levels = c('vegetative', 'heading', 'full')),
                             Vegetative = 'vegetative',
                             `Vegetative & Heading` = 'heading',
                             All = 'full'),
    environment = fct_recode(factor(environment, levels = c('B5IR', 'B2IR', 'DRIP')),
                             'Optimal flat' = 'B5IR',
                             'Moderate drought' = 'B2IR',
                             'Severe drought' = 'DRIP'),
  )

result_df %>%
  filter(CV == 'CV1') %>%
  ggplot(aes(x = date, fill = factor)) +
  geom_bar(position = 'dodge') +
  facet_wrap(~environment + plant_stage) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

result_df %>%
  filter(CV == 'CV1') %>%
  ggplot(aes(x = factor, fill = date)) +
  geom_bar(position = 'dodge') +
  facet_wrap(~plant_stage + environment) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


result_df %>%
  filter(CV == 'CV1') %>%
  ggplot(aes(x = date, fill = factor)) +
  geom_bar(position = 'dodge') +
  facet_wrap(~plant_stage + environment) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# make it into a heatmap of sorts ----
library(magrittr)
library(RColorBrewer)
get_plotdf <- function(occurence_df){
  selected_df <- occurence_df %>%
    pivot_longer(cols = -1, names_to = 'factor', values_to = 'presence') %>%
    group_by(factor) %>%
    mutate(occurence = sum(presence)) %>%
    ungroup()

  # fac_ranks <- selected_df %>%
  #   select(factor, occurence) %>%
  #   distinct() %>%
  #   arrange(desc(occurence)) %>%
  #   mutate(factor_rank = 1:nrow(.),
  #          factor = fct_reorder(.f = as.factor(factor),
  #                               .x = factor_rank))
  # USE levels = levels(fac_ranks$factor) if needed

  plot_df <- selected_df %>%
    mutate(`Timepoint Factor` = factor(factor),
           Dataset = as.numeric(dataset),
           Selected = as.logical(presence), .keep = 'unused')

  return(plot_df)
}

severe_df<- result_df %>%
  filter(plant_stage == 'All' & CV == 'CV1' & environment == 'Severe drought') %>%
  unite(factor, date, col = 'date_fac') %>%
  mutate(n = 1) %>%
  select(dataset, date_fac, n) %>%
  pivot_wider(names_from = date_fac, values_from = n, values_fill = 0, names_sort = TRUE)

optimal_df <- result_df %>%
  filter(plant_stage == 'All' & CV == 'CV1' & environment == 'Optimal flat') %>%
  unite(factor, date, col = 'date_fac') %>%
  mutate(n = 1) %>%
  select(dataset, date_fac, n) %>%
  pivot_wider(names_from = date_fac, values_from = n, values_fill = 0, names_sort = TRUE)

moderate_df <- result_df %>%
  filter(plant_stage == 'All' & CV == 'CV1' & environment == 'Moderate drought') %>%
  unite(factor, date, col = 'date_fac') %>%
  mutate(n = 1) %>%
  select(dataset, date_fac, n) %>%
  pivot_wider(names_from = date_fac, values_from = n, values_fill = 0, names_sort = TRUE)

# generate occurence dfs
res_list <- list(`Severe drought` = get_plotdf(severe_df),
     `Moderate drought` = get_plotdf(moderate_df),
     `Optimal flat` = get_plotdf(optimal_df))

## create large dataframe for facet_wrap ggplot ----
### colors ----
colors <- brewer.pal(3, "Dark2")[c(1,2, 6)]

### visualise ----
plot_df <- bind_rows(res_list, .id = 'Environment') %>%
  separate_wider_delim(`Timepoint Factor`, names = c('Factor', 'Timepoint'), delim = '_') %>%
  # adds rows for timepoints of which no flights where recorded
  complete(Environment, Factor, Timepoint, Dataset) %>%
  # create nice plant stage variable for colouring
  mutate(Environment = factor(Environment, levels = c('Optimal flat', 'Moderate drought', 'Severe drought')),
         `Plant stage` = case_when(
           Environment != 'Optimal flat' & as.numeric(Timepoint) < 150209 & Selected ~ 'Vegetative',
           Environment != 'Optimal flat' & as.numeric(Timepoint) == 150209 & Selected ~ 'Heading',
           Environment != 'Optimal flat' & as.numeric(Timepoint) > 150209 & Selected ~ 'Grain filling',
           Environment == 'Optimal flat' & as.numeric(Timepoint) < 150225 ~ 'Vegetative',
           Environment == 'Optimal flat' & as.numeric(Timepoint) == 150225 ~ 'Heading',
           Environment == 'Optimal flat' & as.numeric(Timepoint) > 150225 ~ 'Grain filling',
           # complete provides NA for thos timepoints without flights in Selected
           is.na(Selected) ~ 'Flight unavailable'
         ),
         Timepoint = gsub("(\\d{2})(\\d{2})(\\d{2})$","20\\1-\\2-\\3", Timepoint),
         Factor = gsub(pattern = 'F', replacement = 'Factor ', x = Factor),
         # change Selected so that they are show in the final plot
         Selected = case_when(is.na(Selected) ~ TRUE,
                              .default = Selected),
         .keep = 'unused') %>%
  filter(Selected)

(plot_timepoint_y <- plot_df %>%
  ggplot(aes(x = Dataset, y = Timepoint)) +
  geom_raster(aes(fill = `Plant stage`)) +
  facet_wrap(~ Factor + Environment) +
  theme_classic() +
  scale_fill_manual(
    name = 'Plant stage',
    labels = c('Flight unavailable', 'Grain filling', 'Heading', 'Vegetative'),
    values = brewer.pal(8, "Dark2")[c(8, 2, 6, 1)]))

ggsave('./figures/subset_selection_y.pdf', plot = plot_timepoint_y, device = 'pdf')

(plot_timepoint_x <- plot_df %>%
  ggplot(aes(x = Timepoint, y = Dataset)) +
  geom_raster(aes(fill = `Plant stage`)) +
  facet_wrap(~Environment + Factor, nrow = 3, ncol = 2) +
  theme_classic() +
  scale_fill_manual(
    name = 'Plant stage',
    labels = c('Flight unavailable', 'Grain filling', 'Heading', 'Vegetative'),
    values = brewer.pal(8, "Dark2")[c(8, 2, 6, 1)]) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8)))

ggsave('./figures/subset_selection_x.pdf', plot = plot_timepoint_x, device = 'pdf')

# heatmap style plot ----
# create adjacency matrix
# adj_mat <- occurence_df %>%
#   select(-dataset) %$%
#   crossprod(as.matrix(.))
#
# # remove diagonal
# adj_mat[lower.tri(adj_mat)] <- NA
# ordering needs to happen on original adjacency matrix
# as.data.frame(adj_mat) %>%
#   rownames_to_column(var = 'factor_x') %>%
#   pivot_longer(cols = -1, names_to = 'factor_y', values_to = 'adjacency') %>%
#   mutate(factor_y = fct_reorder(.f = as.factor(factor_y),
#                                 .x = adjacency),
#          factor_x = fct_reorder(.f = as.factor(factor_x),
#                                 .x = adjacency)) %>%
#   ggplot(aes(x = factor_x, y = factor_y)) +
#   geom_raster(aes(fill = adjacency)) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))




