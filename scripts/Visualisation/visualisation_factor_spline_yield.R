################################################################################
# Visualisation of factor score across time using splines
# with addition of ranking top and bottom yielding genotypes
################################################################################


# libraries ----
library(tidyverse)
library(statgenHTP)
library(wesanderson)
library(corrplot)

# localisation for time plotting ----
Sys.setlocale("LC_ALL", "English")

# data ----
## spline data ----
selected_factor <- 1
mid_timepoint <- 5
selected_env <- 'Bed5IR'
selected_year <- 2014

fname <- paste0('./data/processed/splines/all_data/predsplinedatall_factor',
                selected_factor,
                '_mid_timepoint_',
                mid_timepoint,
                '_',
                selected_year,
                '_',
                selected_env,
                '.RDS')
spline_factor <- read_rds(fname)

## yield data ----
yield_pseudo <- read_csv('./data/processed/2014_yield_pseudo_blue.csv')
yield_pseudo <- yield_pseudo %>%
  mutate(GID = as.character(GID)) %>%
  mutate(env = case_when(env == 'B5IR' ~ 'Bed5IR',
                         env == 'B2IR' ~ 'Bed2IR',
                         env == 'DRIP' ~ 'Drip'))

# yield wrangling ----
check_genos <- table(yield_pseudo$GID)[1:2]
unique(yield_pseudo$env)
# filter yield data for correct year and environment
yield_pseudo_select <- yield_pseudo %>%
  filter(cycle == selected_year & env == selected_env & !GID %in% names(check_genos)) %>%
  group_by(GID) %>%
  summarise(mean_gy = mean(gy)) %>%
  filter(mean_gy <= quantile(probs = 0.1, mean_gy) | mean_gy >= quantile(probs = 0.9, mean_gy))

best_yield <- yield_pseudo %>%
  filter(cycle == selected_year & env == selected_env & !GID %in% names(check_genos)) %>%
  group_by(GID) %>%
  summarise(mean_gy = mean(gy)) %>%
  slice_max(mean_gy, n = 10) |>
  arrange(mean_gy)

least_yield <- yield_pseudo %>%
  filter(cycle == selected_year & env == selected_env & !GID %in% names(check_genos)) %>%
  group_by(GID) %>%
  summarise(mean_gy = mean(gy)) %>%
  slice_min(mean_gy, n = 10) |>
  arrange(desc(mean_gy))

# colouring yield in splines ----
## creating dataframe with yields included ----
geno_spline_df <- spline_factor$genoLevel

geno_yields <- yield_pseudo %>%
  filter(GID %in% geno_spline_df$genotype) %>%
  select(GID, gy) %>%
  group_by(GID) %>%
  summarise(mean_gy = mean(gy)) %>%
  mutate(genotype = as.character(GID))

test_spline_df <- geno_spline_df %>%
  left_join(geno_yields, by = 'genotype')

## plotting ----
### top/bottom 10 plots ----
# factor 1
# check correlations with yield -----
# from df
yield_2014br5ir <- yield_pseudo %>%
  filter(env == 'B5IR' & cycle == 2014) %>%
  select(GID, gy) %>%
  group_by(GID) %>%
  summarise(mean_gy = mean(gy))

yield_cortest <- spline_factor$genoLevel %>%
  mutate(GID = as.character(genotype)) %>%
  left_join(yield_2014br5ir, by = 'GID') %>%
  select(GID, fGeno:fGenoDevDeriv2, mean_gy)

# get spline param estimates from estimateSplineParameters function
spline_params_factor_full <- estimateSplineParameters(spline_factor,
                         estimate = 'predictions',
                         what = c('min', 'max', 'AUC','p75', 'p50', 'p25'),
                         fitLevel = 'geno'
                         )

# Largest variance detection idea ----
test <- test_spline_df %>%
  group_by(timeNumber) %>%
  mutate(splinevariance_t = var(fGeno))

max_pos_x <- test %>%
  select(timeNumber, splinevariance_t) %>%
  distinct() %>%
  filter(splinevariance_t == max(.$splinevariance_t))

test %>%
  ggplot(aes(x = timeNumber, y = splinevariance_t)) +
  geom_point()

as.Date.POSIXct(test_spline_df$timePoint)
as.Date(test_spline_df$timePoint)
str(test_spline_df)

# dates of measurement to indicate in plot ----
measure_dates_df <- tibble(
  date = as_datetime(c('2015-01-10 00:00:00 UTC', '2015-01-19 00:00:00 UTC', '2015-02-04 00:00:00 UTC',
                '2015-02-09 00:00:00 UTC', '2015-02-25 00:00:00 UTC', '2015-03-10 00:00:00 UTC',
                '2015-03-15 00:00:00 UTC', '2015-03-23 00:00:00 UTC', '2015-03-26 00:00:00 UTC',
                '2015-04-07 00:00:00 UTC')),
  physiological_stage = c(rep('Vegetative', 4), 'Heading', rep('Grain filling', 5))
)



# final plot ----
final_plot <- test_spline_df %>%
  mutate(top_bottom_10 = case_when(genotype %in% least_yield$GID ~ 'Bottom 10 genotypes',
                                   genotype %in% best_yield$GID ~ 'Top 10 genotypes',
                                   .default = 'Mid')) |>
  ggplot(aes(x = .data[["timePoint"]])) +
  geom_line(aes(y = .data[["fGeno"]],
                group = .data[["genotype"]],
                colour = .data[['top_bottom_10']],
                alpha = .data[['top_bottom_10']],
                linewidth = .data[['top_bottom_10']]), na.rm = TRUE) +
  geom_hline(yintercept = 0) +
  geom_vline(data = measure_dates_df, aes(xintercept = date), linetype = 2) +
  geom_label(data = measure_dates_df, aes(x = date, y = 6, label = physiological_stage, angle = -90), linewidth = 0) +
  # if you want to visualise the point with most variance
  # geom_vline(xintercept = window_range['mid']) +
  # # geom_rect(aes(xmin = window_range['min'],
  # #           xmax = window_range['max'],
  # #           ymin = -Inf,
  # #           ymax = Inf), alpha = 0.005) +
  scale_color_manual(values = wes_palette('Cavalcanti1')[c(5,1,2)]) +
  scale_alpha_manual(values = c(top = 1,  bottom = 1, Mid = .1), guide = 'none') +
  scale_discrete_manual("linewidth", values =c(1.5, 0.75, 1.5), guide = 'none') +
  scale_x_datetime(date_labels = '%b') +
  labs(y = paste('Genotype specific spline of factor', selected_factor),
       x = 'Time',
       title = NULL,
       color = 'Yield rank') +
  theme_classic()

# store plot ----
ggsave(filename = './figures/factor_spline.jpeg', plot = final_plot, width = 12, height = 10, device = 'jpeg')

