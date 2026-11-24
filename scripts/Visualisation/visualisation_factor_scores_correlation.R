################################################################################
# visualisation of factor score correlation across environments
################################################################################

# libraries ----
library(tidyverse)
library(corrplot)
library(GGally)
library(ggpubr)

# data ----
load('./data/processed/factors/B5IR/timepoint_factors/procrustes/procrustes_1.RData')
HS_B5IR <- procrustes_factor_results$CV2

load('./data/processed/factors/B2IR/timepoint_factors/procrustes/procrustes_1.RData')
HS_B2IR <- procrustes_factor_results$CV2

load('./data/processed/factors/DRIP/timepoint_factors/procrustes/procrustes_1.RData')
HS_DRIP <- procrustes_factor_results$CV2
rm(procrustes_factor_results)

# wrangle data ----
HS <- list(`Optimal Flat` = lapply(HS_B5IR, function(x){
  x$factor_scores %>%
    select(G, F1, F2)}),
     `Moderate Drought` = lapply(HS_B2IR, function(x){
       x$factor_scores %>%
         select(G, F1, F2)}),
     `Severe Drought` = lapply(HS_DRIP, function(x){
       x$factor_scores %>%
         select(G, F1, F2)}))

HS_dfs <- lapply(HS, bind_rows, .id = 'timepoint')


HS_df <- bind_rows(HS_dfs, .id = 'environment') %>%
  pivot_longer(cols = F1:F2, names_to = 'factor', values_to = 'factor_score') %>%
  mutate(factor = gsub('F', 'Factor ', factor)) %>%
  group_by(environment, G, factor) %>%
  mutate(mean_factor_score = mean(factor_score)) %>%
  ungroup()

mean_df <- HS_df %>%
  select(G, environment, factor, mean_factor_score) %>%
  distinct() %>%
  unite('fac_env', environment, factor, sep = ' ') %>%
  pivot_wider(id_cols = G, names_from = fac_env, values_from = mean_factor_score)

# make similar yield plot ----
## adjust df so more readable ----
# remove genotype column and reorder for better clusters in plot
plot_df <- mean_df[, -1]
plot_df <- plot_df[, c(1, 3, 5, 2, 4, 6)]

# rename columns so they are readable in y labels
colnames(plot_df) <- str_replace(colnames(plot_df), 'Factor ', 'F')


# Function to customize single plots within ggpairs
lowerFn <- function(data, mapping, ...) {
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)

  df <- data.frame(x = x,
                   y = y)
  p <- ggpubr::ggscatter(df, x = 'x', y = 'y',
                 color = 'lightgrey',
                 add = 'reg.line',
                 add.params = list(color = "blue", fill = "lightgray"),
                 conf.int = FALSE)

  p + stat_cor(aes(label = after_stat(r.label)),
               cor.coef.name = expression(hat('\u03C1')))
}

scatter_tri <- ggpairs(plot_df,
                       upper = 'blank',
                       diag = 'blank',
                       lower = list(continuous = lowerFn),
                       switch = 'both',
                       xlab ='Genotypic mean of factor scores across timepoints in treatment' ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

scatter_tri

scatter_tri$yAxisLabels = scatter_tri$yAxisLabels[2:6]
scatter_tri$xAxisLabels = scatter_tri$xAxisLabels[1:5]
scatter_tri$plots = scatter_tri$plots[c(7:11, 13:17, 19:23, 25:29, 31:35)]
scatter_tri$nrow = 5
scatter_tri$ncol = 5
scatter_tri

cairo_pdf('./figures/factors_pairs_plot.pdf', width = 15, height = 10)
scatter_tri
dev.off()
