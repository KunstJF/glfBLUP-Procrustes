################################################################################
# yield box plot with and correlation between environments
################################################################################

# library ----
library(tidyverse)
library(lme4)
library(corrplot)
library(GGally)
library(patchwork)
library(ggpubr)
library(GGally)
library(tidyverse)

# data ----
yield_1415 <- readxl::read_xlsx("./data/raw/EYT_all_data_Y13_14_to_Y15_16.xlsx", sheet = 2)

# wrangle ----
dfs <- split(yield_1415, yield_1415$env)

# calculate heritabilities and BLUEs ----
yield_and_h2s <- lapply(dfs, function(df){
  df$GID <- as.factor(df$GID)
  df$trial.no <- as.factor(df$trial.no)
  df$rep <- as.factor(df$rep)
  df$subblock <- as.factor(df$subblock)

  # for some reason ASREML has a difficult time with predicting, (not enough space)
  # so we used lme4
  lmm.yield <- lmer(gy ~ GID + (1|trial.no) + (1|trial.no:rep) + (1|trial.no:rep:subblock), data = df)

  # Storing yield BLUEs for the test set:
  BLUES_yield <- lme4::fixef(lmm.yield)

  # Setting correct names:
  names(BLUES_yield) <- levels(lmm.yield@frame[["GID"]])

  # Adding intercept to all but the first genotype:
  BLUES_yield <- c(BLUES_yield[1],(BLUES_yield[2:length(BLUES_yield)] + BLUES_yield[1]))
  BLUES_df <- as.data.frame(BLUES_yield)

  return(list(BLUES_df = BLUES_df))
})

BLUE_df <- bind_cols(yield_and_h2s$B2IR$BLUES_df,
          yield_and_h2s$B5IR$BLUES_df,
          yield_and_h2s$DRIP$BLUES_df
          )

colMeans(BLUE_df)
colnames(BLUE_df) <- c('Moderate Drought', 'Optimal Flat', 'Severe Drought')
BLUE_df <- BLUE_df[, c(2, 1, 3)]

# store so we don't have to run again
write_csv(BLUE_df, './data/processed/yield_BLUE_df.csv')

BLUE_df <- read_csv('./data/processed/yield_BLUE_df.csv')
BLUE_long_df <- BLUE_df |>
  pivot_longer(cols = everything(),
               names_to = 'treatment',
               values_to = 'BLUE_gy') |>
  mutate(treatment = factor(treatment, levels = c('Optimal Flat', 'Moderate Drought', 'Severe Drought')))

# plot ----
boxplot_gy <- BLUE_long_df |>
ggplot(aes(x = treatment, y = BLUE_gy)) +
  geom_boxplot() +
  theme_classic() +
  # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(y = 'Genotypic means of grain yield (t/ha)',
       x = 'Treatment') +
  scale_x_discrete(labels = c(
    bquote('Optimal Flat'~'(h'^2~'= 0.61)'),
    bquote('Moderate Drought '~'(h'^2~'= 0.55)'),
    bquote('Severe Drought '~'(h'^2~'= 0.39)')))

boxplot_gy

ggsave('./figures/yield_blues_boxplot.pdf', plot = boxplot_gy, device = 'pdf')

# alternatively: ggpairs and patchwork to combine both ----
# okay ggpairs has limits, let's try patchwork together with ggpairs
# lower triangle plot with diag correlation texts
# boxplot with treatments
BLUE_df <- read_csv('./data/processed/yield_BLUE_df.csv')
BLUE_df <- BLUE_df[, c(1, 3, 2)]

# Function to customize single plots within ggpairs
lowerFn <- function(data, mapping, ...) {
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)

  df <- data.frame(x = x,
                   y = y)
  p <- ggscatter(df, x = 'x', y = 'y',
                 color = 'lightgrey',
                 add = 'reg.line',
                 add.params = list(color = "blue", fill = "lightgray"),
                 conf.int = FALSE)

  p + stat_cor(aes(label = after_stat(r.label)),
               cor.coef.name = expression(hat('\u03C1')))
}


scatter_tri <- ggtest <- ggpairs(BLUE_df,
                                 upper = 'blank',
                                 diag = 'blank',
                                 lower = list(continuous = lowerFn),
                                 switch = 'both',
                                 xlab = 'Genotypic means of grain yield (t/ha)') +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
scatter_tri

# some editing of the scatter plot to make it more concise
scatter_tri$yAxisLabels = scatter_tri$yAxisLabels[2:3]
scatter_tri$xAxisLabels = scatter_tri$xAxisLabels[1:2]
scatter_tri$plots = scatter_tri$plots[c(4,5, 7,8)]
scatter_tri$nrow = 2
scatter_tri$ncol = 2
scatter_tri

gy_boxplot <- BLUE_df |>
  pivot_longer(cols = everything(),
               names_to = 'treatment',
               values_to = 'BLUE_gy') |>
  mutate(treatment = factor(treatment, levels = c('Optimal Flat', 'Moderate Drought', 'Severe Drought'))) |>
  ggplot(aes(x = treatment, y = BLUE_gy)) +
  geom_boxplot() +
  theme_classic() +
  # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(y = 'Genotypic means of grain yield (t/ha)',
       x = 'Treatment')

gy_boxplot

cairo_pdf('./figures/yield_blue_pairsplot.pdf', width = 15, height = 10)
gy_boxplot + wrap_elements(ggmatrix_gtable(scatter_tri)) + plot_annotation(tag_levels = c('a', 'b'))
dev.off()


