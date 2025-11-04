# CIMMYT proc factor splines example with statgenHTP

# libraries
library(tidyverse)
library(statgenHTP)
library(statgenGWAS)
library(FMradio)
library(gfBLUP)
library(EFA.dimensions)
library(DandEFA)
library(viridis)

# custom functions ----
# function to generate factors for heading timepoint
generate_factor <- function(timepoint, df, m = NULL){
  df <- select(df, G, contains(timepoint))
  sumsq <- covSS(df)

  data.means <- stats::aggregate(df[, 2:ncol(df)], list(df[, 1]), mean)
  names(data.means)[1] <- names(df)[1]
  data.means$G <- as.character(data.means$G)

  # create factors
  cor_g <- cov2cor(sumsq$Sg)
  cor_e <- cov2cor(sumsq$Se)

  # cor_g <- RF(cor_g)
  HS_factor <- factorModel(df, cormat = cor_g, what = 'genetic', m = m)
  return(HS_factor)
}

# function to perform procrustes rotation
procustes <- function (loadings, target, type = "orthogonal", verbose = TRUE){
  if (type == "orthogonal") {
    s = t(loadings) %*% target
    w1 = s %*% t(s)
    v1 = t(s) %*% s
    w <- eigen(w1)$vectors
    ew <- diag(eigen(w1)$values)
    v <- eigen(v1)$vectors
    ev <- diag(eigen(v1)$values)
    o = t(w) %*% s %*% v
    k = diag(((diag(o))/abs(diag(o))), nrow = nrow(o), ncol = nrow(o))
    ww = w %*% k
    t1 = ww %*% t(v)
    procrust = loadings %*% t1
  }
  if (type == "oblique") {
    cat1 = (solve(t(loadings) %*% loadings)) %*% t(loadings) %*%
      target
    cat2 = sqrt(colSums(cat1 * cat1))
    numrows = nrow(cat1) - 1
    cat3 = cat2
    for (t in 1:numrows) {
      cat3 = rbind(cat3, cat2)
    }
    procrust = loadings %*% (cat1/cat3)
  }
  colnames(procrust) <- colnames(loadings)
  rtproc = sum(target * procrust)/sqrt(sum(target * target) *
                                         sum(procrust * procrust))
  resid = target - procrust
  nvars = nrow(target)
  pertarget = sum(target * target)/nvars
  perresid = sum((target - procrust) * (target - procrust))/nvars
  rmsr = sqrt(sum(resid * resid)/(nvars * ncol(target)))
  if (verbose == TRUE) {
    if (type == "orthogonal") {
      message("\n\nOrthogonal Procrustes Rotation")
    }
    if (type == "oblique") {
      message("\n\nOblique Procrustes Rotation:")
    }
    message("\nOriginal Loadings:\n")
    print(round(loadings, 3), print.gap = 3)
    message("\nTarget:\n")
    print(round(target, 3), print.gap = 3)
    message("\nProcrustes-Rotated Loadings:\n")
    print(round(procrust, 3), print.gap = 3)
    message("\nCongruence Between the Target & Procrustes-Rotated Loadings = ",
            round(rtproc, 3))
    message("\nProportion of Variance in Target = ", round(pertarget,
                                                           3))
    message("\nProportion of Variance in Residual = ", round(perresid,
                                                             3))
    message("\nRoot Mean Square Residual = ", round(rmsr,
                                                    3), "\n")
  }
  procrustesOutput <- list(loadingsPROC = procrust, congruence = rtproc,
                           rmsr = rmsr, residmat = resid, rotmatrix = t1)
  return(invisible(procrustesOutput))
}

# rotation matrix smoother (from Carel)
procustes_smoother <- function(rotation_matrix){
  m <- dim(rotation_matrix)[1]
  perm_mat = matrix(0, nrow = m, ncol = m)
  for (k in 1:m){
    k_plus <- which.max(rotation_matrix[, k])
    k_minus <- which.min(rotation_matrix[, k])
    r_plus <- max(rotation_matrix[k, ])
    r_minus <- min(rotation_matrix[k, ])
    if (r_minus < r_plus & r_minus > 0 & r_plus > 0){
      perm_mat[k_plus, k] = 1
    } else if (abs(r_minus) < r_plus & r_minus < 0 & r_plus > 0){
      perm_mat[k_plus, k]  = 1
    } else if (abs(r_minus) > r_plus & r_minus < 0 & r_plus > 0){
      perm_mat[k_minus, k] = -1
    } else{
      perm_mat[k_minus, k] = -1
    }
  }
  return(perm_mat)
}

# data ----
HS_files <- str_subset(list.files('./data/processed/', full.names = TRUE), 'Y(1415|1516)')
HS_files

# MANUAL AND DESCRIPTION ----
# description of environment identifiers
# Bed2IR: Moderate Drought, two furrow irrigations
# Melga: Optimal Bed, five furrow irrigations
# Bed5IR: Optimal Flat, Five Flood irrigations
# Drip: Severe Drought, Three minimal irrigations through drip
# Heat: Heat, late february planting five furrow irrigations

# Optimal Flat: 25th of February is heading date (5th timepoint)
# Severe/Moderate Drought: 9th of February is heading date (4th timepoint)

# STILL MANUAL: check environment and year
selected_year <- 2014
selected_env <- 'Drip'
selected_factor <- 'factor1' # which factor to fit a spline through and save
mid_time <- 4 # which timepoint to use



## HS wrangling ----

fname <- str_detect(HS_files, paste0('./data/processed/Y', selected_year-2000, '15_', selected_env, '_pseudo_blue.csv'))

temp_df <- as.data.frame(read_csv(HS_files[fname]))
HS_data_cols <- str_detect(colnames(temp_df), 'nm')
temp_df[, HS_data_cols] <- scale(temp_df[, HS_data_cols])
colnames(temp_df)[str_detect(colnames(temp_df), 'gid')] <- 'G'

sec_df <- temp_df

# factor creation only requires secondary traits to only contain g + traits, no other columns
HS_df <- sec_df %>%
  select('G', contains('nm'))

## wrangling for different timepoint factors ----
# get dates
timepoints <- unique(str_split(colnames(HS_df)[-1], pattern = '_', simplify = TRUE)[,2])
timepoints

# first fit factor for best timepoint (here heading date 25th of february)
factor_target <- generate_factor(timepoint = timepoints[mid_time], df = HS_df, m = 2)

# target matrix has m factors, so we generate m factors for the other timepoints as well
factors <- lapply(timepoints, HS_df = HS_df, FUN = function(x, HS_df, m = factor_target$m){
  df <- select(HS_df, G, contains(x))
  sumsq <- covSS(df)

  data.means <- stats::aggregate(df[, 2:ncol(df)], list(df[, 1]), mean)
  names(data.means)[1] <- names(df)[1]
  data.means$G <- as.character(data.means$G)

  # create factors ----
  cor_g <- cov2cor(sumsq$Sg)
  cor_e <- cov2cor(sumsq$Se)

  # cor_g <- RF(cor_g)
  HS_factor <- factorModel(df, cormat = cor_g, what = 'genetic', m = m)
  return(HS_factor)
})

# create procustes rotated factors for each timepoint
proc_factors <- lapply(factors[-mid_time], function(x){
  procustes(x$loadings, factor_target$loadings, type = 'orthogonal', verbose = TRUE)
})

# Procrustean smoothing ----
rotmat_origin <- proc_factors[[mid_time]]$rotmatrix

smoothed <- lapply(proc_factors, function(x){
  procustes_smoother(x$rotmatrix)
})

# last timepoint got columns switched
smoothed_factors <- map2(factors[-mid_time], smoothed, function(factor_object, rotation_matrix){
  factor_object$loadings %*% rotation_matrix
})

origin_mat <- as.matrix(factors[[mid_time]]$loadings[,])
colnames(origin_mat) <- gsub(x = colnames(origin_mat), pattern = 'F', replacement =  'V')

smoothed_factors <- append(smoothed_factors, list(origin_mat), after = mid_time-1)
names(smoothed_factors) <- timepoints

# visualisation of factors ----
# create overal df ----
smoothed_factor_df <- lapply(smoothed_factors, as.data.frame) %>%
  bind_rows(.id = 'timePoint') %>%
  rownames_to_column('trait') %>%
  separate_wider_regex(trait, patterns = c('nm', wavelength = '[0-9]+', '_', '[0-9]+')) %>%
  rename_with(~ gsub('V', 'Factor_', .x)) %>%
  pivot_longer(contains('Factor'), names_to = 'factor', values_to = 'loading') %>%
  separate(factor, into = c(NA, 'factor')) %>%
  mutate(wavelength = as.numeric(wavelength))

## add congruence ----
names(proc_factors) <- timepoints[-mid_time]

congruences <- lapply(proc_factors, function(x){
  x$congruence
}) %>%
  bind_rows() %>%
  pivot_longer(everything(), names_to = 'timePoint', values_to = 'congruence')

smoothed_factor_df <- left_join(smoothed_factor_df, congruences, by = 'timePoint') %>%
  mutate(timePoint_label = case_when(is.na(congruence) ~ paste0(timePoint, ' (', 1, ')'),
                                     .default = paste0(timePoint, ' (', round(congruence, 2), ')')))


# visualise ----
smoothed_factor_df %>%
  ggplot(aes(x = wavelength, y = loading, colour = factor)) +
  geom_point() +
  geom_line() +
  facet_wrap(~timePoint_label) +
  theme_classic()

# create factor score ----
proc_scores <- vector('list', length = length(timepoints))

for(timepoint_i in 1:length(timepoints)){
  temp_matrix <- select(HS_df, contains(timepoints[timepoint_i])) %>%
    as.matrix(.) %>%
    scale(.)

  rownames(temp_matrix) <- HS_df$G

  proc_scores[[timepoint_i]] <- facScore(temp_matrix, smoothed_factors[[timepoint_i]], factors[[timepoint_i]]$uniquenesses)
}

names(proc_scores) <- timepoints

# something is off with the number of genotypes generated
proc_score_df <- bind_cols(proc_scores) %>%
  rename_with( ~paste('factor', 1:factor_target$m, rep(timepoints, each = factor_target$m), sep = '_')) %>%
  bind_cols(sec_df[, 1:8]) %>%
  pivot_longer(cols = contains('factor'), values_to = 'factor_score', names_to = 'factor') %>%
  separate_wider_delim(factor, delim = '_', names = c(NA, 'factor', 'timePoint')) %>%
  mutate(timePoint = as.numeric(timePoint))

# store procrustes factor scores ----
selected_year
selected_env
mid_time

write_csv(proc_score_df, paste('./data/processed/proc_factor_score',
                               'mid_time',
                               mid_time,
                               selected_year,
                               selected_env,
                               '.csv', sep = '_'))
# statgenHTP ----
# plot for the 10% worst and best yielding genotypes
yield_pseudo <- read_csv('./data/processed/2014_yield_pseudo_blue.csv')

# create statgenhtp compatible df
statgen_df <- proc_score_df %>%
  pivot_wider(names_from = factor, values_from = factor_score, names_prefix = 'factor') %>%
  mutate(G = as.factor(G),
         timePoint = as.POSIXct(ymd(timePoint)),
         pop = 'population1',
         colId = col,
         rowId = row,
         plotId = plot)

# for only data up to and including flowering
unique(proc_score_df$timePoint)[5]
statgen_flower_df <- proc_score_df %>%
  filter(timePoint <= unique(proc_score_df$timePoint)[5]) %>%
  pivot_wider(names_from = factor, values_from = factor_score, names_prefix = 'factor') %>%
  mutate(G = as.factor(G),
         timePoint = as.POSIXct(ymd(timePoint)),
         pop = 'population1',
         colId = col,
         rowId = row,
         plotId = plot)

# for only data up to veg
unique(proc_score_df$timePoint)[4]
statgen_veg_df <- proc_score_df %>%
  filter(timePoint <= unique(proc_score_df$timePoint)[4]) %>%
  pivot_wider(names_from = factor, values_from = factor_score, names_prefix = 'factor') %>%
  mutate(G = as.factor(G),
         timePoint = as.POSIXct(ymd(timePoint)),
         pop = 'population1',
         colId = col,
         rowId = row,
         plotId = plot)


check_genos <- sort(table(statgen_df$G), decreasing = TRUE)[1:2]
write_csv(statgen_df, './data/processed/CIMMYT_statgenhtp_example_df.csv')
# filter yield data for correct year and environment
## spline fitting ----
# fit spline for factors ------
# genotypes = c(levels(statgen_df$G)[1:n_geno], as.character(unique(yield_pseudo_select$GID))),
# n_geno = 750
spline_all_fac <- fitSplineHDM(inDat = statgen_df,
                           trait = selected_factor,
                           genotype = 'G',
                           pop = 'pop'
)

# up to including flowering date
spline_flower_fac <- fitSplineHDM(inDat = statgen_flower_df,
                           trait = selected_factor,
                           genotype = 'G',
                           pop = 'pop'
)

# up to flowering date
spline_veg_fac <- fitSplineHDM(inDat = statgen_veg_df,
                                  trait = selected_factor,
                                  genotype = 'G',
                                  pop = 'pop'
)

# interpolate across timepoints for errors and spline estimates
pred_all.psHDM <- predict(object = spline_all_fac,
                      newtimes = seq(min(spline_all_fac$time[["timeNumber"]]),
                                     max(spline_all_fac$time[["timeNumber"]]),
                                     length.out = 100),
                      pred = list(pop = TRUE, geno = TRUE, plot = TRUE),
                      se = list(pop = FALSE, geno = FALSE, plot = FALSE),
                      trace = FALSE)

# up to flowering date
pred_flowering.psHDM <- predict(object = spline_flower_fac,
                      newtimes = seq(min(spline_flower_fac$time[["timeNumber"]]),
                                     max(spline_flower_fac$time[["timeNumber"]]),
                                     length.out = 100),
                      pred = list(pop = TRUE, geno = TRUE, plot = TRUE),
                      se = list(pop = FALSE, geno = FALSE, plot = FALSE),
                      trace = FALSE)

# up to veg date
pred_veg.psHDM <- predict(object = spline_veg_fac,
                                newtimes = seq(min(spline_veg_fac$time[["timeNumber"]]),
                                               max(spline_veg_fac$time[["timeNumber"]]),
                                               length.out = 100),
                                pred = list(pop = TRUE, geno = TRUE, plot = TRUE),
                                se = list(pop = FALSE, geno = FALSE, plot = FALSE),
                                trace = FALSE)

# store splines of factors ----
selected_factor
write_rds(pred_all.psHDM, paste0(paste('./data/processed/splines/all_data/predsplinedatall',
                            selected_factor,
                            'mid_timepoint',
                            mid_time,
                            selected_year,
                            selected_env, sep = '_'),'.RDS'))

selected_factor
write_rds(pred_flowering.psHDM, paste0(paste('./data/processed/splines/all_data/predsplinedatheading',
                                   selected_factor,
                                   'mid_timepoint',
                                   mid_time,
                                   selected_year,
                                   selected_env, sep = '_'),'.RDS'))
selected_factor
write_rds(pred_veg.psHDM, paste0(paste('./data/processed/splines/all_data/predsplinedatveg',
                                             selected_factor,
                                             'mid_timepoint',
                                             mid_time,
                                             selected_year,
                                             selected_env, sep = '_'),'.RDS'))

# LmmSolver with other penalty comparison -----
library(LMMsolver)

statgen_df$xF <- as.factor(statgen_df$timePoint)
## Convert time point to time number with the first time point as 0.
minTime <- min(statgen_df[["timePoint"]], na.rm = TRUE)
statgen_df[["timeNumber"]] <- as.numeric(statgen_df[["timePoint"]] - minTime) / 1000

nseg_10 <- LMMsolve(factor1 ~ 1,
                    spline = ~spl1D(timeNumber, nseg = 10),
                    random = ~xF + G,
                    data = statgen_df)
summary(nseg_10)

# 2. logL as function of number of segments
nseg <- 1:10 *100
nseg
logL <- rep(NA, length(nseg))
for (i in 1:length(nseg)) {
  obj <- LMMsolve(factor1 ~ 1,
                  spline = ~spl1D(timeNumber, nseg = i),
                  random = ~xF + G,
                  data = statgen_df)
  logL[i] <- obj$logL
}
plot(nseg,y=logL,type='b')

obj <- LMMsolve(factor1 ~ 1,
                spline = ~spl1D(timeNumber, nseg = 100),
                random = ~xF,
                data = statgen_df)
summary(obj)

obj <- LMMsolve(factor1 ~ 1,
                spline = ~spl1D(timeNumber, nseg = 100),
                random = ~xF + G,
                data = statgen_df)
summary(obj)
summary(spline_all_fac)
plot(spline_all_fac)
plot(obj$yhat, obj$residuals)
