# CIMMYT proc factor splines example with statgenHTP

# libraries
library(tidyverse)
library(statgenHTP)
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

  eigen_g <- eigen(cor_g)

  # cor_g <- RF(cor_g)
  HS_factor <- factorModel(df, cormat = cor_g, what = 'genetic', m = m)
  return(list(factor_model = HS_factor, eigen_g = eigen_g))
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
load('./data/processed/hyper_folds/B5IR/hyper_dataset_1.RData')
HS_B5IR <- x
load('./data/processed/hyper_folds/B2IR/hyper_dataset_1.RData')
HS_B2IR <- x
load('./data/processed/hyper_folds/DRIP/hyper_dataset_1.RData')
HS_DRIP <- x
rm(x)
## HS wrangling ----
# description of environment identifiers
# Bed2IR: Moderate Drought, two furrow irrigations
# Melga: Optimal Bed, five furrow irrigations
# Bed5IR: Optimal Flat, Five Flood irrigations
# Drip: Severe Drought, Three minimal irrigations through drip
# Heat: Heat, late february planting five furrow irrigations
temp_df <- HS_DRIP$data

# factor creation only requires secondary traits to only contain g + traits, no other columns
HS_list <- list(`Severe Drought` = select(HS_DRIP$data, 'G', contains('nm')),
                `Moderate Drought` = select(HS_B2IR$data, 'G', contains('nm')),
                `Optimal Flat` = select(HS_B5IR$data, 'G', contains('nm'))
)


procrustes_loadings_generator <- function(df, focal_timepoint){
  ## wrangling for different timepoint factors ----
  # get dates
  timepoints <- unique(str_split(colnames(df)[-1], pattern = '_', simplify = TRUE)[,2])

  # set sepecific timepoint
  # Optimal Flat/Bed: 25th of February is heading date
  # Severe/Moderate Drought: 9th of February is heading date
  timepoints[focal_timepoint]

  # first fit factor for best timepoint (here heading date 25th of february)
  factor_target <- generate_factor(timepoint = timepoints[focal_timepoint], df = df, m = 2)

  # target matrix has m factors, so we generate m factors for the other timepoints as well
  factors <- lapply(timepoints, df = df, FUN = function(x, df, m = factor_target$factor_model$m){
    df <- select(df, G, contains(x))
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
  proc_factors <- lapply(factors[-focal_timepoint], function(x){
    procustes(x$loadings, factor_target$factor_model$loadings, type = 'orthogonal', verbose = TRUE)
  })

  # Procrustean smoothing ----
  rotmat_origin <- proc_factors[[focal_timepoint]]$rotmatrix

  smoothed <- lapply(proc_factors, function(x){
    procustes_smoother(x$rotmatrix)
  })

  # last timepoint got columns switched
  smoothed_factors <- map2(factors[-focal_timepoint], smoothed, function(factor_object, rotation_matrix){
    factor_object$loadings %*% rotation_matrix
  })

  origin_mat <- as.matrix(factors[[focal_timepoint]]$loadings[,])
  colnames(origin_mat) <- gsub(x = colnames(origin_mat), pattern = 'F', replacement =  'V')

  smoothed_factors <- append(smoothed_factors, list(origin_mat), after = focal_timepoint-1)

  # create better timepoints for labels
  timepoints <-  gsub("(\\d{2})(\\d{2})(\\d{2})$","20\\1-\\2-\\3", timepoints)
  names(smoothed_factors) <- timepoints

  # create overal df ----
  smoothed_factor_df <- lapply(smoothed_factors, as.data.frame) %>%
    bind_rows(.id = 'timePoint') %>%
    rownames_to_column('trait') %>%
    separate_wider_regex(trait, patterns = c('nm', wavelength = '[0-9]+', '_', '[0-9]+')) %>%
    rename_with(~ gsub('V', 'Factor_', .x)) %>%
    pivot_longer(contains('Factor'), names_to = 'factor', values_to = 'loading') %>%
    separate(factor, into = c(NA, 'factor')) %>%
    mutate(wavelength = as.numeric(wavelength))

  return(smoothed_factor_df)
}

varimax_loadings_generator <- function(df){
  ## wrangling for different timepoint factors ----
  # get dates
  timepoints <- unique(str_split(colnames(df)[-1], pattern = '_', simplify = TRUE)[,2])

  # set sepecific timepoint
  # Optimal Flat/Bed: 25th of February is heading date
  # Severe/Moderate Drought: 9th of February is heading date

  # target matrix has m factors, so we generate m factors for the other timepoints as well
  factors <- lapply(timepoints, df = df, FUN = function(x, df, m = 2){
    df <- select(df, G, contains(x))
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

  # wrap up varimax loadings result
  vari_load <- lapply(factors, function(x){data.frame(trait = rownames(x$loadings),
                                                      F1 = x$loadings[, 1],
                                                      F2 = x$loadings[, 2])})
  # create better timepoints for labels
  timepoints <-  gsub("(\\d{2})(\\d{2})(\\d{2})$","20\\1-\\2-\\3", timepoints)

  # add timepoints labels
  names(vari_load) <- timepoints

  vari_df <- bind_rows(vari_load, .id = 'timePoint') %>%
    separate_wider_regex(trait, patterns = c('nm', wavelength = '[0-9]+', '_', '[0-9]+')) %>%
    pivot_longer(contains('F'), names_to = 'factor', values_to = 'loading') %>%
    separate_wider_regex(factor, patterns = c('NA' = 'F', factor = '[0-9]+')) %>%
    mutate(wavelength = as.numeric(wavelength)) %>%
    select(timePoint, wavelength, factor, loading)

  return(vari_df)
}

HS_proc_loadings <- mapply(procrustes_loadings_generator, HS_list, c(3, 3, 4), SIMPLIFY = FALSE)
HS_vari_loadings <- lapply(HS_list, varimax_loadings_generator)

# visualisation of factors ----
# single environment plot
HS_vari_loadings$`Severe Drought` %>%
  ggplot(aes(x = wavelength, y = loading, colour = factor)) +
  geom_point() +
  geom_line() +
  facet_wrap(~timePoint) +
  theme_classic()



# multi environment plot ----
## varimax rotated ----
(loadings_plot_vari <- bind_rows(HS_vari_loadings, .id = 'environment') %>%
   mutate(Factor = as.factor(factor),
          Treatment = factor(environment, levels = c('Optimal Flat', 'Moderate Drought', 'Severe Drought')),
          `Factor Loading` = loading,
          Wavelength = wavelength,
          new_fact = fct_cross(environment, factor, sep = ' Factor '),
          .keep = 'unused') %>%
   ggplot(aes(x = Wavelength, y = `Factor Loading`, colour = Factor, shape = Treatment)) +
   geom_point(size = 0.8) +
   # geom_line() +
   facet_wrap(~timePoint) +
   theme_classic())


ggsave('./figures/loadings_plot_varimax.pdf', plot = loadings_plot_vari, device = pdf)

## procrustes rotated ----
(loadings_plot_proc <- bind_rows(HS_proc_loadings, .id = 'environment') %>%
  mutate(Factor = as.factor(factor),
         Treatment = factor(environment, levels = c('Optimal Flat', 'Moderate Drought', 'Severe Drought')),
         `Factor Loading` = loading,
         Wavelength = wavelength,
         new_fact = fct_cross(environment, factor, sep = ' Factor '),
         .keep = 'unused') %>%
  ggplot(aes(x = Wavelength, y = `Factor Loading`, colour = Factor, shape = Treatment)) +
  geom_point(size = 0.8) +
  # geom_line() +
  facet_wrap(~timePoint) +
  theme_classic())

ggsave('./figures/loadings_plot_procrustes.pdf', plot = loadings_plot_proc, device = 'pdf')

# presentation figures ----
## timepoint 2015-03-23 for example of rotation
(loadings_20150323 <- bind_rows(HS_vari_loadings, .id = 'environment') %>%
   mutate(Factor = as.factor(factor),
          Treatment = factor(environment, levels = c('Optimal Flat', 'Moderate Drought', 'Severe Drought')),
          `Factor Loading` = loading,
          Wavelength = wavelength,
          new_fact = fct_cross(environment, factor, sep = ' Factor '),
          .keep = 'unused') %>%
   filter(timePoint %in% c('2015-03-15', '2015-03-23'))  |>
   ggplot(aes(x = Wavelength, y = `Factor Loading`, colour = Factor, shape = Treatment)) +
   geom_point() +
   # geom_line() +
   facet_wrap(~timePoint) +
   theme_classic())


(loadings_plot_proc <- bind_rows(HS_proc_loadings, .id = 'environment') %>%
    mutate(Factor = as.factor(factor),
           Treatment = factor(environment, levels = c('Optimal Flat', 'Moderate Drought', 'Severe Drought')),
           `Factor Loading` = loading,
           Wavelength = wavelength,
           new_fact = fct_cross(environment, factor, sep = ' Factor '),
           .keep = 'unused') %>%
    filter(timePoint %in% c('2015-03-15', '2015-03-23')) |>
    ggplot(aes(x = Wavelength, y = `Factor Loading`, colour = Factor, shape = Treatment)) +
    geom_point() +
    # geom_line() +
    facet_wrap(~timePoint) +
    theme_classic())

(loadings_plot_proc <- bind_rows(HS_proc_loadings, .id = 'environment') %>%
    mutate(Factor = as.factor(factor),
           Treatment = factor(environment, levels = c('Optimal Flat', 'Moderate Drought', 'Severe Drought')),
           `Factor Loading` = loading,
           Wavelength = wavelength,
           new_fact = fct_cross(environment, factor, sep = ' Factor '),
           .keep = 'unused') %>%
    filter(timePoint %in% '2015-03-15' & Treatment == 'Optimal Flat') |>
    ggplot(aes(x = Wavelength, y = `Factor Loading`, colour = Factor, shape = Treatment)) +
    geom_point(size = 3) +
    # geom_line() +
    facet_wrap(~timePoint) +
    theme_classic())

