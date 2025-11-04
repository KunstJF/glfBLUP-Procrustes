# gblup
# libraries ----
library(doParallel)
library(tictoc)
library(dplyr)
library(magrittr)
library(rrBLUP)

# Loading functions:
for (file.name in list.files("functions")) {
  source(paste0("functions/", file.name))
}

# for each environment we need to perform gBLUP on the first 100 datasets
n.cores <- parallel::detectCores() - 2
cl <- parallel::makeCluster(n.cores, outfile = "logs/hyper_data_generation.txt")
doParallel::registerDoParallel(cl)

environments <- c('B5IR', 'B2IR', 'DRIP')
load('./data/processed/kinship.RData')

n_datasets <- 100
selected_env = environments[1]
i = 1
tic()
invisible(
  foreach::foreach(i = 1:n_datasets, .packages = c("magrittr", 'dplyr', 'rrBLUP', 'rlist')) %dopar% {
    load(sprintf('./data/processed/hyper_folds/%s/hyper_dataset_%d.Rdata', selected_env, i))

    df <- x$data %>%
      select(G, Y) %>%
      distinct()

    kinship <- kinship[rownames(kinship) %in% df$G, colnames(kinship) %in% df$G]
    gblup_result <- kin.blup(df, geno = 'G', pheno = 'Y', K = kinship)
    gblup_df <- data.frame(G = names(gblup_result$pred),
                           Y_pred = as.numeric(gblup_result$pred))

    x$pred.target %>%
      left_join(gblup_df, by = 'G') %$%
      plot(pred.target, Y_pred)

    gblup_result$rho <-   x$pred.target %>%
      left_join(gblup_df, by = 'G') %$%
      cor(pred.target, Y_pred)

    list.save(gblup_result, file = sprintf('./data/processed/gblup/%s/gblup_%d.RData', selected_env, i))
  })
toc()

