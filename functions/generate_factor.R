generate_factor <- function(df, timepoint = NULL, m = NULL){
  df <- dplyr::select(df, G, contains(timepoint))
  sumsq <- gfBLUP::covSS(df)

  data.means <- stats::aggregate(df[, 2:ncol(df)], list(df[, 1]), mean)
  names(data.means)[1] <- names(df)[1]
  data.means$G <- as.character(data.means$G)

  # create factors
  cor_g <- cov2cor(sumsq$Sg)

  if(is.null(m)){
    scree_analysis <- nFactors::nScree(x = eigen(cor_g)$values)
    m <- scree_analysis$Components$naf
  }

  factor_structure <- gfBLUP::factorModel(df, cormat = cor_g, what = 'genetic', m = m)

  # factor score
  factor_scores <- gfBLUP::factorScores(df,
                                        loadings = factor_structure$loadings,
                                        uniquenesses = factor_structure$uniquenesses,
                                        m = factor_structure$m)

  factor_result <- list(factor_structure = factor_structure,
                        factor_scores = factor_scores)
  return(factor_result)
}
