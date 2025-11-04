custom_gfBLUP <- function(data, selection = NULL, K, sepExp = FALSE, verbose = TRUE){
  # Determining training and test set as well as the scenario (CV1 or CV2):
  train.set <- as.character(unique(data$G[!is.na(data$Y)]))
  test.set <- as.character(setdiff(unique(data$G), train.set))

  if (all(is.na(data[which(data$G %in% test.set), 2:(ncol(data) - 1)]))) {
    CV <- "CV1"
  } else if (!any(is.na(data[which(data$G %in% test.set), 2:(ncol(data) - 1)]))) {
    CV <- "CV2"
  } else {
    stop("Either all secondary traits must be present for the test set, or missing...")
  }

  # Number of replicates:
  n.rep.vector <- as.integer(table(data$G))
  reps <- (sum(n.rep.vector) - sum(n.rep.vector^2) / sum(n.rep.vector)) / (length(n.rep.vector) - 1)

  # Subsetting data to the selected factors if specified:
  if (!is.null(selection)) {
    data.subset <- data[, c("G", selection, "Y")]
  } else {
    selection <- names(data)[2:(length(names(data)) - 1)]
    data.subset <- data[, c("G", selection, "Y")]
  }

  # Determining covariances between factors and focal trait:
  if (verbose) {cat("Determining covariances between factors and focal trait...\n")}
  covmats <- gfBLUP::covSS(data.subset[which(data.subset$G %in% train.set),], sepExp = sepExp, verbose = TRUE)

  # Reducing pseudoCRD data to genotypic means:
  data.subset.means <- gfBLUP:::genotypeMeans(data.subset)
  data.subset.means <- data.subset.means[match(c(test.set, train.set), data.subset.means$G),]
  rownames(data.subset.means) <- c(test.set, train.set)

  # Calculating training set BLUPs:
  BLUPs_train <- tryCatch({
    BLUPs_train <- gfBLUP:::gBLUPo(Y = as.matrix(data.subset.means[train.set, -1]),
                                   K = K,
                                   Sg = covmats$Sg,
                                   Se = covmats$Se / reps,
                                   train.set = train.set)
  },error = function(e){
    print(e)
    print('Residual covariance matrix likely singular, nearPD used')
    BLUPs_train <- gfBLUP:::gBLUPo(Y = as.matrix(data.subset.means[train.set, -1]),
                                   K = K,
                                   Sg = covmats$Sg,
                                   Se = Matrix::nearPD(covmats$Se / reps, base.matrix = TRUE)$mat,
                                   train.set = train.set)
    return(BLUPs_train)
  })


  if (CV == "CV1") {
    # Calculating focal trait test set BLUPs:
    focalBLUP_test <- train2test(BLUPs_train$focal_gBLUP, K = K, train.set = train.set, test.set = test.set)

  } else if (CV == "CV2") {
    # Calculating all test set BLUPs:
    allBLUPs_test <- train2test(BLUPs_train$all_gBLUPs, K = K, train.set = train.set, test.set = test.set)

    BLUP_Y_test <- allBLUPs_test[, "Y"]
    BLUP_S_test <- allBLUPs_test[, setdiff(colnames(allBLUPs_test), "Y")]
    Ktt <- K[test.set, test.set]
    Kto <- K[test.set, train.set]
    Koo <- K[train.set, train.set]
    Vc <- kronecker(covmats$Sg[selection, selection], solve(Ktt)) + kronecker(covmats$Se[selection, selection] / reps, diag(dim(Ktt)[1]))

    # Final CV2 test set focal trait BLUP calculation following two-step approach (Runcie & Cheng, 2019):
    focalBLUP_test <- as.numeric(matrix(BLUP_Y_test) + kronecker(t(as.matrix(covmats$Sg[selection, "Y"])), solve(Ktt)) %*%
                                   solve(Vc) %*% (matrix(as.matrix(data.subset.means[test.set, selection])) - matrix(BLUP_S_test)))

    names(focalBLUP_test) <- test.set
  }

  h2s <- diag(covmats$Sg) / diag(covmats$Sg + covmats$Se)
  names(h2s) <- colnames(covmats$Sg)

  return(list(preds = focalBLUP_test,
              factorSelection = selection,
              Sg = covmats$Sg,
              Se = covmats$Se,
              h2s = h2s))
}
