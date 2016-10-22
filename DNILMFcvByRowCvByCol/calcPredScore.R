calcPredScore <- function(U,
                          V,
                          simDrug,
                          simTarget,
                          knownDrugIndex,
                          knownTargetIndex,
                          testIndexRow,
                          testIndexCol,
                          K = 5,
                          testLabel,
                          thisAlpha,
                          thisBeta,
                          thisGamma) {
  # INPUT
  # U: row latent matrix
  # V: col latent matrix
  # simDrug: similarity matrix for drug, but diagonal elements are zeros
  # simTarget: similarity matrix for target, but diagonal elements are zeros
  # testIndexRow: row index for test set
  # testIndexCol: col index for test set
  # K: number of nearest neighbor for prediction
  # testLabel: labels for the test set
  
  # OUTPUT 
  # a list of AUC and AUPR
  
  if (K < 0) {
    stop("K MUST be '>=' 0! \n")
  }
  
  if (K > 0) {
    ## cat("with K smoothing! \n")
    ## for drug
    indexTestD <- unique(testIndexRow)
    testD <- U[indexTestD, ]
    testD <- cbind(indexTestD, testD)
    numTest <- length(indexTestD)
    numColTestD <- ncol(testD)
    simDrugKnown <- simDrug[, knownDrugIndex]
    numDrugKnown <- length(knownDrugIndex)
    
    for (i in 1:numTest) {
      indexCurr <- indexTestD[i]
      isNewDrug <- !(indexCurr %in% knownDrugIndex)
      if (isNewDrug) {
        simDrugNew <- simDrugKnown[indexCurr, ] # vector
        indexRank <- rank(simDrugNew) # vector
        indexNeig <- which(indexRank > (numDrugKnown - K))
        simCurr <- simDrugNew[indexNeig] # vector
        # index for U
        index4U <- knownDrugIndex[indexNeig]
        U_Known <- U[index4U, , drop = FALSE] # force to matrix
        # vec %*% matrix => matrix
        testD[i, 2:numColTestD] <- (simCurr %*% U_Known) / sum(simCurr)
      }
    }
    
    Unew <- U
    Unew[indexTestD, ] <- testD[, -1]
    
    ## for target
    # unique index for test target
    indexTestT <- unique(testIndexCol)
    testT <- V[indexTestT, ]
    # add first column as labels
    testT <- cbind(indexTestT, testT) # 1st column is unique test label
    # number of unique test set
    numTest <- length(indexTestT)
    # number of column for testT
    numColTestT <- ncol(testT)
    # known similarity matrix for targets
    simTargetKnown <- simTarget[, knownTargetIndex]
    # number of known targets
    numTargetKnown <- length(knownTargetIndex)
    
    for (i in 1:numTest) {
      indexCurr <- indexTestT[i]
      isNewTarget <- !(indexCurr %in% knownTargetIndex)
      if (isNewTarget) {
        simTargetNew <- simTargetKnown[indexCurr, ] # vector
        indexRank <- rank(simTargetNew) # vector
        # selected neighbor index with top K neighbor
        indexNeig <- which(indexRank > (numTargetKnown - K))
        # get similarity value of K
        simCurr <- simTargetNew[indexNeig] # vector
        # index for V
        index4V <- knownTargetIndex[indexNeig]
        V_Known <- V[index4V, , drop = FALSE] # force to matrix
        # vec %*% matrix => matrix
        testT[i, 2:numColTestT] <- (simCurr %*% V_Known) / sum(simCurr)
      }
    }

    Vnew <- V
    Vnew[indexTestT, ] <- testT[, -1]
    
    Vnewt <- t(Vnew)
    UnewVnewt <- Unew %*% Vnewt

    val <- thisAlpha * UnewVnewt + thisBeta * (simDrug %*% UnewVnewt) + thisGamma * (UnewVnewt %*% simTarget) 
    testSetIndex <- cbind(testIndexRow, testIndexCol)
    val <- val[testSetIndex]
    
    # score from val
    score <- exp(val) / (1 + exp(val))
    result <- calcAUPR(testLabel, score)
  } else {  # K = 0 condition
    # cat("without K smoothing! \n")
    # flush.console()
    Vt <- t(V)
    UVt <- U %*% Vt
    val <- thisAlpha * UVt + thisBeta * (simDrug %*% UVt) + thisGamma * (UVt %*% simTarget) 
    testSetIndex <- cbind(testIndexRow, testIndexCol)
    val <- val[testSetIndex]
    # score
    score <- exp(val) / (1 + exp(val))
    result <- calcAUPR(testLabel, score)
  }
  return(result)
}
