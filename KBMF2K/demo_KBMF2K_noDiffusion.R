
rm(list = ls())
setwd("C:/Users/kevin/Desktop/KBMF2k/")


# current data set name
db <- "nr"

switch (db,
        en = {
          cat("en data\n")
          flush.console()
          sd <- read.table("e_simmat_dc.txt")
          sd <- as.matrix(sd)
          st <- read.table("e_simmat_dg.txt")
          st <- as.matrix(st)
          Y <- read.table("e_admat_dgc.txt")
          Y <- as.matrix(Y) 
          Y <- t(Y)         
        },
        ic = {
          cat("ic data\n")
          flush.console()
          sd <- read.table("ic_simmat_dc.txt")
          sd <- as.matrix(sd)
          st <- read.table("ic_simmat_dg.txt")
          st <- as.matrix(st)
          Y <- read.table("ic_admat_dgc.txt")
          Y <- as.matrix(Y)
          Y <- t(Y)
        },
        gpcr = {
          cat("gpcr data\n")
          flush.console()
          sd <- read.table("gpcr_simmat_dc.txt")
          sd <- as.matrix(sd)
          st <- read.table("gpcr_simmat_dg.txt")
          st <- as.matrix(st)
          Y <- read.table("gpcr_admat_dgc.txt")
          Y <- as.matrix(Y)
          Y <- t(Y)
        },
        nr = {
          cat("nr data\n")
          flush.console()
          sd <- read.table("nr_simmat_dc.txt")
          sd <- as.matrix(sd)
          st <- read.table("nr_simmat_dg.txt")
          st <- as.matrix(st)
          Y <- read.table("nr_admat_dgc.txt")
          Y <- as.matrix(Y)
          Y <- t(Y)
        },
        srep = {
          cat("Scientific Reports data\n")
          flush.console()
          sd <- read.table("drug ChemSimilarity.txt", sep = ",")
          sd <- data.matrix(sd)
          st <- read.table("target SeqSimilarity.txt", sep = ",")
          st <- data.matrix(st)
          Y <- read.table("adjacent Matrix.txt", sep = ",")
          Y <- data.matrix(Y)
        },
        stop("db should be one of the follows: 
             {en, ic, gpcr, nr or srep}\n")
        )


numDrug <- nrow(sd)
numTarget <- nrow(st)

## load required packages
pkgs <- c("matrixcalc", "data.table", "Rcpp", "ROCR", "Bolstad2", "MESS")
rPkgs <- lapply(pkgs, require, character.only = TRUE)

## source required R files
rSourceNames <- c("doCrossValidation.R",
                  "constrNeig.R", 
                  "inferZeros.R",
                  "calcAUPR.R",
                  "kbmf_regression_train.R",
                  "kbmf_regression_test.R",
                  "kbmf1k1k_supervised_regression_variational_train.R",
                  "kbmf1k1k_supervised_regression_variational_test.R")
rSN <- lapply(rSourceNames, source, verbose = FALSE)

## sourceCPP required C++ files
cppSourceNames <- c("fastKF.cpp", "fastKgipMat.cpp")
cppSN <- lapply(cppSourceNames, sourceCpp, verbose = FALSE)


## convert to kernel
isKernel <- TRUE
if (isKernel) {
  if (!isSymmetric(sd)) {
    sd <- (sd + t(sd)) / 2
  }
  epsilon <- 0.1
  while (!is.positive.semi.definite(sd)) {
    sd <- sd + epsilon * diag(nrow(sd))
  }
  if (!isSymmetric(st)) {
    st <- (st + t(st)) / 2
  }
  epsilon <- 0.1
  while (!is.positive.semi.definite(st)) {
    st <- st + epsilon * diag(nrow(st))
  }
}

## do cross-validation
kfold <- 10
numSplit <- 5
K1 <- 5

## split training and test sets
savedFolds <- doCrossValidation(Y, kfold = kfold, numSplit = numSplit)

## hyper-parameters

numLat <- 20
runDiffusion <- FALSE

## for saving results
AUPRVec <- vector(length = kfold)
AUCVec <- vector(length = kfold)
finalResult <- matrix(NA, nrow = numSplit, ncol = 2)
colnames(finalResult) <- c("AUPR", "AUC")

# main loop
for (i in 1:numSplit) {
  for (j in 1:kfold) {
    cat("numSplit:", i, "/", numSplit, ";", "kfold:", j, "/", kfold, "\n")
    flush.console()
    Y <- savedFolds[[i]][[j]][[7]]
    if (runDiffusion) {
      Yr <- inferZeros(Y, sd, K = K1)
      Yc <- inferZeros(t(Y), st, K = K1)
      KgipD <- fastKgipMat(Yr, 1)
      KgipT <- fastKgipMat(Yc, 1)
      # nNeig = 3, nIter = 2
      sd <- fastKF(KgipD, sd, 3, 2)
      st <- fastKF(KgipT, st, 3, 2)
      lap <- constrNeig(sd, st, K = K1)
      lapD <- lap$lapD
      lapT <- lap$lapT
      simD <- lap$simD
      simT <- lap$simT
    }
    ## KBMF2K 
    Kx <- array(0, c(numDrug, numDrug, 1))
    Kx[, , 1] <- sd
    Kz <- array(0, c(numTarget, numTarget, 1))
    Kz[, , 1] <- st
    state <- kbmf_regression_train(Kx = Kx, Kz = Kz, Y = Y, R = numLat)
    prediction <- kbmf_regression_test(Kx = Kx, Kz = Kz, state = state)
    score <- prediction$Y$mu
    ####################################################################
    
    testIndexRow = savedFolds[[i]][[j]][[3]]      
    testIndexCol = savedFolds[[i]][[j]][[4]]      
    testLabel = savedFolds[[i]][[j]][[1]]         
    # test set
    testSetIndex <- cbind(testIndexRow, testIndexCol)
    score <- score[testSetIndex]
    result <- calcAUPR(testLabel, score)
    AUPRVec[j] <- result[1, "aupr"]
    AUCVec[j] <- result[1, "auc"]
  }
  AUPR <- mean(AUPRVec)
  AUC <- mean(AUCVec)
  finalResult[i, "AUPR"] <- AUPR
  finalResult[i, "AUC"] <- AUC
}


## print the result
cat(
  "\n======================\n\n",
  "db is: ", db, "\n",
  ## hyper-parameters
  "numLat = ", numLat, "\n",
  "runDiffusion = ", runDiffusion, "\n",
  "\n=====================\n")
cat(numSplit, "trails 10-fold CV", "\n")
cat("KBMF2K: no diffusion kernels \n")
cat("aupr:", mean(finalResult[, "AUPR"]), "+/-", sd(finalResult[, "AUPR"]), "\n")
cat("auc:", mean(finalResult[, "AUC"]), "+/-", sd(finalResult[, "AUC"]), "\n")

# save to file
curDate <- format(Sys.time(), format = "%Y-%m-%d")
curTime <- format(Sys.time(), format =  "%H.%M.%S")
savedFileName <- paste0(db, "-", "KBMF2K", "_", "runDiffusion", "-", runDiffusion, "-", curDate, "_", curTime, ".RData")
cat("\n\n")
print(savedFileName)
save.image(file = savedFileName)