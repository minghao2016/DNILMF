## The following code is used for reproducing
## the results in the following paper:
## Predicting drug-target interactions by dual-network
## integrated logistica matrix factorization

## rm(list = ls())
## setwd("YourDir\\DNILMF")


if (TRUE) {stop("Please do not use this demo file due to information leak in the loop {sd, st}!", "\n")}



# current data set name
db <- "gpcr"

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

## load required packages
pkgs <- c("matrixcalc", "data.table", "Rcpp", "ROCR", "Bolstad2", "MESS")
rPkgs <- lapply(pkgs, require, character.only = TRUE)

## source required R files
rSourceNames <- c("doCrossValidation.R",
                  "constrNeig.R", 
                  "inferZeros.R",
                  "calcLogLik.R",
                  "calcDeriv.R",
                  "updateUV.R",
                  "calcPredScore.R",
                  "calAUPR.R",
                  "calcPredScore.R")
rSN <- lapply(rSourceNames, source, verbose = FALSE)

## sourceCPP required C++ files
cppSourceNames <- c("fastKF.cpp", "fastKgipMat.cpp", "log1pexp.cpp", "sigmoid.cpp")
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

## split training and test sets
savedFolds <- doCrossValidation(Y, kfold = kfold, numSplit = numSplit)

## hyper-parameters

isDefaultPara <- TRUE
if (isDefaultPara) {
  numLat <- 50
  cc <- 5
  thisAlpha <- 0.5
  lamU <- 5
  lamV <- 1
  K1 <- 5
} else {
  # best for gpcr data
  numLat <- 90
  cc <- 6
  thisAlpha <- 0.4
  lamU <- 2
  lamV <- 2
  K1 <- 2
}


## values according to hyper-parameters
thisBeta <- (1 - thisAlpha)/2
thisGamma <- 1 - thisAlpha - thisBeta

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
    Yr <- inferZeros(Y, sd, K = K1)
    Yc <- inferZeros(t(Y), st, K = K1)
    KgipD <- fastKgipMat(Yr, 1)
    KgipT <- fastKgipMat(Yc, 1)
    
    #########################################################################
    if (TRUE) {stop("Below sd and st are causing information leak!", "\n")}
    # nNeig = 3, nIter = 2
    sd <- fastKF(KgipD, sd, 3, 2)
    st <- fastKF(KgipT, st, 3, 2)
    lap <- constrNeig(sd, st, K = K1)
    if (TRUE) {stop("Above sd and st are causing information leak!", "\n")}
    #########################################################################
    
    lapD <- lap$lapD
    lapT <- lap$lapT
    simD <- lap$simD
    simT <- lap$simT
    ## use AdaGrid to update U and V
    UV <- updateUV(
      cc = cc,
      inMat = Y,
      thisAlpha = thisAlpha,
      thisBeta = thisBeta,
      Sd = simD,
      thisGamma = thisGamma,
      St = simT,
      lamU = lamU,
      lamV = lamV,
      numLat = numLat,
      initMethod = "useNorm",
      thisSeed = 123,
      maxIter = 100)
    
    U <- UV$U
    V <- UV$V
    
    knownDrugIndex <- savedFolds[[i]][[j]][[5]]
    knownTargetIndex <- savedFolds[[i]][[j]][[6]] 
    testIndexRow = savedFolds[[i]][[j]][[3]]      
    testIndexCol = savedFolds[[i]][[j]][[4]]      
    testLabel = savedFolds[[i]][[j]][[1]]         
    # result
    result <- calcPredScore(
      U = U,
      V = V,
      simDrug = simD,
      simTarget = simT,
      knownDrugIndex = knownDrugIndex,
      knownTargetIndex = knownTargetIndex,
      testIndexRow = testIndexRow,
      testIndexCol = testIndexCol,
      K = K1,
      testLabel = testLabel,
      thisAlpha = thisAlpha, 
      thisBeta = thisBeta,   
      thisGamma = thisGamma  
    )
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
  "cc = ", cc, "\n",
  "thisAlpha = ", thisAlpha, "\n",
  "lamU = ", lamU, "\n",
  "lamV = ", lamV, "\n",
  "K1 = ", K1, "\n",
  "\n=====================\n")
cat(numSplit, "trails 10-fold CV", "\n")
print(summary(finalResult))
cat("\n\n mean values:\n")
print(apply(finalResult, 2, mean))
cat("\n\n sd values:\n")
print(apply(finalResult, 2, sd))

# save to file
curDate <- format(Sys.time(), format = "%Y-%m-%d")
curTime <- format(Sys.time(), format =  "%H.%M.%S")
savedFileName <- paste0(db, "_", curDate, "_", curTime, ".RData")
cat("\n\n")
print(savedFileName)
#save.image(file = savedFileName)
