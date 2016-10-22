updateUV <- function(cc = 5,
                     inMat,
                     numLat = 100,
                     lamU = 0.125,
                     lamV = 0.125,
                     lapD,
                     thisAlpha = 0.25,
                     lapT,
                     thisBeta = 0.125,
                     thisTheta = 0.5,
                     initMethod = "useNorm",
                     thisSeed = 7771,
                     maxIter = 100) {
  
  
  # INPUT
  # cc:
  # inMat:
  # numLat:
  # lamU:
  # lamV:
  # lapD
  # thisAlpha:
  # lapT
  # thisBeta:
  # thisTheta:
  # initMethod:
  # thisSeed:
  # maxIter:
  
  # OUTPUT:
  # a list with two elements: U and V
  
  
  numRow <- nrow(inMat)
  numCol <- ncol(inMat)
  
  if (initMethod == "useNorm") {
    U <- matrix(NA, nrow = numRow, ncol = numLat)
    U <- apply(U, 2, function(x) {
      rnorm(x, mean = 0, sd = 1)
    })
    U <- base::sqrt(1 / numLat) * U
    V <- matrix(NA, nrow = numCol, ncol = numLat)
    V <- apply(V, 2, function(x) {
      rnorm(x, mean = 0, sd = 1)
    })
    V <- base::sqrt(1 / numLat) * V
  } else if (initMethod == "useSeed") {
    set.seed(thisSeed)
    U <- matrix(NA, nrow = numRow, ncol = numLat)
    U <- apply(U, 2, function(x) {
      rnorm(x, mean = 0, sd = 1)
    })
    U <- base::sqrt(1 / numLat) * U 
    V <- matrix(NA, nrow = numCol, ncol = numLat)
    V <- apply(V, 2, function(x) {
      rnorm(x, mean = 0, sd = 1)
    })
    V <- base::sqrt(1 / numLat) * V 
  } else {
    stop("initMethod should be one of {useNorm, useSeed}\n")  
  }
  
  sumGradU <- matrix(0, nrow = numRow, ncol = numLat)
  sumGradV <- matrix(0, nrow = numCol, ncol = numLat)

  # last log-likelihood
  lastLog <- calcLogLik(
      cc = cc,
      inMat = inMat,
      U = U,
      lamU = lamU,
      V = V,
      lamV = lamV,
      lapD = lapD,
      thisAlpha = thisAlpha,
      lapT = lapT,
      thisBeta = thisBeta)
  
  currDeltaLL <- 1000
  # main loop
  for (i in 1:maxIter) {
    # gradU
    gradU <- calcDeriv(
      cc = cc,
      inMat = inMat,
      U = U,
      lamU = lamU,
      V = V,
      lamV = lamV,
      lapD = lapD,
      thisAlpha = thisAlpha,
      lapT = lapT,
      thisBeta = thisBeta,
      isGradU = TRUE)
    sumGradU <- sumGradU + (gradU ^ 2)
    stepSize <- thisTheta / sqrt(sumGradU)
    U <- U + stepSize * gradU
    
    # gradV
    gradV <- calcDeriv(
      cc = cc,
      inMat = inMat,
      U = U,
      lamU = lamU,
      V = V,
      lamV = lamV,
      lapD = lapD,
      thisAlpha = thisAlpha,
      lapT = lapT,
      thisBeta = thisBeta,
      isGradU = FALSE
    )
    sumGradV <- sumGradV + (gradV ^ 2)
    stepSize <- thisTheta / sqrt(sumGradV)
    V <- V + stepSize * gradV
    
    currLog <- calcLogLik(
      cc = cc,
      inMat = inMat,
      U = U,
      lamU = lamU,
      V = V,
      lamV = lamV,
      lapD = lapD,
      thisAlpha = thisAlpha,
      lapT = lapT,
      thisBeta = thisBeta)
    
    # delta log-likelihood
    deltaLog <- (currLog - lastLog) / abs(lastLog)
    
    # stop earlier
    if (abs(deltaLog) < 1e-5) {
      break
    }
    
    if ((i > 50) & (deltaLog > currDeltaLL)) {
      break
    }
    
    currDeltaLL <- deltaLog
    lastLog <- currLog
  }
  
  UV <- list(U = U, V = V)
  return(UV)
}
