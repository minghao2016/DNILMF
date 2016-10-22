
# see python code of NRLMF
# this log-likelihood for NRLMF
calcLogLik <- function(
  cc = 5,
  inMat,
  U,
  lamU = 0.125,
  V,
  lamV = 0.125,
  lapD,
  thisAlpha = 0.25,
  lapT,
  thisBeta = 0.125
  ) {
  
  # INPUT
  # cc: default 5
  # inMat: input interaction matrix
  # U: latent matrix for rows
  # lamU: lambda for U
  # V: latent matrix for cols
  # lamV: lambda for V
  # lapD: Laplacian matrix for drugs
  # thisAlpha: coefficient of lapD
  # lapT: Laplacian matrix for targets
  # thisBeta: coefficient of lapT
  
  # OUTPUT
  # a scalar of log-likelihood
    
  
  Y <- inMat
  cY <- cc * Y
  
  Ut <- t(U)
  Vt <- t(V)
  
  UVt <- U %*% Vt
  
  LL <- sum((1 + cY - Y) * log(1 + exp(UVt))) - sum(cY * UVt) + 
    0.5 * lamU * (base::norm(U, "F") ^ 2) + 0.5 * lamV * (base::norm(V, "F") ^ 2) +
    0.5 * thisAlpha * sum(base::diag(Ut %*% lapD %*% U)) + 0.5 * thisBeta * sum(base::diag(Vt %*% lapT %*% V)) # eq (12) of NRLMF
  LL <- -LL # ascend order
  
  return(LL)
}
