# look at Python code of NRLMF
# this derivative is from NRLMF, see eq (13), but ascend order
calcDeriv <- function(
  cc = 5,
  inMat,
  U,
  lamU = 0.125,
  V,
  lamV = 0.125,
  lapD,
  thisAlpha = 0.25,
  lapT,
  thisBeta = 0.125,
  isGradU = TRUE) {
  
  # INPUT
  # cc: augment fold for positive samples
  # inMat: input interaction matrix
  # U: latent matrix for rows
  # lamU: lambda for U
  # V: latent matrix for cols
  # lamV: lambda for V
  # lapD: Laplacian matrix for drugs
  # thisAlpha: coefficient of lapD
  # lapT: Laplacian matrix for target
  # thisBeta: weight coefficient of lapT
  # isGradU: TRUE or FALSE
  
  # OUTPUT
  # gradient of U or V
  
  UVt <- U %*% t(V)
  AA <- exp(UVt)
  P <- AA / (1 + AA)
  
  Y <- inMat
  
  if (isGradU) {
    gradU <- P %*% V + ((cc - 1) * (Y * P)) %*% V - (cc * Y) %*% V + lamU * U + (thisAlpha * lapD) %*% U # eq (13) of NRLMF
    gradU <- -gradU # ascend
    # return
    return(gradU)
  } else {
    Pt <- t(P)
    Yt <- t(Y)
    gradV <- Pt %*% U + ((cc - 1) * (Yt * Pt)) %*% U - (cc * Yt) %*% U + lamV * V + (thisBeta * lapT) %*% V # eq (13) of NRLMF
    gradV <- -gradV # ascend
    # return
    return(gradV)
  }
}