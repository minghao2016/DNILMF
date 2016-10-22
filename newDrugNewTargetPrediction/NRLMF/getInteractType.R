getInteractType <- function(Y) {
  # INPUT
  # Y: matrix for input interaction matrix with {0, 1}
  # row refers to drugs and column refers to targets
  
  # OUTPUT
  # typeComb: data.table for interaction type combination
  

  nr <- nrow(Y)
  nc <- ncol(Y)
  
  rs <- rowSums(Y) # get sum by row
  cs <- colSums(Y) # get sum by col
  
  # find new drug, known drug, new target, known target
  rr <- rep("knownDrug", nr)     # initialize rows
  cc <- rep("knownTarget", nc)   # initialize cols
  
  idxZeroRows <- which(rs == 0)
  if (length(idxZeroRows) > 0) {
    rr[idxZeroRows] <- "newDrug"
  }
  
  idxZeroCols <- which(cs == 0)
  if (length(idxZeroCols) > 0) {
    cc[idxZeroCols] <- "newTarget"
  }
  
  # stack by column in R
  rr <- rep(rr, times = nc) # stack row
  cc <- rep(cc, each = nr)  # stack col
  
  # paste together
  interactType <- paste(rr, cc, sep = "-")
  
  # dictionary
  library(data.table)
  dict <- data.table(
    typeName = c(
        "newDrug-newTarget",
        "newDrug-knownTarget",
        "knownDrug-newTarget",
        "knownDrug-knownTarget"
      ),
      typeId = c(1, 2, 3, 4)
    )
                     
  # matched positions
  matchedPos <- match(interactType, dict[, typeName])
  typeComb <- dict[matchedPos, ]
  
  # cbind
  typeComb <- cbind(interactType, typeComb)

  # transform Y to i-j-v format
  rownames(Y) <- NULL
  colnames(Y) <- NULL
  
  # width -> length, stack by column using "melt()"
  # i, j, v format
  triplet <- as.data.table(reshape2::melt(Y))
  setnames(triplet, c("ii", "jj", "vv"))
  
  # cbind agaion
  typeComb <- cbind(triplet, typeComb)
  
  return(typeComb)
}

