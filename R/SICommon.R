# C. File containing function used several time for different part of seqimpute
#
#
###############################################################################
# Consider time-dependent covariates

COtselection <- function(COtselected, COt, ncot, t1, t2, nr, nc, j, shifted) {
  COttemp <- as.data.frame(matrix(nrow = nr, ncol = 0))
  for (d in 1:(ncot / nc)) {
    COttemp <- cbind(COttemp, COt[, (j + shifted) + (d - 1) * nc])
  }
  COtselected[t1:t2, ] <- COttemp

  return(COtselected)
}


#############################################################################
# Concatenating CDi with COtselected_i (the matrix containing the current
# time-dependent covariates) Checking if COt is NOT completely empty

COtselectionSpe <- function(CDi, COt, ncot, nc, i, j, k) {
  if (ncot > 0) {
    COtselected_i <- as.data.frame(matrix(nrow = 1, ncol = 0))
    for (d in 1:(ncot / nc)) {
      COtselected_i <- cbind(COtselected_i, COt[i, (j) + (d - 1) * nc])
    }
    COtselected_i <- do.call(rbind, replicate(k, 
      as.data.frame(COtselected_i[1, ]), simplify = FALSE))
    # Concatenating CDi and COtselected_i
    # into CDi
    CDi <- cbind(CDi, COtselected_i)
    # Transformation of the names of the columns
    # of CDi (called V1, V2, ..., "VtotV")
    colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
  }
  return(CDi)
}


#############################################################################
# Creation of matrices REFORD with ORDER

REFORDInit <- function(ORDER, nr, nc) {
  # The purpose of this function is to accelerate part 3.3 in which we
  # initially (i.e. with the first versions of seqimpute3.R) go "order"
  # times through the matrix ORDER
  #
  # Going one single time through the matrix ORDER, we create MaxGap
  # matrices REFORD (numbered from REFORD_1 to "REFORD_MaxGap") which
  # collect the coordinates of each corresponding values greater than 0
  # It will create MaxGap lookup matrices that will be used in point 3.3
  # to directly pinpoint the NA that have to be currently imputed
  # according to the value of the variable "order"
  # Updating MaxGap

  MaxGap <- max(ORDER[ORDER != 0]) - (min(ORDER[ORDER != 0]) - 1)

  # Initialization of the REFORD matrices
  REFORD_L <- lapply(1:MaxGap, matrix, data = NA, nrow = 0, ncol = 2)

  # Return the matrix of coordinates of value of ORDER bigger than 0
  non_zero <- which(ORDER > 0, arr.ind = TRUE)
  non_zero <- non_zero[(non_zero[, 1] <= nr) & 
      (non_zero[, 2] <= nc),, drop = FALSE]
  # Updating ORDER so that it becomes a matrix with positive
  # values going from 1 to MaxGap
  ORDER[non_zero] <- ORDER[non_zero] - (min(ORDER[ORDER != 0]) - 1)

  # Collecting the coordinates for each value of "order"
  non_zero <- non_zero[order(non_zero[, 1]), , drop = FALSE]
  ord_cord <- ORDER[non_zero]

  for (i in 1:MaxGap) {
    REFORD_L[[i]] <- non_zero[which(ord_cord == i), , drop = FALSE]
  }
  return(list(MaxGap=MaxGap, REFORD_L=REFORD_L, ORDER=ORDER))
}

############################################################################
# Creation of matrices REFORD with GapSize

REFORDInit_TI <- function(GapSize, nr, ORDER, GapSizelist) {
  # Initialization of the REFORD matrices
  REFORD_L <- lapply(1:GapSize, matrix, data = NA, nrow = 0, ncol = 2)

  # Return the matrix of coordinates of value of ORDER bigger than 0
  non_zero <- which(ORDER > 0, arr.ind = TRUE)
  non_zero <- non_zero[(non_zero[, 1] <= nr) & 
    (non_zero[, 2] %in% GapSizelist), , drop = FALSE]

  # Collecting the coordinates for each value of "order"
  ord_cord <- ORDER[non_zero, drop = FALSE]

  for (i in 1:GapSize) {
    REFORD_L[[i]] <- non_zero[which(ord_cord == i), ]
  }

  return(REFORD_L)
}



###############################################################
# Past VI computation

PastVICompute <- function(CD, CO, OD, ncot, frameSize, iter, nr, nc, ud, np, 
  COtsample, COt, pastDistrib, futureDistrib, k)
{
  CDp <- matrix(NA, nrow = nr * ud, ncol = np)

  if (ncot > 0) {
    # initialisation of matrix COtselected
    COtselected <- do.call(rbind, 
      replicate(ud, COtsample, simplify = FALSE))
  }

  # initialisation of matrix CDdb (for past
  # distribution analysis) (Distribution Before)
  if (pastDistrib) {
    CDdb <- matrix(NA, nrow = nr * ud, ncol = k)
    db <- matrix(NA, nrow = nr, ncol = k)
    # submatrix of CDdb:
    # CDdb is composed of
    # ud matrix db on top
    # of each other
  }

  # initialisation of matrix CDda (for future
  # distribution analysis) (Distribution After)
  if (futureDistrib) {
    CDda <- matrix(NA, nrow = nr * ud, ncol = k)
    # CDda has same dimensions as CDdb
    da <- matrix(NA, nrow = nr, ncol = k)
    # da has same dimensions as db
  }
  for (j in frameSize:nc) {
    # /!\ j is initialised at
    # the very end (utmost right) of the
    # frame
    t1 <- (nr * (iter - 1) + 1)
    # Determining the locations
    # of the time span (always nr) of
    # the piled up blocks of CD
    t2 <- nr * iter
    # VD
    CD[t1:t2, 1] <- OD[, j - frameSize + np + 1]


    # past VIs
    CDp[t1:t2, ] <- OD[, (j - frameSize + 1):(j - frameSize + np)]


    # /!\ current
    # pointer on column
    # is thus:
    # "j-frameSize

    # Eventually considering time-dependent
    # covariates
    if (ncot > 0) {
      COtselected <- COtselection(COtselected, COt, ncot, t1, t2, nr, 
        nc, j, shifted = -frameSize + np + 1)
    }
    # Past distribution (i.e. Before)
    if (pastDistrib) {
      ODt <- t(OD)
      ODt <- as.data.frame(ODt)
      tempOD <- lapply(ODt[(1:(j - frameSize + np)), ], factor, 
        levels = c(1:k, NA), exclude = NULL)
      # because:
      # j-frameSize+np+1 - 1 = j-frameSize+np

      db_list <- lapply(tempOD, summary)
      db_matrix <- do.call(rbind, db_list)
      CDdb[t1:t2, ] <- db_matrix[, 1:k] / length(1:(j - frameSize + np))
    }
    # Future distribution (i.e. After)
    if (futureDistrib) {
      if ((j - frameSize + np + 2) <= nc) {
        ODt <- t(OD)
        ODt <- as.data.frame(ODt)
        tempOD <- lapply(ODt[((j - frameSize + np + 2):nc), ], factor, 
          levels = c(1:k, NA), exclude = NULL)
        # because:
        # j-frameSize+np+1 + 1
        # = j-frameSize+np+2

        da_list <- lapply(tempOD, summary)
        da_matrix <- do.call(rbind, da_list)
        CDda[t1:t2, ] <- da_matrix[, 1:k] / length((j - frameSize + np + 2):nc)
      } else {
        # if index in OD exceeds OD number of
        # columns, the future distribution of
        # the k categorical variables is simply
        # null for everyone of them
        CDda[t1:t2, ] <- matrix(nrow = nr, ncol = k, 0)
      }
    }

    iter <- iter + 1
  }


  # past VIs
  CD <- cbind(CD, CDp)


  if (pastDistrib) {
    CD <- cbind(CD, CDdb)
  }

  if (futureDistrib) {
    CD <- cbind(CD, CDda)
  }

  # Conversion of CD into a data frame
  CD <- as.data.frame(CD)

  # Eventually concatenating CD with COs (the matrix
  # containing the covariates)
  if (all(is.na(CO)) == FALSE) {
    if (is.null(dim(CO))) {
      CO <- matrix(CO, nrow = nrow(OD), ncol = 1)
      COs <- do.call("rbind", rep(list(CO), ud))
      # Concatenating CD and COs into CD
      CD <- cbind(CD, COs)
    } else {
      # Checking if CO is NOT
      # completely empty
      # Creation of the stacked covariates
      # matrix for 3.1
      COs <- do.call("rbind", rep(list(CO), ud))
      # Concatenating CD and COs into CD
      CD <- cbind(CD, COs)
    }
  }
  # Else, in case CO is empty (i.e. we don't consider
  # any covariate) simply continue with the current CD

  # Eventually concatenating CD with COtselected (the
  # matrix containing the current time-dependent
  # covariates)
  # Checking if COt is NOT completely empty
  if (ncot > 0) {
    # Concatenating CD and COtselected into CD
    CD <- cbind(CD, as.data.frame(COtselected))
  }
  return(CD)
}


###############################################################################
# Future VI computation

FutureVICompute <- function(CD, CO, OD, ncot, frameSize, iter, nr, nc, ud, np, 
  COtsample, COt, pastDistrib, futureDistrib, k, nf)
{
  CDf <- matrix(NA, nrow = nr * ud, ncol = nf)

  if (ncot > 0) {
    # initialisation of matrix COtselected
    COtselected <- do.call(rbind, replicate(ud, COtsample, simplify = FALSE))
  }

  # initialisation of matrix CDdb (for past
  # distribution analysis) (Distribution Before)
  if (pastDistrib) {
    CDdb <- matrix(NA, nrow = nr * ud, ncol = k)
    db <- matrix(NA, nrow = nr, ncol = k)
    # submatrix of CDdb:
    # CDdb is composed of
    # ud matrix db on top
    # of each other
  }

  # initialisation of matrix CDda (for future
  # distribution analysis) (Distribution After)
  if (futureDistrib) {
    CDda <- matrix(NA, nrow = nr * ud, ncol = k)
    # CDda has same dimensions as CDdb
    da <- matrix(NA, nrow = nr, ncol = k)
    # da has same dimensions as db
  }
  for (j in frameSize:nc) {
    # /!\ j is initialised at
    # the very end (utmost right) of the
    # frame
    t1 <- (nr * (iter - 1) + 1)
    # Determining the locations
    # of the time span (always nr) of
    # the piled up blocks of CD
    t2 <- nr * iter
    # VD
    CD[t1:t2, 1] <- OD[, j - frameSize + np + 1]



    # future VIs
    CDf[t1:t2, ] <- OD[, (j - nf + 1):j]

    # /!\ current
    # pointer on column
    # is thus:
    # "j-frameSize

    # Eventually considering time-dependent
    # covariates
    if (ncot > 0) {
      COtselected <- COtselection(COtselected, COt, ncot, t1, t2, nr, nc, j, 
        shifted = -frameSize + np + 1)
    }
    # Past distribution (i.e. Before)
    if (pastDistrib) {
      ODt <- t(OD)
      ODt <- as.data.frame(ODt)
      tempOD <- lapply(ODt[(1:(j - frameSize + np)), ], factor,
        levels = c(1:k, NA), exclude = NULL)
      # because:
      # j-frameSize+np+1 - 1 = j-frameSize+np

      db_list <- lapply(tempOD, summary)
      db_matrix <- do.call(rbind, db_list)
      CDdb[t1:t2, ] <- db_matrix[, 1:k] / length(1:(j - frameSize + np))
    }
    # Future distribution (i.e. After)
    if (futureDistrib) {
      ODt <- t(OD)
      ODt <- as.data.frame(ODt)
      tempOD <- lapply(ODt[((j - frameSize + np + 2):nc), ], factor, 
        levels = c(1:k, NA), exclude = NULL)
      # because:
      # j-frameSize+np+1 + 1
      # = j-frameSize+np+2

      da_list <- lapply(tempOD, summary)
      da_matrix <- do.call(rbind, da_list)
      CDda[t1:t2, ] <- da_matrix[, 1:k] / length((j - frameSize + np + 2):nc)
    }
    iter <- iter + 1
  }


  # future VIs
  CD <- cbind(CD, CDf)

  if (pastDistrib) {
    CD <- cbind(CD, CDdb)
  }

  if (futureDistrib) {
    CD <- cbind(CD, CDda)
  }

  # Conversion of CD into a data frame
  CD <- as.data.frame(CD)

  # Eventually concatenating CD with COs (the matrix
  # containing the covariates)
  if (all(is.na(CO)) == FALSE) {
    if (is.null(dim(CO))) {
      CO <- matrix(CO, nrow = nrow(OD), ncol = 1)
      COs <- do.call("rbind", rep(list(CO), ud))
      # Concatenating CD and COs into CD
      CD <- cbind(CD, COs)
    } else {
      # Checking if CO is NOT
      # completely empty
      # Creation of the stacked covariates
      # matrix for 3.1
      COs <- do.call("rbind", rep(list(CO), ud))
      # Concatenating CD and COs into CD
      CD <- cbind(CD, COs)
    }
  }
  # Else, in case CO is empty (i.e. we don't consider
  # any covariate) simply continue with the current CD

  # Eventually concatenating CD with COtselected (the
  # matrix containing the current time-dependent
  # covariates)
  # Checking if COt is NOT completely empty
  if (ncot > 0) {
    # Concatenating CD and COtselected into CD
    CD <- cbind(CD, as.data.frame(COtselected))
  }
  return(CD)
}


###########################################################################
# Past or future VI computation

PastFutureVICompute <- function(CD, CO, OD, ncot, frameSize, iter, nr, nc, 
  ud, np, COtsample, COt, pastDistrib, futureDistrib, k, nf, shift)
{
  CDp <- matrix(NA, nrow = nr * ud, ncol = np)
  CDf <- matrix(NA, nrow = nr * ud, ncol = nf)

  if (ncot > 0) {
    # initialisation of matrix COtselected
    COtselected <- do.call(rbind, replicate(ud, COtsample, simplify = FALSE))
  }

  # initialisation of matrix CDdb (for past
  # distribution analysis) (Distribution Before)
  if (pastDistrib) {
    CDdb <- matrix(NA, nrow = nr * ud, ncol = k)
    db <- matrix(NA, nrow = nr, ncol = k)
    # submatrix of CDdb:
    # CDdb is composed of
    # ud matrix db on top
    # of each other
  }
  # initialisation of matrix CDda (for future
  # distribution analysis) (Distribution After)
  if (futureDistrib) {
    CDda <- matrix(NA, nrow = nr * ud, ncol = k)
    # CDda has same dimensions as CDdb
    da <- matrix(NA, nrow = nr, ncol = k)
    # da has same dimensions as db
  }
  for (j in frameSize:nc) {
    # /!\ j is initialised at
    # the very end (utmost right) of the
    # frame
    t1 <- (nr * (iter - 1) + 1)
    # Determining the locations
    # of the time span (always nr) of
    # the piled up blocks of CD
    t2 <- nr * iter
    # VD
    CD[t1:t2, 1] <- OD[, j - frameSize + np + 1 + shift]


    # past VIs
    CDp[t1:t2, ] <- OD[, (j - frameSize + 1):(j - frameSize + np)]

    # future VIs
    CDf[t1:t2, ] <- OD[, (j - nf + 1):j]

    # /!\ current
    # pointer on column
    # is thus:
    # "j-frameSize

    # Eventually considering time-dependent
    # covariates
    if (ncot > 0) {
      COtselected <- COtselection(COtselected, COt, ncot, t1, t2, nr, nc, j, 
        shifted = -frameSize + np + 1 + shift)
    }
    # Past distribution (i.e. Before)
    if (pastDistrib) {
      ODt <- t(OD)
      ODt <- as.data.frame(ODt)
      tempOD <- lapply(ODt[(1:(j - frameSize + np + shift)), ], factor, 
        levels = c(1:k, NA), exclude = NULL)
      # because:
      # j-frameSize+np+1 - 1 = j-frameSize+np

      db_list <- lapply(tempOD, summary)
      db_matrix <- do.call(rbind, db_list)
      CDdb[t1:t2, ] <- db_matrix[, 1:k] / length(1:(j-frameSize + np + shift))
    }
    # Future distribution (i.e. After)
    if (futureDistrib) {
      ODt <- t(OD)
      ODt <- as.data.frame(ODt)
      tempOD <- lapply(ODt[((j - frameSize + np + shift + 2):nc), ], factor, 
        levels = c(1:k, NA), exclude = NULL)
      # because:
      # j-frameSize+np+1 + 1
      # = j-frameSize+np+2
      da_list <- lapply(tempOD, summary)
      da_matrix <- do.call(rbind, da_list)
      CDda[t1:t2,] <- da_matrix[,1:k]/length((j-frameSize+np+shift+2):nc)
    }

    iter <- iter + 1
  }

  # past and future VIs
  CD <- cbind(CD, CDp, CDf)



  if (pastDistrib) {
    CD <- cbind(CD, CDdb)
  }

  if (futureDistrib) {
    CD <- cbind(CD, CDda)
  }

  # Conversion of CD into a data frame
  CD <- as.data.frame(CD)
  # Eventually concatenating CD with COs (the matrix
  # containing the covariates)
  if (all(is.na(CO)) == FALSE) {
    if (is.null(dim(CO))) {
      CO <- matrix(CO, nrow = nrow(OD), ncol = 1)
      COs <- do.call("rbind", rep(list(CO), ud))
      # Concatenating CD and COs into CD
      CD <- cbind(CD, COs)
    } else {
      # Checking if CO is NOT
      # completely empty
      # Creation of the stacked covariates
      # matrix for 3.1
      COs <- do.call("rbind", rep(list(CO), ud))
      # Concatenating CD and COs into CD
      CD <- cbind(CD, COs)
    }
  }
  # Else, in case CO is empty (i.e. we don't consider
  # any covariate) simply continue with the current CD
  # Eventually concatenating CD with COtselected (the
  # matrix containing the current time-dependent
  # covariates)
  # Checking if COt is NOT completely empty
  if (ncot > 0) {
    # Concatenating CD and COtselected into CD
    CD <- cbind(CD, as.data.frame(COtselected))
  }
  return(CD)
}



###########################################################################
# Imputation where only PAST VIs  exist

ODiImputePAST <- function(CO, ODi, CD, COt, REFORD, nr_REFORD, pastDistrib, 
  futureDistrib, k, np, nf, nc, ncot, reglog, LOOKUP, regr, noise)
{
  for (u in 1:nr_REFORD) {
    i <- REFORD[u, 1]
    # taking out the first coordinate
    # (row number in ORDER) from REFORD
    j <- REFORD[u, 2]
    # taking out the second coordinate
    # (column number in ORDER) from REFORD


    CDi <- matrix(NA, nrow = k, ncol = 1)


    # Matrix for past values
    vect <- LOOKUP[i, (j - np):(j - 1)]
    # /!\ current pointer
    # on olumn is thus: "j"
    CDpi <- matrix(vect, nrow = k, ncol = length(vect), byrow = TRUE)

    # Matrix for past distribution
    if (pastDistrib) {
      dbi <- summary(factor(LOOKUP[i, 1:(j - 1)], levels = c(1:k), 
        exclude = NULL)) / length(1:(j - 1))
      CDdbi <- matrix(dbi[1:k], nrow = k, ncol = k, byrow = TRUE)
    }

    # Matrix for future distribution
    if (futureDistrib) {
      dai <- summary(factor(LOOKUP[i, (j + 1):nc], levels = c(1:k), 
        exclude = NULL)) / length((j + 1):nc)
      CDdai <- matrix(dai[1:k], nrow = k, ncol = k, byrow = TRUE)
    }

    # Concatenating CDi
    CDi <- cbind(CDi, CDpi)

    if (pastDistrib) {
      CDi <- cbind(CDi, CDdbi)
    }

    if (futureDistrib) {
      CDi <- cbind(CDi, CDdai)
    }

    # Conversion of CDi into a data frame
    CDi <- as.data.frame(CDi)
    # Type transformation of the columns of CDi
    # The first values of CDi must be of type factor
    # (categorical values)



    # if (regr == "lm" | regr == "lrm") {
    #   for (v in 1:(1 + np + nf)) {
    #     CDi[, v] <- factor(CDi[, v], levels = levels(CD[, v]), exclude = NULL)
    #   }
    # } else if (regr == "rf") {
    #   for (v in 2:(1 + np + nf)) {
    #     CDi[, v] <- factor(CDi[, v], levels = c(1:(k + 1)))
    #     CDi[, v][is.na(CDi[, v])] <- k + 1
    #   }
    #   CDi[, 1] <- factor(CDi[, 1], levels = levels(CD[, 1]))
    # } else {
    #   # multinom
    #   CDi[, 1] <- factor(CDi[, 1], levels = c(1:k))
    #   for (v in 2:(1 + np + nf)) {
    #     CDi[, v] <- factor(CDi[, v], levels = levels(CD[, v]), exclude = NULL)
    #   }
    # }
    
    if (regr == "rf") {
      for (v in 2:(1 + np + nf)) {
        CDi[, v] <- factor(CDi[, v], levels = c(1:(k + 1)))
        CDi[, v][is.na(CDi[, v])] <- k + 1
      }
      CDi[, 1] <- factor(CDi[, 1], levels = levels(CD[, 1]))
    } else {
      # multinom
      CDi[, 1] <- factor(CDi[, 1], levels = c(1:k))
      for (v in 2:(1 + np + nf)) {
        CDi[, v] <- factor(CDi[, v], levels = levels(CD[, v]), exclude = NULL)
      }
    }
    
    
    # The last values of CDi must be of type numeric
    # (distributions)
    if (pastDistrib | futureDistrib) {
      CDi[, (1 + np + nf + 1):ncol(CDi)] <- lapply(
        CDi[, (1 + np + nf + 1):ncol(CDi)], as.numeric)
    }

    # Eventually concatenating CDi with COi
    # (the matrix containing the covariates)
    if (all(is.na(CO)) == FALSE) {
      # Checking if CO is NOT
      # completely empty
      # Creation of the matrix COi used in 3.3
      if (is.null(dim(CO))) {
        COi <- do.call(rbind, 
          replicate(k, as.data.frame(CO[i]), simplify = FALSE))
      } else {
        COi <- do.call(rbind, 
          replicate(k, as.data.frame(CO[i, ]), simplify = FALSE))
      }
      # Concatenating CDi and COi into CDi
      CDi <- cbind(CDi, COi)
      # Transformation of the names of the columns
      # of CDi (called V1, V2, ..., "VtotV")
      colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
    }
    # Else, in case CO is empty (i.e. we don't
    # consider any covariate)
    # simply continue with the current CDi

    # Eventually concatenating CDi with
    # COtselected_i (the matrix containing the
    # current time-dependent covariates)
    # Checking if COt is NOT completely empty
    if (ncot > 0) {
      COtselected_i <- as.data.frame(matrix(nrow = 1, ncol = 0))
      for (d in 1:(ncot / nc)) {
        COtselected_i <- cbind(COtselected_i, 
          COt[i, (j) + (d - 1) * nc])
      }
      COtselected_i <- do.call(rbind, replicate(k, 
          as.data.frame(COtselected_i[1, ]), simplify = FALSE))
      # Concatenating CDi and COtselected_i
      # into CDi
      CDi <- cbind(CDi, COtselected_i)
      # Transformation of the names of the columns
      # of CDi (called V1, V2, ..., "VtotV")
      colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
    }


    # Check for missing-values among predictors
    if (max(is.na(CDi[1, 2:ncol(CDi)])) == 0) {
      ODi <- RegrImpute(ODi, CDi, regr, reglog, noise, i, j, k)
    }
  }
  return(ODi)
}



##################################################################
# Imputation where past and future VIs exist

ODiImputePF <- function(CO, ODi, CD, COt, REFORD, nr_REFORD, pastDistrib, 
  futureDistrib, k, np, nf, nc, ncot, reglog, LOOKUP, regr, noise, 
  shift, MaxGap, order)
{
  for (u in 1:nr_REFORD) {
    i <- REFORD[u, 1]
    # taking out the first coordinate
    # (row number in ORDER) from REFORD
    j <- REFORD[u, 2]

    # taking out the second coordinate
    # (column number in ORDER) from REFORD
    CDi <- matrix(NA, nrow = k, ncol = 1)

    # Matrix for past values
    shift <- as.numeric(shift)

    vect <- LOOKUP[i, (j - shift - np):(j - shift - 1)]

    CDpi <- matrix(vect, nrow = k, ncol = length(vect), byrow = TRUE)

    # Matrix for future values
    vect <- LOOKUP[i, 
      (j - shift + MaxGap - order + 1):(j - shift + MaxGap - order + nf)]
    CDfi <- matrix(vect, nrow = k, ncol = length(vect), byrow = TRUE)

    # Matrix for past distribution
    if (pastDistrib) {
      dbi <- summary(factor(LOOKUP[i, 1:(j - 1)], 
        levels = c(1:k), exclude = NULL)) / length(1:(j - 1))
      CDdbi <- matrix(dbi[1:k], nrow = k, ncol = k, byrow = TRUE)
    }
    # Matrix for future distribution
    if (futureDistrib) {
      dai <- summary(factor(LOOKUP[i, (j + 1):nc], 
        levels = c(1:k), exclude = NULL)) / length((j + 1):nc)
      CDdai <- matrix(dai[1:k], nrow = k, ncol = k, byrow = TRUE)
    }
    # Concatenating CDi
    CDi <- cbind(CDi, CDpi, CDfi)

    if (pastDistrib) {
      CDi <- cbind(CDi, CDdbi)
    }

    if (futureDistrib) {
      CDi <- cbind(CDi, CDdai)
    }

    CDi <- as.data.frame(CDi)
   
    
    if (regr == "rf") {
      for (v in 2:(1 + np + nf)) {
        CDi[, v] <- factor(CDi[, v], levels = c(1:(k + 1)))
        CDi[, v][is.na(CDi[, v])] <- k + 1
      }
      CDi[, 1] <- factor(CDi[, 1], levels = levels(CD[, 1]))
    } else {
      # multinom
      CDi[, 1] <- factor(CDi[, 1], levels = c(1:k))
      for (v in 2:(1 + np + nf)) {
        CDi[, v] <- factor(CDi[, v], levels = levels(CD[, v]), exclude = NULL)
      }
    }
    
    
    if (pastDistrib | futureDistrib) {
      CDi[, (1 + np + nf + 1):ncol(CDi)] <- lapply(
        CDi[, (1 + np + nf + 1):ncol(CDi)], as.numeric
        )
    }


    if (all(is.na(CO)) == FALSE) {
      # Checking if CO is NOT
      # completely empty
      # Creation of the matrix COi used in 3.3
      if (is.null(dim(CO))) {
        COi <- do.call(rbind, replicate(k, as.data.frame(CO[i]), 
          simplify = FALSE))
      } else {
        COi <- do.call(rbind, replicate(k, as.data.frame(CO[i, ]), 
          simplify = FALSE))
      }
      # Concatenating CDi and COi into CDi
      CDi <- cbind(CDi, COi)
      # Transformation of the names of the columns
      # of CDi (called V1, V2, ..., "VtotV")
      colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
    }
    # Else, in case CO is empty (i.e. we don't
    # consider any covariate)
    # simply continue with the current CDi

    # Eventually concatenating CDi with
    # COtselected_i (the matrix containing the
    # current time-dependent covariates)
    # Checking if COt is NOT completely empty
    if (ncot > 0) {
      COtselected_i <- as.data.frame(matrix(nrow = 1, ncol = 0))
      for (d in 1:(ncot / nc)) {
        COtselected_i <- cbind(COtselected_i, COt[i, (j) + (d - 1) * nc])
      }
      COtselected_i <- do.call(rbind, replicate(k, 
        as.data.frame(COtselected_i[1, ]), simplify = FALSE))
      # Concatenating CDi and COtselected_i into
      # CDi
      CDi <- cbind(CDi, COtselected_i)
      # Transformation of the names of the columns
      # of CDi (called V1, V2, ..., "VtotV")
      colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
    }


    # Check for missing-values among predictors
    # (i.e. we won't impute any value on the current
    # MD if there is any NA among the VIs)
    if (max(is.na(CDi[1, 2:ncol(CDi)])) == 0) {
      # checking that
      # there is no NA
      # among the current
      # VI (otherwise no
      # data will be
      # imputed for the
      # current NA)
      ODi <- RegrImpute(ODi, CDi, regr, reglog, noise, i, j, k)
    }
  }
  return(ODi)
}


RegrImpute <- function(ODi, CDi, regr, reglog, noise, i, j, k) {
  if (regr == "multinom") {
    if (length(reglog$lev) > 2) {
      pred <- predict(reglog, CDi, type = "probs")[1, ]
      pred <- cumsum(pred)

      names_saved <- names(pred)


      alea <- runif(1)

      sel <- as.numeric(names_saved[which(pred >= alea)])
    } else {
      pred <- predict(reglog, CDi, type = "probs")
      alea <- runif(1)
      if (alea > pred[1]) {
        sel <- as.numeric(reglog$lev[1])
      } else {
        sel <- as.numeric(reglog$lev[2])
      }
    }
  } else if (regr == "rf") { # regr=="rf" randomForest
    # pred <- predict(reglog,newdata=CDi,type="prob")[1,]
    pred <- predict(reglog, data = CDi, predict.all = TRUE)$predictions[1, ]
    pred <- factor(pred, levels = c(1:k))
    tab <- table(pred)
    tab <- tab / sum(tab)
    
    pred <- cumsum(tab)

    
    alea <- runif(1)
    
    sel <- levels(CDi[, 1])[which(pred >= alea)[1]]
    
  } else {


  }


  ODi[i, j] <- sel[1]
  return(ODi)
}


################################################################################
# Compute model with the chosen regression model

ComputeModel <- function(CD, regr, np, nf, k, ...) {
  npfi <- np + nf
  
  CD <- as.data.frame(CD)


  colnames(CD) <- paste("V", 1:ncol(CD), sep = "")

  if (regr == "rf") { 

    CD[, (1:(1 + npfi))] <- lapply(CD[, (1:(1 + npfi))], 
      factor, levels = c(1:k))
    for (v in 2:(1 + npfi)) {
      CD[, v] <- factor(CD[, v], levels = c(1:(k + 1)))
      CD[, v][is.na(CD[, v])] <- k + 1
    }
    CD[, 1] <- factor(CD[, 1], levels = c(1:k))
    CD[, 1] <- droplevels(CD[, 1])
    CD <- CD[!is.na(CD[, 1]), ]


    factors_character <- paste("V", 2:ncol(CD), sep = "")
    
    factors <- as.vector(factors_character)
    
    fmla <- as.formula(paste("V1~", paste(factors, collapse = "+")))

    if ("num.trees" %in% names(list(...))) {
      reglog <- ranger(fmla, data = CD, ...)
    } else {
      reglog <- ranger(fmla, data = CD, num.trees = 10, ...)
    }
  } else if (regr == "multinom") {
    CD[, 1] <- factor(CD[, 1], levels = c(1:k))
    CD[, 1] <- droplevels(CD[, 1])



    if (npfi > 1) {
      CD[, (2:(1 + npfi))] <- lapply(CD[, (2:(1 + npfi))], factor, 
        levels = c(1:k, NA), exclude = NULL)
    } else {
      CD[, 2] <- factor(CD[, 2], levels = c(1:k, NA), exclude = NULL)
    }


    factors_character <- paste("V", 2:ncol(CD), sep = "")
   
    factors <- as.vector(factors_character)
    
    fmla <- as.formula(paste("V1~", paste(factors, collapse = "+")))

    reglog <- nnet::multinom(CD, maxit = 100, trace = FALSE, 
      MaxNwts = 1500, ...)
  } else {
  }

  return(list(reglog, CD))
}
