# 4. Imputation of initial NAs
#
#
################################################################################
# Impute initial NAs

ImputingInitialNAs <- function(OD, covariates, time.covariates, ODi, totVi, 
  COtsample, futureDistrib, InitGapSize, MaxInitGapSize, nr, nc, ud, nco, ncot, 
  nfi, regr, k, available, noise, ...) 
{
  # 4.1.-2. Creation of ORDERI -------------------------------------------------
  REFORDI_L <- REFORDICreation(nr, nc, InitGapSize, MaxInitGapSize)


  # Case when there is not enough observations after the longest initial gap
  # -> nfi is therefore reduced
  if (nfi > nc - MaxInitGapSize) {
    nfi <- nc - MaxInitGapSize
    totVi <- 1 + nfi + nco + (ncot / nc)
    if (futureDistrib) {
      totVi <- totVi + k
    }
  }
  # 4.3. Imputation using a specific model -------------------------------------

  # 4.3.1 Building of the data matrix CD for the computation of the model ------
  CD <- CDMatCreate(covariates, COtsample, OD, time.covariates, nfi, nr, nc, 
    ncot, futureDistrib, k)

  # 4.3.2 Computation of the model (Dealing with the LOCATIONS of imputation) --
  log_CD <- list()
  log_CD[c("reglog", "CD")] <- ComputeModel(CD, regr, totVi, 0, nfi, k, ...)

  # 4.3.3 Imputation using the just created model 
  ODi <- Init_NA_CreatedModelImputation(OD, ODi, covariates, log_CD$CD, 
    time.covariates, MaxInitGapSize, REFORDI_L, futureDistrib, totVi, nc, k, 
    nfi, ncot, regr, log_CD$reglog, noise, available)
  return(ODi)
}


################################################################################
# Create ORDERI matrix

REFORDICreation <- function(nr, nc, InitGapSize, MaxInitGapSize) {
  # Creation of matrix ORDERI
  ORDERI <- matrix(0, nr, nc)
  for (i in 1:nr) {
    if (InitGapSize[i] != 0) {
      ORDERI[i, 1:InitGapSize[i]] <- c(MaxInitGapSize:
        (MaxInitGapSize + 1 - InitGapSize[i]))
    } else {
      next
    }
  }

  # In a similar manner to part 2.4., we go here one single
  # time through a reduced version of ORDERI and we create
  # MaxInitGapSize REFORDI matrices collecting the coordinates
  # of each corresponding values in ORDERI which are greater
  # than 0
  REFORDI_L <- REFORDInit_TI(MaxInitGapSize, nr, ORDERI, MaxInitGapSize:1)

  return(REFORDI_L)
}



################################################################################
# Building of the data matrix CD for the computation of the model

CDMatCreate <- function(CO, COtsample, OD, COt, nfi, nr, nc, 
  ncot, futureDistrib, k) {
  # For Initial Gaps
  # we will impute single NA after single NA going from the
  # center of OD towards its very left border
  # We here only take observations from the FUTURE into
  # account --> nfi
  # But we will create one single regression
  # model used for the imputation of all
  # the NAs belonging to an "initial gap"
  frameSize <- 1 + nfi
  ud <- nc - frameSize + 1

  CD <- matrix(NA, nrow = nr * ud, ncol = 1)

  # initialisation of matrix CDf
  CDf <- matrix(NA, nrow = nr * ud, ncol = nfi)
  if (ncot > 0) {
    # initialisation of matrix COtselected
    COtselected <- do.call(rbind, replicate(ud, COtsample, simplify = FALSE))
  }

  # initialisation of matrix CDda (for future distribution
  # analysis)
  # (Distribution After)
  if (futureDistrib) {
    CDda <- matrix(NA, nrow = nr * ud, ncol = k)
    da <- matrix(NA, nrow = nr, ncol = k)
    # submatrix of CDda: CDda is
    # composed of ud matrix da on
    # top of each other
  }

  iter <- 1

  for (j in frameSize:nc) {
    t1 <- (nr * (iter - 1) + 1)
    t2 <- nr * iter

    # VD
    CD[t1:t2, 1] <- OD[, j - frameSize + 1]
    # /!\ current pointer on
    # column is thus:
    # "j-frameSize+1"

    # Future VIs
    CDf[t1:t2, ] <- OD[, (j - nfi + 1):j]

    # Eventually considering time-dependent covariates
    if (ncot > 0) {
      COtselected <- COtselection(COtselected, COt, ncot, t1, t2, nr, nc, j, 
        shifted = -frameSize + 1)
    }

    # Future distribution (i.e. After)
    if (futureDistrib) {
      ODt <- t(OD)
      ODt <- as.data.frame(ODt)
      tempOD <- lapply(ODt[((j - frameSize + 2):nc), ], factor, 
        levels = c(1:k, NA), exclude = NULL)
      da_list <- lapply(tempOD, summary)
      da_matrix <- do.call(rbind, da_list)
      CDda[t1:t2, ] <- da_matrix[, 1:k] / length((j - frameSize + 2):nc)
    }

    iter <- iter + 1
  }

  # Concatenating CD
  CD <- cbind(CD, CDf)

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
  # Else, in case CO is empty (i.e. we don't consider any
  # covariate)
  # simply continue with the current CD

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



################################################################################
# Imputation using the just created model 

Init_NA_CreatedModelImputation <- function(OD, ODi, CO, CD, COt, 
  MaxInitGapSize, REFORDI_L, futureDistrib, totVi, nc, k, nfi, ncot, regr, 
  reglog, noise, available)
{
  # Conversion of ODi from data.frame to matrix
  ODi <- as.matrix(ODi)

  # Only FUTURE VIs are useful
  for (order in 1:MaxInitGapSize) {
    # Analysing the value of parameter available
    if (available == TRUE) {
      # we take the previously imputed data
      # into account
      LOOKUP <- ODi
    } else {
      # that is available == FALSE and thus we don't take
      # the previously imputed data into account
      LOOKUP <- OD
    }

    # Assigning the current "REFORDI_order" matrix to the
    # variable matrix REFORDI
    # (according to the current value of "order")
    # tempObject = get(paste0("REFORDI_",order))
    REFORDI <- as.matrix(REFORDI_L[[order]])
    if (ncol(REFORDI) == 1) {
      REFORDI <- t(REFORDI)
    }
    nr_REFORDI <- nrow(REFORDI)

    for (u in 1:nr_REFORDI) {
      i <- REFORDI[u, 1]
      # taking out the first coordinate (row
      # number in ORDER) from REFORDI
      j <- REFORDI[u, 2]
      # taking out the second coordinate
      # (column number in ORDER) from REFORDI

      CDi <- matrix(NA, nrow = k, ncol = 1)

      # Matrix for future values
      vect <- LOOKUP[i, (j + 1):(j + nfi)]
      CDfi <- matrix(vect, nrow = k, ncol = length(vect), byrow = TRUE)

      # Matrix for future distribution
      if (futureDistrib) {
        dai <- summary(factor(LOOKUP[i, (j + 1):nc], levels = c(1:k), 
          exclude = NULL)) / length((j + 1):nc)
        CDdai <- matrix(dai[1:k], nrow = k, ncol = k, byrow = TRUE)
      }

      # Concatenating CDi
      CDi <- cbind(CDi, CDfi)

      if (futureDistrib) {
        CDi <- cbind(CDi, CDdai)
      }

      # Conversion of CDi into a data frame
      CDi <- as.data.frame(CDi)

      # Type transformation of the columns of CDi
      # The first values of CDi must be of type factor
      # (categorical values)

      if (regr != "rf") {
        for (v in 1:(1 + nfi)) {
          CDi[, v] <- factor(CDi[, v], levels = levels(CD[, v]), exclude = NULL)
        }
      } else {
        for (v in 2:(1 + nfi)) {
          CDi[, v] <- factor(CDi[, v], levels = c(1:(k + 1)))
          CDi[, v][is.na(CDi[, v])] <- k + 1
        }
        CDi[, 1] <- factor(CDi[, 1], levels = levels(CD[, 1]))
      }
      # The last values of CDi must be of type numeric
      # (distributions)
      if (futureDistrib) {
        CDi[, (1+nfi+1):ncol(CDi)] <- lapply(CDi[, (1 + nfi + 1):ncol(CDi)], 
          as.numeric)
      }
      # Eventually concatenating CDi with COi (the matrix
      # containing the covariates)
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
      # Else, in case CO is empty (i.e. we don't consider
      # any covariate) simply continue with the current
      # CDi

      # Eventually concatenating CDi with COtselected_i
      # (the matrix containing the current time-dependent
      # covariates)
      # Checking if COt is NOT completely empty
      CDi <- COtselectionSpe(CDi, COt, ncot, nc, i, j, k)

      # Check for missing-values among predictors
      if (max(is.na(CDi[1, 2:totVi])) == 0) {
        ODi <- RegrImpute(ODi, CDi, regr, reglog, noise, i, j, k)
      }
    }
  }
  return(ODi)
}
