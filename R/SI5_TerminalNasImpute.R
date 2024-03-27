# 5. Imputation of terminal NAs
#
#
################################################################################
# Impute terminal NAs

ImputingTerminalNAs <- function(OD, covariates, time.covariates, ODi, 
  COtsample, MaxTermGapSize, TermGapSize, pastDistrib, regr, npt, nco, ncot, 
  totVt, nr, nc, ud, available, k, noise, ...)
{
  # 5.1.-2. Creation of ORDERT -------------------------------------------------
  REFORDT_L <- REFORDTCreation(nr, nc, TermGapSize, MaxTermGapSize)

  # Case when there is not enough observations before the longest terminal gap
  # -> npt is therefore reduced
  if (npt > nc - MaxTermGapSize) {
    npt <- nc - MaxTermGapSize
    totVt <- 1 + npt + nco + (ncot / nc)
    if (pastDistrib) {
      totVt <- totVt + k
    }
  }
  # 5.3. Imputation using a specific model -------------------------------------

  # 5.3.1 Building of the data matrix CD for the computation of the model ------
  CD <- TerminalCDMatCreate(covariates, OD, time.covariates, COtsample, 
    pastDistrib, npt, nr, nc, ncot, k)
  # 5.3.2 Computation of the model (Dealing with the LOCATIONS of imputation) --
  log_CD <- list()
  log_CD[c("reglog", "CD")] <- ComputeModel(CD, regr, totVt, npt, 0, k, ...)
  # 5.3.3 Imputation using the just created model 
  ODi <- TerminalCreatedModelImputation(covariates, OD, log_CD$CD, ODi, 
    time.covariates, nc, ncot, totVt, REFORDT_L, pastDistrib, MaxTermGapSize, 
    available, regr, log_CD$reglog, k, npt, noise)

  return(ODi)
}



################################################################################
# Create EFORDT matrix

REFORDTCreation <- function(nr, nc, TermGapSize, MaxTermGapSize) {
  # Creation of matrix ORDERT
  ORDERT <- matrix(0, nr, nc)
  for (i in 1:nr) {
    if (TermGapSize[i] != 0) {
      ORDERT[i,(nc-TermGapSize[i]+1):nc] <- c((MaxTermGapSize+1-TermGapSize[i]):
        MaxTermGapSize)
    } else {
      next
    }
  }

  # In a similar manner to part 2.4., we go here one single
  # time through a reduced version of ORDERT and we create
  # MaxTermGapSize REFORDT matrices collecting the coordinates
  # of each corresponding values in ORDERT which are greater
  # than 0

  REFORDT_L <- REFORDInit_TI(MaxTermGapSize,nr,ORDERT, (nc-MaxTermGapSize+1):nc)
  return(REFORDT_L)
}



################################################################################
# Building of the data matrix CD for the computation of the model

TerminalCDMatCreate <- function(CO, OD, COt, COtsample, pastDistrib, npt, nr, 
  nc, ncot, k) {
  # For Terminal Gaps
  # we will impute single NA after single NA going from the
  # center of OD towards its very right border
  # We here only take observations from the PAST
  # into account --> npt
  # But we will create one single regression
  # model used for the imputation of all the NAs belonging to
  # a "terminal gap"
  frameSize <- npt + 1
  ud <- nc - frameSize + 1

  CD <- matrix(NA, nrow = nr * ud, ncol = 1)

  # initialisation of matrix CDp
  CDp <- matrix(NA, nrow = nr * ud, ncol = npt)
  if (ncot > 0) {
    # initialisation of matrix COtselected
    COtselected <- do.call(rbind, replicate(ud, COtsample, simplify = FALSE))
  }

  # initialisation of matrix CDdb
  # (for past distribution analysis)
  # (Distribution Before)
  if (pastDistrib) {
    CDdb <- matrix(NA, nrow = nr * ud, ncol = k)
    db <- matrix(NA, nrow = nr, ncol = k)
    # submatrix of CDdb: CDdb is
    # composed of ud matrix db on
    # top of each other
  }

  iter <- 1

  for (j in frameSize:nc) {
    t1 <- (nr * (iter - 1) + 1)
    t2 <- nr * iter

    # VD
    CD[t1:t2, 1] <- OD[, j]
    # /!\ current pointer on column is thus: "j"

    # Past VIs
    CDp[t1:t2, ] <- OD[, (j - npt):(j - 1)]

    # Eventually considering time-dependent covariates
    if (ncot > 0) {
      COtselected <- COtselection(COtselected, COt, ncot, t1, t2, nr, nc, j, 
        shifted = 0)
    }

    # Past distribution (i.e. Before)
    if (pastDistrib) {
      ODt <- t(OD)
      ODt <- as.data.frame(ODt)
      tempOD <- lapply(ODt[(1:(j - 1)), ], factor, levels = c(1:k, NA),
        exclude = NULL)
      db_list <- lapply(tempOD, summary)
      db_matrix <- do.call(rbind, db_list)
      CDdb[t1:t2, ] <- db_matrix[, 1:k] / length(1:(j - 1))
    }

    iter <- iter + 1
  }

  # Concatening CD
  CD <- cbind(CD, CDp)

  if (pastDistrib) {
    CD <- cbind(CD, CDdb)
  }

  # Conversion of CD into a data frame
  CD <- as.data.frame(CD)

  # Eventually concatenating CD with COs
  # (the matrix containing the covariates)
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
  # covariate) simply continue with the current CD

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
# Imputation of the terminal created model

TerminalCreatedModelImputation <- function(CO, OD, CD, ODi, COt, nc, ncot, 
  totVt, REFORDT_L, pastDistrib, MaxTermGapSize, available, regr, reglog, k, 
  npt, noise)
{
  # Conversion of ODi from data.frame to matrix
  ODi <- as.matrix(ODi)

  # Only PAST VIs are useful
  for (order in 1:MaxTermGapSize) {
    # Analysing the value of parameter available
    if (available == TRUE) {
      # we take the previously
      # imputed data into account
      LOOKUP <- ODi
    } else { # that is available == FALSE and thus we
      # don't take the previously imputed
      # data into account
      LOOKUP <- OD
    }

    # Assigning the current "REFORDT_order" matrix to the
    # variable matrix REFORDT
    # (according to the current value of "order")
    REFORDT <- as.matrix(REFORDT_L[[order]])
    if (ncol(REFORDT) == 1) {
      REFORDT <- t(REFORDT)
    }
    nr_REFORDT <- nrow(REFORDT)

    for (u in 1:nr_REFORDT) {
      i <- REFORDT[u, 1] # taking out the first coordinate
      # (row number in ORDER) from REFORDT
      j <- REFORDT[u, 2] # taking out the second coordinate
      # (column number in ORDER) from REFORDT

      CDi <- matrix(NA, nrow = k, ncol = 1)

      # Matrix for past values
      vect <- LOOKUP[i, (j - npt):(j - 1)]
      CDpi <- matrix(vect, nrow = k, ncol = length(vect), byrow = TRUE)

      # Matrix for past distribution
      if (pastDistrib) {
        dbi <- summary(factor(LOOKUP[i, 1:(j - 1)], levels = c(1:k), 
          exclude = NULL)) / length(1:(j - 1))
        CDdbi <- matrix(dbi[1:k], nrow = k, ncol = k, byrow = TRUE)
      }

      # Concatenating CDi
      CDi <- cbind(CDi, CDpi)

      if (pastDistrib) {
        CDi <- cbind(CDi, CDdbi)
      }

      # Conversion of CDi into a data frame
      CDi <- as.data.frame(CDi)

      # Type transformation of the columns of CDi
      # The first values of CDi must be of type factor
      # (categorical values)


      if (regr != "rf") {
        for (v in 1:(1 + npt)) {
          CDi[, v] <- factor(CDi[, v], levels = levels(CD[, v]), exclude = NULL)
        }
      } else {
        for (v in 2:(1 + npt)) {
          CDi[, v] <- factor(CDi[, v], levels = c(1:(k + 1)))
          CDi[, v][is.na(CDi[, v])] <- k + 1
        }
        CDi[, 1] <- factor(CDi[, 1], levels = levels(CD[, 1]))
      }
      # The
      # The last values of CDi must be of type numeric
      # (distributions)
      if (pastDistrib) {
        CDi[, (1+npt+1):ncol(CDi)] <- lapply(CDi[, (1 + npt + 1):ncol(CDi)], 
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
      # any covariate)
      # simply continue with the current CDi

      # Eventually concatenating CDi with COtselected_i
      # (the matrix containing the current time-dependent
      # covariates)
      # Checking if COt is NOT completely empty
      CDi <- COtselectionSpe(CDi, COt, ncot, nc, i, j, k)

      # Check for missing-values among predictors
      if (max(is.na(CDi[1, 2:totVt])) == 0) {
        ODi <- RegrImpute(ODi, CDi, regr, reglog, noise, i, j, k)
      }
    }
  }
  return(ODi)
}
