seqimpute_timing <- function(data, regr = "multinom", np = 1, nf = 0, nfi = 1, 
  npt = 1, available = TRUE, covariates = matrix(NA, nrow = 1, ncol = 1),
  time.covariates = matrix(NA, nrow = 1, ncol = 1), pastDistrib = FALSE,
  futureDistrib = FALSE, m = 1, noise = 0, SetRNGSeed = FALSE, ParExec = TRUE, 
  ncores = NULL, frame.radius = 0, verbose = TRUE, ...)
{

  if (sum(is.na(data)) == 0) {
    if (verbose == TRUE) {
      message("This dataset has no missing values!")
    }
    return(data)
  }

  rownamesDataset <- rownames(data)
  nrowsDataset <- nrow(data)

  # 0. Initial tests and manipulations on parameters --------------------------
  dataOD <- preliminaryChecks(data, covariates, time.covariates, np, nf, nfi, 
    npt, pastDistrib, futureDistrib)
  dataOD[c("pastDistrib", "futureDistrib", "totV", "totVi", "totVt", 
    "noise")] <- InitCorectControl(regr, dataOD$ODClass, dataOD$OD, dataOD$nr, 
        dataOD$nc, dataOD$k, np, nf, dataOD$nco, dataOD$ncot, nfi, npt, 
        pastDistrib, futureDistrib, dataOD$totV, dataOD$totVi, 
        dataOD$totVt, noise)


  # 1. Analysis of OD and creation of matrices ORDER, ORDER2 and ORDER3 
  dataOD[c("MaxInitGapSize", "InitGapSize", "MaxTermGapSize", "TermGapSize",
    "MaxGap", "ORDER", "ORDER2", "ORDER3")] <- OrderCreation(dataOD$OD, 
        dataOD$nr, dataOD$nc)



  # 2. Computation of the order of imputation of each MD 
  if (max(dataOD$ORDER) != 0) {
    dataOD[c("ORDERSLGLeft", "ORDERSLGRight", "ORDERSLGBoth", "LongGap", 
      "MaxGap", "REFORD_L", "ORDER")] <- ImputeOrderComputation(dataOD$ORDER, 
        dataOD$ORDER3, dataOD$MaxGap, np, nf, dataOD$nr, dataOD$nc)
  } else {
    dataOD$ORDERSLGLeft <- matrix(nrow = dataOD$nr, ncol = dataOD$nc, 0)
    dataOD$ORDERSLGRight <- matrix(nrow = dataOD$nr, ncol = dataOD$nc, 0)
    dataOD$ORDERSLGBoth <- matrix(nrow = dataOD$nr, ncol = dataOD$nc, 0)
    dataOD$LongGap <- FALSE
  }



  # Setting parallel or sequential backend and  random seed
  if (ParExec & (parallel::detectCores() > 2 & m > 1)) {
    if (is.null(ncores)) {
      Ncpus <- min(m, parallel::detectCores() - 1)
    } else {
      Ncpus <- min(ncores, parallel::detectCores() - 1)
    }
    cl <- parallel::makeCluster(Ncpus)
    doSNOW::registerDoSNOW(cl) 
    if (SetRNGSeed) {
      doRNG::registerDoRNG(SetRNGSeed)
    }
    # set progress bar for parallel processing
    pb <- txtProgressBar(max = m, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    # condition used to run code part needed for parallel processing
    ParParams <- TRUE
  } else {
    if (ParExec & m == 1) {
      if (verbose == TRUE) {
        message(paste("/!\\ The number of multiple imputations is 1, 
          parallel processing is only available for m > 1."))
      }
    } else if (ParExec) {
      if (verbose == TRUE) {
        message(paste("/!\\ The number of cores of your processor does not 
        allow paralell processing, at least 3 cores are needed."))
      }
    }

    if (SetRNGSeed) {
      set.seed(SetRNGSeed)
    }

    foreach::registerDoSEQ()
    opts <- NULL

    # condition used to run code part needed for sequential processing
    ParParams <- FALSE
  }


  # Beginning of the multiple imputation (imputing "mi" times)
  o <- NULL
  RESULT <- foreach(o = 1:m, .inorder = TRUE, 
    .options.snow = opts) %dopar% {
    if (!ParParams) {
      if (verbose == TRUE) {
        cat("iteration :", o, "/", m, "\n")
      }
    }
    # Trying to catch the potential singularity error (part 1/2)
    # (comment this part of code (as well as the one at the very end of
    # seqimpute.R) to debug more easily and see the full message of any
    # occuring error)
    
    # 3. Imputation using a specific model ---------------
    if (max(dataOD$ORDER) != 0) {
      # Otherwise if there is only 0 in ORDER,
      # there is no need to impute internal gaps
      # and we directly jump to the imputation of
      # external gaps (i.e. points 4. and 5.)

      if (verbose == TRUE) {
        print("Imputation of the internal gaps...")
      }

      dataOD[["ODi"]] <- ModelImputationTiming(dataOD$OD, dataOD$CO, 
        dataOD$COt, dataOD$ODi, dataOD$MaxGap, dataOD$totV, dataOD$totVi, 
        regr, dataOD$nc, np, nf, dataOD$nr, dataOD$ncot, dataOD$COtsample, 
        dataOD$pastDistrib, dataOD$futureDistrib, dataOD$k, available, 
        dataOD$REFORD_L, dataOD$noise, frame.radius, verbose, ...)
    }
    # 4. Imputing initial NAs -------------------------------------------
    if ((nfi != 0) & (dataOD$MaxInitGapSize != 0)) {
      if (verbose == TRUE) {
        print("Imputation of the initial gaps...")
      }
      # # we only impute the initial gaps if nfi > 0
      dataOD[["ODi"]] <- ImputingInitialNAsTiming(
        OD = dataOD$OD, covariates = dataOD$CO, time.covariates = dataOD$COt, 
        ODi = dataOD$ODi, COtsample = dataOD$COtsample,
        futureDistrib = dataOD$futureDistrib, InitGapSize = dataOD$InitGapSize, 
        MaxInitGapSize = dataOD$MaxInitGapSize, nr = dataOD$nr, nc = dataOD$nc,
        ud = dataOD$ud, nco = dataOD$nco, ncot = dataOD$ncot, nfi = nfi,
        regr = regr, k = dataOD$k, available = available, noise = dataOD$noise, 
        frame.radius = frame.radius, ...)
    }
    # 5. Imputing terminal NAs --------------------------------------------
    if ((npt != 0) & (dataOD$MaxTermGapSize != 0)) {
      # we only impute the terminal
      # gaps if npt > 0
      if (verbose == TRUE) {
        print("Imputation of the terminal gaps...")
      }
      dataOD[["ODi"]] <- ImputingTerminalNAsTiming(
        OD = dataOD$OD, covariates = dataOD$CO, time.covariates = dataOD$COt, 
        ODi = dataOD$ODi, COtsample = dataOD$COtsample,
        pastDistrib = dataOD$pastDistrib, TermGapSize = dataOD$TermGapSize, 
        MaxTermGapSize = dataOD$MaxTermGapSize, nr = dataOD$nr, nc = dataOD$nc,
        ud = dataOD$ud, nco = dataOD$nco, ncot = dataOD$ncot, npt = npt, 
        regr = regr, k = dataOD$k, available = available, noise = dataOD$noise, 
        frame.radius = frame.radius, ...)
    }
    if (max(dataOD$ORDERSLGLeft) != 0 & !is.null(dataOD$ORDERSLGLeft)) {
      # Checking if we have to impute
      # left-hand side SLG
      if (verbose == TRUE) {
        print("Imputation of the left-hand side SLG...")
      }
      dataOD[["ODi"]] <- LSLGNAsImputeTiming(dataOD$OD, dataOD$ODi, dataOD$CO, 
        dataOD$COt, dataOD$COtsample, dataOD$ORDERSLGLeft, dataOD$pastDistrib, 
        dataOD$futureDistrib, regr, np, dataOD$nr, nf, dataOD$nc, dataOD$ud, 
        dataOD$ncot, dataOD$nco, dataOD$k, dataOD$noise, available, 
        frame.radius, ...)
    }
    # right-hand side SLG
    if (max(dataOD$ORDERSLGRight) != 0 & !is.null(dataOD$ORDERSLGRight)) {
      # Checking if we have to impute right-hand
      # side SLG
      if (verbose == TRUE) {
        print("Imputation of the right-hand side SLG...")
      }

      dataOD[["ODi"]] <- RSLGNAsImputeTiming(dataOD$OD, dataOD$ODi, dataOD$CO, 
        dataOD$COt, dataOD$COtsample, dataOD$ORDERSLGRight, dataOD$pastDistrib, 
        dataOD$futureDistrib, regr, np, dataOD$nr, nf, dataOD$nc, dataOD$ud, 
        dataOD$ncot, dataOD$nco, dataOD$k, dataOD$noise, 
        available, frame.radius, ...)
    }
    # Checking if we have to impute
    # Both-hand side SLG
    if (dataOD$LongGap) {
      if (verbose == TRUE) {
        print("Imputation of the both-hand side SLG...")
      }
      for (h in 2:np) {
        if (sum(dataOD$ORDERSLGBoth[, h - 1] == 0 & 
          dataOD$ORDERSLGBoth[, h] != 0) > 0) {
          tt <- which(dataOD$ORDERSLGBoth[, h - 1] == 0 & 
            dataOD$ORDERSLGBoth[, h] != 0)
          tmpORDER <- matrix(0, nrow(dataOD$ORDERSLGBoth), 
            ncol(dataOD$ORDERSLGBoth))
          tmpORDER[tt, h:ncol(dataOD$ORDERSLGBoth)] <- dataOD$ORDERSLGBoth[tt, 
            h:ncol(dataOD$ORDERSLGBoth)]
          dataOD[["ODi"]] <- RSLGNAsImputeTiming(dataOD$OD, dataOD$ODi, 
            dataOD$CO, dataOD$COt, dataOD$COtsample, tmpORDER, 
            dataOD$pastDistrib, dataOD$futureDistrib, regr, h - 1, dataOD$nr, 
            nf, dataOD$nc, dataOD$ud, dataOD$ncot, dataOD$nco, dataOD$k, 
            dataOD$noise, available, frame.radius, ...)
        }
      }
    }


    # Updating the matrix RESULT used to store the multiple imputations
      return(dataOD$ODi)
    }
  if (ParParams) {
    parallel::stopCluster(cl)
  }
  
  names(RESULT) <- paste0("imp",1:m)
  
  # RESULT <- rbind(cbind(replicate(dataOD$nr, 0), dataOD$OD), RESULT)
  
  # X. Final conversions -----------------------------------------------------
  RESULT <- lapply(RESULT,FinalResultConvert, ODClass = dataOD$ODClass,
      ODlevels = dataOD$ODlevels, rownamesDataset = rownamesDataset, 
      nrowsDataset = nrowsDataset, nr = dataOD$nr, nc = dataOD$nc, 
      rowsNA = dataOD$rowsNA, mi = m)
  
  return(RESULT)
  
  
}

ModelImputationTiming <- function(OD, covariates, time.covariates, ODi, MaxGap,
  totV, totVi, regr, nc, np, nf, nr, ncot, COtsample, pastDistrib, 
  futureDistrib, k, available, REFORD_L, noise, frame.radius, verbose, ...)
{
  for (order in 1:MaxGap) {
    if(verbose == TRUE){
      print(paste0("Step ", order, "/", MaxGap))
    }
    ncol_imp <- length(unique(REFORD_L[[order]][, 2]))
    col_to_imp <- unique(sort(unique(REFORD_L[[order]])[, 2]))
    for (i in 1:ncol_imp) {
      CD_shift <- CDComputeTiming(covariates, OD, time.covariates, MaxGap, 
        order, np, nc, nr, nf, COtsample, pastDistrib, futureDistrib, ncot, k, 
        col_to_imp[i], frame.radius)
      if (length(table(CD_shift$CD[, 1])) > 1) {
        log_CD <- list()
        log_CD[c("reglog", "CD")] <- ComputeModel(CD_shift$CD, regr, 
          totV, np, nf, k, ...)
        row_to_imp <- REFORD_L[[order]][
          which(REFORD_L[[order]][, 2] == col_to_imp[[i]]), 1]
        ODi <- CreatedModelImputationTiming(order, covariates, log_CD$CD, 
          time.covariates, OD, ODi, pastDistrib, futureDistrib, available, 
          col_to_imp[i], row_to_imp, ncot, nc, np, nf, k, totV, regr, 
          log_CD$reglog, noise, CD_shift$shift, MaxGap)
      } else {
        lev <- names(table(CD_shift$CD[, 1]))
        REFORD <- as.matrix(REFORD_L[[order]])
        if (ncol(REFORD) == 1) {
          REFORD <- t(REFORD)
        }
        nr_REFORD <- nrow(REFORD)

        for (u in 1:nr_REFORD) {
          i <- REFORD[u, 1]
          # taking out the first coordinate (row
          # number in ORDER) from REFORDI
          j <- REFORD[u, 2]
          ODi[i, j] <- lev
        }
      }
    }
  }

  return(ODi)
}


CDComputeTiming <- function(CO, OD, COt, MaxGap, order, np, nc, nr, nf, 
  COtsample, pastDistrib, futureDistrib, ncot, k, col, frame.radius)
{
  # Building of a data matrix for the computation of the
  # model
  # number of usable data in past and futur


  # for each row of

  # frameSize <- MaxGap-order+np+nf+1 # size of the current

  # mobile caracteristic frame (that
  # changes according to
  # "order") which is equal to
  # nc-ud+1
  # Structure and building of the data matrix CD
  # The first column of CD is the dependent variable (VD,
  # response variable)
  # The following columns are the independent variables
  # (VIs, explanatory variables) coming from the past
  # (np>0) or the future (nf>0) ordered by time and the
  # distribution of the possible values (i.e. all possible
  # categorical variables numbered from 1 to k and of
  # course the value NA) respectively Before and After the
  # NA to impute.
  #

  # We are then going to create a specific CD according to
  # the value of np and nf
  # initialization of
  # the current very left part of
  # the predictive model matrix
  # ("matrice de modele de
  # prediction") with NA
  # everywhere (/!\ for each
  # "order", we are going to
  # build such a CD)

  if ((np > 0 & nf > 0) & ((MaxGap %% 2 == 0 & order %% 2 == 0) | 
    (MaxGap %% 2 != 0 & order %% 2 != 0))) {
    shift <- MaxGap - order # jumping at the end of
    udp <- min(frame.radius, col - (MaxGap - order) - np - 1) #data in the past
    udf <- min(frame.radius, nc - col - nf) # number of usable data in the futur
    col_to_use <- (col - udp):(col + udf)
    ud <- udp + udf + 1
    # the gap
  } else {
    shift <- 0 # no shift is needed (note that
    # no shift is needed for the
    # case of model 2 (only past)
    # and model 3 (only future))

    if (np > 0 & nf > 0) {
      udp <- min(frame.radius, col - np - 1)
      udf <- min(frame.radius, nc - col - (MaxGap - order) - nf)
      col_to_use <- (col - udp):(col + udf)
      ud <- udp + udf + 1
    }

    if (np > 0 & nf == 0) {
      udp <- min(frame.radius, col - np - 1)
      udf <- min(frame.radius, nc - col)
      col_to_use <- (col - udp):(col + udf)
      ud <- udp + udf + 1
    }

    if (nf > 0 & np == 0) {
      udf <- min(frame.radius, nc - col - nf)
      udp <- min(col - 1, frame.radius)
      col_to_use <- (col - udp):(col + udf)
      ud <- udp + udf + 1
    }
  }


  CD <- matrix(NA, nrow = nr * ud, ncol = 1)


  iter <- 1 # initialisation of the number
  # of iterations of the
  # following for loops
  # Only PAST
  if (np > 0 & nf == 0) {
    CD <- PastVIComputeTiming(CD, CO, OD, ncot, iter, nr, nc, ud, np, 
      COtsample, COt, pastDistrib, futureDistrib, k, col_to_use)
    # Only FUTURE
  } else if (np == 0 & nf > 0) {
    # only FUTURE VIs do exist
    CD <- FutureVIComputeTiming(CD, CO, OD, ncot, iter, nr, nc, ud, np, 
      COtsample, COt, pastDistrib, futureDistrib, k, nf, col_to_use)

    # PAST and FUTURE
  } else {
    # meaning np>0 and nf>0 and that, thus,
    # PAST as well as FUTURE VIs do exist
    
  
    CD <- PastFutureVIComputeTiming(CD, CO, OD, ncot, iter, nr, nc, ud, np, 
    COtsample, COt, pastDistrib, futureDistrib, k, nf, shift, col_to_use, 
    MaxGap, order)
  }

  CD_shift <- list()
  CD_shift[c("CD", "shift")] <- list(CD, shift)
  return(CD_shift)
}

################################################################################
FutureVIComputeTiming <- function(CD, CO, OD, ncot, iter, nr, nc, ud, np, 
  COtsample, COt, pastDistrib, futureDistrib, k, nf, col_to_use) {
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
  for (j in col_to_use) {
    t1 <- (nr * (iter - 1) + 1)
    # Determining the locations
    # of the time span (always nr) of
    # the piled up blocks of CD
    t2 <- nr * iter
    # VD
    CD[t1:t2, 1] <- OD[, j]

    # future VIs
    CDf[t1:t2, ] <- OD[, (j + 1):(j + nf)]

    # /!\ current
    # pointer on column
    # is thus:
    # "j-frameSize

    # Eventually considering time-dependent
    # covariates


    if (ncot > 0) {
      COtselected[t1:t2,] <- COt[t1:t2, j + (1:(ncot / nc) - 1) * nc]
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
    # Future distribution (i.e. After)
    if (futureDistrib) {
      ODt <- t(OD)
      ODt <- as.data.frame(ODt)
      tempOD <- lapply(ODt[((j + 1):nc), ], factor, levels = c(1:k, NA), 
        exclude = NULL)
      # because:
      # j-frameSize+np+1 + 1
      # = j-frameSize+np+2

      da_list <- lapply(tempOD, summary)
      da_matrix <- do.call(rbind, da_list)
      CDda[t1:t2, ] <- da_matrix[, 1:k] / length((j + 1):nc)
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


PastVIComputeTiming <- function(CD, CO, OD, ncot, iter, nr, nc, ud, np, 
  COtsample, COt, pastDistrib, futureDistrib, k, col_to_use) {
  CDp <- matrix(NA, nrow = nr * ud, ncol = np)

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
  for (j in col_to_use) {
    # /!\ j is initialised at
    # the very end (utmost right) of the
    # frame
    t1 <- (nr * (iter - 1) + 1)
    # Determining the locations
    # of the time span (always nr) of
    # the piled up blocks of CD
    t2 <- nr * iter
    # VD
    CD[t1:t2, 1] <- OD[, j]


    # past VIs
    CDp[t1:t2, ] <- OD[, (j - np):(j - 1)]


    # /!\ current
    # pointer on column
    # is thus:
    # "j-frameSize

    # Eventually considering time-dependent
    # covariates
    if (ncot > 0) {
      COtselected[t1:t2,] <- COt[t1:t2, j + (1:(ncot / nc) - 1) * nc]
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
    # Future distribution (i.e. After)
    if (futureDistrib) {
      ODt <- t(OD)
      ODt <- as.data.frame(ODt)
      tempOD <- lapply(ODt[((j + 1):nc), ], factor, levels = c(1:k, NA), 
        exclude = NULL)
      # because:
      # j-frameSize+np+1 + 1
      # = j-frameSize+np+2

      da_list <- lapply(tempOD, summary)
      da_matrix <- do.call(rbind, da_list)
      CDda[t1:t2, ] <- da_matrix[, 1:k] / length((j + 1):nc)
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

# Past or future VI computation

PastFutureVIComputeTiming <- function(CD, CO, OD, ncot, iter, nr, nc, ud, np, 
  COtsample, COt, pastDistrib, futureDistrib, k, nf, shift, col_to_use, 
  MaxGap, order)
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
  
  for (j in col_to_use) {
    # /!\ j is initialised at
    # the very end (utmost right) of the
    # frame
    t1 <- (nr * (iter - 1) + 1)
    # Determining the locations
    # of the time span (always nr) of
    # the piled up blocks of CD
    t2 <- nr * iter
    # VD
    CD[t1:t2, 1] <- OD[, j]


    if (shift == 0) {
      CDp[t1:t2, ] <- OD[, (j - np):(j - 1)]
      CDf[t1:t2, ] <- OD[, (j + MaxGap - order + 1):(j + MaxGap - order + nf)]
    } else {
      CDp[t1:t2, ] <- OD[, (j - shift - np):(j - shift - 1)]
      CDf[t1:t2, ] <- OD[, (j + 1):(j + nf)]
    }
    #
    # # past VIs
    # CDp[t1:t2,] <- OD[,(j-frameSize+1):(j-frameSize+np)]
    #
    # # future VIs
    # CDf[t1:t2,] <- OD[,(j-nf+1):j]
    #
    # /!\ current
    # pointer on column
    # is thus:
    # "j-frameSize

    # Eventually considering time-dependent
    # covariates
    
    
    if (ncot > 0) {
      COtselected[t1:t2,] <- COt[t1:t2, j + (1:(ncot / nc) - 1) * nc]
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
    # Future distribution (i.e. After)
    if (futureDistrib) {
      ODt <- t(OD)
      ODt <- as.data.frame(ODt)
      tempOD <- lapply(ODt[((j + 1):nc), ], factor, levels = c(1:k, NA), 
        exclude = NULL)
      # because:
      # j-frameSize+np+1 + 1
      # = j-frameSize+np+2

      da_list <- lapply(tempOD, summary)
      da_matrix <- do.call(rbind, da_list)
      CDda[t1:t2, ] <- da_matrix[, 1:k] / length((j + 1):nc)
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


################################################################################
# Imputation using the just created model 

CreatedModelImputationTiming <- function(order, CO, CD, COt, OD, ODi, 
  pastDistrib, futureDistrib, available, col, row_to_imp, ncot, nc, np, nf, k, 
  totV, regr, reglog, noise, shift, MaxGap) {
  # Structure and building of the data matrix CDi
  # The first column of CDi is the dependent variable (VD,
  # response variable) that we have to implement during
  # the current iteration (i.e. automatically a NA)
  # The following columns are the corresponding
  # independent variables (VIs, explanatory variables)
  # coming from the past (np>0) or the future (nf>0)
  # (ordered by time and corresponding to the current
  # predictive pattern) and the distribution of the
  # possible values (i.e. all possible categorical
  # variables numbered from 1 to k and of course
  # the value NA) respectively Before and After the NA to
  # impute.
  # (The fact that every lines of CDi are identical is
  # related to the working of the function mlogit that has
  # to have as much lines in CDi as there are categories
  # of the variable (see the parameter "k") --> so, CDi is
  # composed of k identical lines)
  #

  # We are then going to create a specific CDi according
  # to the value of np and nf


  # Analysing the value of parameter available
  if (available == TRUE) { # we take the previously imputed
    # data into account
    LOOKUP <- ODi
  } else { # that is available == FALSE and thus we
    # don't take the previously imputed
    # data into account
    LOOKUP <- OD
  }


  # Assigning the current "REFORD_order" matrix to the
  # variable matrix REFORD
  # (according to the current value of "order")
  # REFORD <- as.matrix(REFORD_L[[order]])
  # if (ncol(REFORD) == 1) {
  #   REFORD <- t(REFORD)
  # }
  # nr_REFORD <- nrow(REFORD)


  if (np > 0 & nf == 0) { # only PAST VIs do exist
    ODi <- ODiImputePASTTiming(CO, ODi, CD, COt, col, row_to_imp, pastDistrib, 
      futureDistrib, k, np, nf, nc, ncot, totV, reglog, LOOKUP, regr, noise)

    
  } else if (np == 0 & nf > 0) { # only FUTURE VIs do exist
    ODi <- ODiImputeFUTURETiming(CO, ODi, CD, COt, col, row_to_imp, 
      pastDistrib, futureDistrib, k, np, nf, nc, ncot, totV, reglog, 
      LOOKUP, regr, noise)
  } else { # meaning np>0 and nf>0 and that,
    # thus, PAST as well as FUTURE VIs
    # do exist
    ODi <- ODiImputePFTiming(CO, ODi, CD, COt, col, row_to_imp, pastDistrib, 
      futureDistrib, k, np, nf, nc, ncot, totV, reglog, LOOKUP, regr, noise, 
      shift, MaxGap, order)
  }
  return(ODi)
}

ODiImputePASTTiming <- function(CO, ODi, CD, COt, col, row_to_imp, pastDistrib, 
  futureDistrib, k, np, nf, nc, ncot, totV, reglog, LOOKUP, regr, noise) {
  for (u in 1:length(row_to_imp)) {
    i <- row_to_imp[u]
    # taking out the first coordinate
    # (row number in ORDER) from REFORD
    j <- col
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

    # We also account for the fact that levels that do not appear at
    # all in a given variable of CD were discarded with droplevels before
    # the fit of the mlogit model


    if (regr != "rf") {
      for (v in 1:(1 + np + nf)) {
        CDi[, v] <- factor(CDi[, v], levels = levels(CD[, v]), exclude = NULL)
      }
    } else {
      for (v in 2:(1 + np + nf)) {
        CDi[, v] <- factor(CDi[, v], levels = c(1:(k + 1)))
        CDi[, v][is.na(CDi[, v])] <- k + 1
      }
      CDi[, 1] <- factor(CDi[, 1], levels = levels(CD[, 1]))
    }
    # The last values of CDi must be of type numeric
    # (distributions)
    if (pastDistrib | futureDistrib) {
      CDi[, (1 + np + nf + 1):ncol(CDi)] <- lapply(CDi[, 
        (1 + np + nf + 1):ncol(CDi)], as.numeric)
    }

    # Eventually concatenating CDi with COi
    # (the matrix containing the covariates)
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
      # Concatenating CDi and COtselected_i
      # into CDi
      CDi <- cbind(CDi, COtselected_i)
      # Transformation of the names of the columns
      # of CDi (called V1, V2, ..., "VtotV")
      colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
    }


    # Check for missing-values among predictors
    if (max(is.na(CDi[1, 2:totV])) == 0) {
      ODi <- RegrImpute(ODi, CDi, regr, reglog, noise, i, j, k)
    }
  }
  return(ODi)
}


ODiImputeFUTURETiming <- function(CO, ODi, CD, COt, col, row_to_imp, 
  pastDistrib, futureDistrib, k, np, nf, nc, ncot, totV, reglog, 
  LOOKUP, regr, noise) {
  for (u in 1:length(row_to_imp)) {
    i <- row_to_imp[u]
    # taking out the first coordinate
    # (row number in ORDER) from REFORD
    j <- col
    # taking out the second coordinate
    # (column number in ORDER) from REFORD
    CDi <- matrix(NA, nrow = k, ncol = 1)


    # Matrix for future values
    vect <- LOOKUP[i, (j + 1):(j + nf)]
    CDfi <- matrix(vect, nrow = k, ncol = length(vect), byrow = TRUE)

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
    CDi <- cbind(CDi, CDfi)

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
    # We also account for the fact that levels that do not appear at
    # all in a given variable of CD were discarded with droplevels before
    # the fit of the mlogit model

    if (regr != "rf") {
      for (v in 1:(1 + np + nf)) {
        CDi[, v] <- factor(CDi[, v], levels = levels(CD[, v]), exclude = NULL)
      }
    } else {
      for (v in 2:(1 + np + nf)) {
        CDi[, v] <- factor(CDi[, v], levels = c(1:(k + 1)))
        CDi[, v][is.na(CDi[, v])] <- k + 1
      }
      CDi[, 1] <- factor(CDi[, 1], levels = levels(CD[, 1]))
    }
    # The last values of CDi must be of type numeric
    # (distributions)
    if (pastDistrib | futureDistrib) {
      CDi[, (1 + np + nf + 1):ncol(CDi)] <- lapply(CDi[,
        (1 + np + nf + 1):ncol(CDi)], as.numeric)
    }

    # Eventually concatenating CDi with COi
    # (the matrix containing the covariates)
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
    if (max(is.na(CDi[1, 2:totV])) == 0) {
      ODi <- RegrImpute(ODi, CDi, regr, reglog, noise, i, j, k)
    }
  }
  return(ODi)
}


ODiImputePFTiming <- function(CO, ODi, CD, COt, col, row_to_imp, pastDistrib, 
  futureDistrib, k, np, nf, nc, ncot, totV, reglog, LOOKUP, regr, noise, shift, 
  MaxGap, order) {
  for (u in 1:length(row_to_imp)) {
    i <- row_to_imp[u]
    # taking out the first coordinate
    # (row number in ORDER) from REFORD
    j <- col
    # taking out the second coordinate
    # (column number in ORDER) from REFORD
    CDi <- matrix(NA, nrow = k, ncol = 1)

    # Matrix for past values
    shift <- as.numeric(shift)

    if (shift == 0) {
      vect <- LOOKUP[i, (j - np):(j - 1)]
    } else {
      vect <- LOOKUP[i, (j - shift - np):(j - shift - 1)]
    }

    CDpi <- matrix(vect, nrow = k, ncol = length(vect), byrow = TRUE)

    # Matrix for future values
    if (shift == 0) {
      vect <- LOOKUP[i, (j + MaxGap - order + 1):(j + MaxGap - order + nf)]
    } else {
      vect <- LOOKUP[i, (j + 1):(j + nf)]
    }
    CDfi <- matrix(vect, nrow = k, ncol = length(vect), byrow = TRUE)

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
    CDi <- cbind(CDi, CDpi, CDfi)

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
    # We also account for the fact that levels that do not appear at
    # all in a given variable of CD were discarded with droplevels before
    # the fit of the mlogit model

    # if (regr == "lm" | regr == "lrm" | regr == "mlogit") {
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
    #   # glmnet
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
      CDi[, (1 + np + nf + 1):ncol(CDi)] <- lapply(CDi[, 
        (1 + np + nf + 1):ncol(CDi)], as.numeric)
    }

    # Eventually concatenating CDi with COi
    # (the matrix containing the covariates)
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
    if (max(is.na(CDi[1, 2:totV])) == 0) {
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


# 4. Imputation of initial NAs
#
#
###############################################################################
# Impute initial NAs
ImputingInitialNAsTiming <- function(OD, covariates, time.covariates, ODi,
  COtsample, futureDistrib, InitGapSize, MaxInitGapSize, nr, nc, ud,
  nco, ncot, nfi, regr, k, available, noise, frame.radius, ...)
{
  # 4.1.-2. Creation of ORDERI -------------------------------------------------
  REFORDI_L <- REFORDICreation(nr, nc, InitGapSize, MaxInitGapSize)

  # 4.3. Imputation using a specific model -------------------------------------

  for (order in 1:MaxInitGapSize) {
    col <- MaxInitGapSize - order + 1
    CD <- CDMatCreateTiming(covariates, COtsample, OD, time.covariates, 
      min(nfi, nc - col), nr, nc, ncot, futureDistrib, k, col, frame.radius)

    # 4.3.2 Computation of the model (Dealing with the LOCATIONS of imputation)
    if (length(table(CD[, 1])) > 1) {
      log_CD <- list()

      totVi <- 1 + min(nfi, nc - col) + nco + (ncot / nc)
      log_CD[c("reglog", "CD")] <- ComputeModel(CD, regr, totVi, 0, 
        min(nfi, nc - col), k, ...)

      # 4.3.3 Imputation using the just created model 
      ODi <- Init_NA_CreatedModelImputationTiming(OD, ODi, covariates, 
        log_CD$CD, time.covariates, MaxInitGapSize, REFORDI_L, futureDistrib, 
        totVi, nc, k, min(nfi, nc - col), ncot, regr, log_CD$reglog, 
        noise, available, order)
    } else {
      lev <- names(table(CD[, 1]))
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
        ODi[i, j] <- lev
      }
    }
  }


  return(ODi)
}

Init_NA_CreatedModelImputationTiming <- function(OD, ODi, CO, CD, COt, 
  MaxInitGapSize, REFORDI_L, futureDistrib, totVi, nc, k, nfi, ncot, regr, 
  reglog, noise, available, order)
{
  # Conversion of ODi from data.frame to matrix
  ODi <- as.matrix(ODi)

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
    # We also account for the fact that levels that do not appear at
    # all in a given variable of CD were discarded with droplevels before
    # the fit of the mlogit model
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
      CDi[, (1 + nfi + 1):ncol(CDi)] <- lapply(CDi[, (1 + nfi + 1):ncol(CDi)], 
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
  return(ODi)
}



CDMatCreateTiming <- function(CO, COtsample, OD, COt, nfi, nr, nc, ncot, 
  futureDistrib, k, col, frame.radius) {
  # For Initial Gaps
  # we will impute single NA after single NA going from the
  # center of OD towards its very left border
  # We here only take observations from the FUTURE into
  # account --> nfi
  # But we will create one single regression
  # model used for the imputation of all
  # the NAs belonging to an "initial gap"
  udf <- min(frame.radius, nc - col - nfi)
  udp <- min(frame.radius, col - 1)
  ud <- udf + udp + 1
  col_to_use <- (col - udp):(col + udf)

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

  for (j in col_to_use) {
    t1 <- (nr * (iter - 1) + 1)
    t2 <- nr * iter

    # VD
    CD[t1:t2, 1] <- OD[, j]

    # Future VIs
    CDf[t1:t2, ] <- OD[, (j + 1):(j + nfi)]

    # Eventually considering time-dependent covariates
    if (ncot > 0) {
      COtselected[t1:t2,] <- COt[t1:t2, j + (1:(ncot / nc) - 1) * nc]
    }
    # Future distribution (i.e. After)
    if (futureDistrib) {
      ODt <- t(OD)
      ODt <- as.data.frame(ODt)
      tempOD <- lapply(ODt[((j + 1):nc), ], factor, levels = c(1:k, NA), 
        exclude = NULL)
      # because:
      # j-frameSize+np+1 + 1
      # = j-frameSize+np+2

      da_list <- lapply(tempOD, summary)
      da_matrix <- do.call(rbind, da_list)
      CDda[t1:t2, ] <- da_matrix[, 1:k] / length((j + 1):nc)
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


# 5. Imputation of terminal NAs
#
#
################################################################################
# Impute terminal NAs

ImputingTerminalNAsTiming <- function(OD, covariates, time.covariates, ODi, 
  COtsample, MaxTermGapSize, TermGapSize, pastDistrib, regr, npt, nco, ncot, 
  totVt, nr, nc, ud, available, k, noise, frame.radius, ...) {
  # 5.1.-2. Creation of ORDERT -------------------------------------------------
  REFORDT_L <- REFORDTCreation(nr, nc, TermGapSize, MaxTermGapSize)

  for (order in 1:MaxTermGapSize) {
    col <- nc - MaxTermGapSize + order
    CD <- TerminalCDMatCreateTiming(covariates, OD, time.covariates, COtsample, 
      pastDistrib, min(npt, col - 1), nr, nc, ncot, k, col, frame.radius)
    # 5.3.2 Computation of the model (Dealing with the LOCATIONS of imputation)
    if (length(table(CD[, 1])) > 1) {
      log_CD <- list()
      totVt <- 1 + min(npt, col - 1) + nco + (ncot / nc)
      log_CD[c("reglog", "CD")] <- ComputeModel(CD, regr, totVt, 
        min(npt, col - 1), 0, k, ...)


      # 5.3.3 Imputation using the just created model 
      ODi <- TerminalCreatedModelImputationTiming(covariates, OD, log_CD$CD, 
        ODi, time.covariates, nc, ncot, totVt, REFORDT_L, pastDistrib, 
        MaxTermGapSize, available, regr, log_CD$reglog, k, min(npt, col - 1), 
        noise, order)
    } else {
      lev <- names(table(CD[, 1]))
      REFORDT <- as.matrix(REFORDT_L[[order]])
      if (ncol(REFORDT) == 1) {
        REFORDT <- t(REFORDT)
      }
      nr_REFORDT <- nrow(REFORDT)

      for (u in 1:nr_REFORDT) {
        i <- REFORDT[u, 1]
        # taking out the first coordinate (row
        # number in ORDER) from REFORDT
        j <- REFORDT[u, 2]
        ODi[i, j] <- lev
      }
    }
  }

  # 5.3. Imputation using a specific model -------------------------------------

  # 5.3.1 Building of the data matrix CD for the computation of the model ------

  return(ODi)
}



###############################################################################
# Imputation of the terminal created model

TerminalCreatedModelImputationTiming <- function(CO, OD, CD, ODi, COt, nc, 
  ncot, totVt, REFORDT_L, pastDistrib, MaxTermGapSize, available, regr, 
  reglog, k, npt, noise, order)
{
  # Conversion of ODi from data.frame to matrix
  ODi <- as.matrix(ODi)

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
    # We also account for the fact that levels that do not appear at
    # all in a given variable of CD were discarded with droplevels before
    # the fit of the mlogit model
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
      CDi[, (1 + npt + 1):ncol(CDi)] <- lapply(CDi[, (1 + npt + 1):ncol(CDi)], 
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
  return(ODi)
}


TerminalCDMatCreateTiming <- function(CO, OD, COt, COtsample, pastDistrib, npt, 
  nr, nc, ncot, k, col, frame.radius) {
  # For Terminal Gaps
  # we will impute single NA after single NA going from the
  # center of OD towards its very right border
  # We here only take observations from the PAST
  # into account --> npt
  # But we will create one single regression
  # model used for the imputation of all the NAs belonging to
  # a "terminal gap"
  udf <- min(frame.radius, nc - col)
  udp <- min(frame.radius, col - 1 - npt)
  ud <- udf + udp + 1
  col_to_use <- (col - udp):(col + udf)
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
  for (j in col_to_use) {
    t1 <- (nr * (iter - 1) + 1)
    t2 <- nr * iter
    # VD
    CD[t1:t2, 1] <- OD[, j]
    # /!\ current pointer on column is thus: "j"

    # Past VIs
    CDp[t1:t2, ] <- OD[, (max(j - npt, 1)):(j - 1)]

    # Eventually considering time-dependent covariates
    if (ncot > 0) {
      COtselected[t1:t2,] <- COt[t1:t2, j + (1:(ncot / nc) - 1) * nc]
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


# 6. Imputation of SLG NAs
#
#
################################################################################
# Left-hand side SLG imputation

LSLGNAsImputeTiming <- function(OD, ODi, CO, COt, COtsample, ORDERSLG, 
  pastDistrib, futureDistrib, regr, np, nr, nf, nc, ud, ncot, nco, k, noise, 
  available, frame.radius, ...)
  { # Checking if we have to impute
  # left-hand side SLG


  # 6.2.LEFT Computation of the order of imputation of each MD ----

  # Initialization of a list to take all the variable returned by the functions
  ParamList <- list()
  # Creation of the temporary SLG matrices for the left-hand
  # side of OD

  for (h in 2:np) {
    if (max(ORDERSLG[, h]) > 0) {
      # Checking if a gap begins
      # somewhere on the current column
      # If that is not the case, there is
      # no need to perform
      # the entire following procedure
      # and we simply can go
      # to the next column of ORDERSLGLeft

      ParamList[c("ORDERSLG_temp","totV_temp","np_temp")] <- SLGMatrix_temp(nr, 
        nc, nf, h, ORDERSLG, nco, ncot, pastDistrib, futureDistrib, k)

      if (max(ParamList$ORDERSLG_temp) == 0) {
        next
      }

      # In a similar manner to part 2.4., we go here one
      # single time through the reduced version
      # ORDERSLGLeft_temp of ORDERSLG and we create
      # MaxGapSLGLeft "REFORDSLGRLeft_" matrices
      # collecting the coordinates of each corresponding
      # values in ORDERSLGLeft_temp which are greater
      # than 0


      # REFORDSLGLeft matrices
      # Initialization of the REFORDSLGLeft matrices
      ParamList[c("MaxGap", "REFORD_L", 
        "ORDERSLG_temp")] <- REFORDInit(ParamList$ORDERSLG_temp, nr, nc)

     
      # 6.3.LEFT Imputation of the missing data listed by ORDERSLGLeft_temp a

      # 6.3.1.LEFT Building of the different data matrices CD ----
      # for the computation of the model for every SLG
      # on the left-hand side of OD
      for (order in 1:ParamList$MaxGap) {
        ncol_imp <- length(unique(ParamList$REFORD_L[[order]][, 2]))
        col_to_imp <- unique(sort(unique(ParamList$REFORD_L[[order]])[, 2]))
        for (i in 1:ncol_imp) {
          CD_shift <- CDComputeTiming(CO, OD, COt, ParamList$MaxGap, order, 
            ParamList$np_temp, nc, nr, nf, COtsample, pastDistrib, 
            futureDistrib, ncot, k, col_to_imp[i], frame.radius)

          if (length(table(CD_shift$CD[, 1])) > 1) {
            log_CD <- list()
            log_CD[c("reglog", "CD")] <- ComputeModel(CD_shift$CD, regr, 
                ParamList$totV_temp, ParamList$np_temp, nf, k, ...)

            row_to_imp <- ParamList$REFORD_L[[order]][
              which(ParamList$REFORD_L[[order]][, 2] == col_to_imp[i]), 1]
            # 3.3. Imputation using the just created model 
            ODi <- CreatedModelImputationTiming(order, CO, log_CD$CD, COt, OD, 
              ODi, pastDistrib, futureDistrib, available, col_to_imp[i], 
              row_to_imp, ncot, nc, ParamList$np_temp, nf, k, 
              ParamList$totV_temp, regr, log_CD$reglog, noise, 
              CD_shift$shift, ParamList$MaxGap)
          } else {
            lev <- names(table(CD_shift$CD[, 1]))
            REFORDI <- as.matrix(ParamList$REFORD_L[[order]])
            if (ncol(REFORDI) == 1) {
              REFORDI <- t(REFORDI)
            }
            nr_REFORDI <- nrow(REFORDI)

            for (u in 1:nr_REFORDI) {
              i <- REFORDI[u, 1]
              # taking out the first coordinate (row
              # number in ORDER) from REFORDI
              j <- REFORDI[u, 2]
              ODi[i, j] <- lev
            }
          }
        }
      }
    }
  }
  return(ODi)
}

##############################################################################
# Right-hand side SLG imputation

RSLGNAsImputeTiming <- function(OD, ODi, CO, COt, COtsample, ORDERSLGRight, 
  pastDistrib, futureDistrib, regr, np, nr, nf, nc, ud, ncot, nco, k, 
  noise, available, frame.radius, ...)
{
 
  # 6.2.RIGHT Computation of the order of imputation of each MD ----------------

  # Initialization of a list to take all the variable returned by the functions
  ParamList <- list()

  # Creation of the temporary SLG matrices for the right-hand
  # side of OD
  for (h in (nc - 1):(nc - nf + 1)) {
    if (max(ORDERSLGRight[, h]) > 0) {
      # Checking if a gap begins
      # somewhere on the current
      # column.
      # If that is not the case, there is no need to
      # perform the entire following procedure and we
      # simply can go to the next column of ORDERSLGRight

      ParamList[c("ORDERSLGRight_temp", "totV_temp", 
        "nf_temp")] <- SLGMatrixRight_temp(nr, nc, np, h, ORDERSLGRight, 
          nco, ncot, pastDistrib, futureDistrib, k)

      if (max(ParamList$ORDERSLGRight_temp) == 0) {
        next
      }

      # In a similar manner to part 2.4., we go here one
      # single time through the reduced version
      # ORDERSLGRight_temp of ORDERSLG and we create
      # MaxGapSLGLRight "REFORDSLGRight_" matrices
      # collecting the coordinates of
      # each corresponding values in
      # ORDERSLGRight_temp which are
      # greater than 0


      # REFORDSLGRight matrices
      # Initialization of the REFORDSLGRight matrices


      ParamList[c("MaxGap", "REFORD_L", 
        "ORDERSLGRight_temp")] <- REFORDInit(ParamList$ORDERSLGRight_temp,
          nr, nc)


      # 6.3.RIGHT Imputation of the missing data

      for (order in 1:ParamList$MaxGap) {
        ncol_imp <- length(unique(ParamList$REFORD_L[[order]][, 2]))
        col_to_imp <- unique(sort(unique(ParamList$REFORD_L[[order]])[, 2]))
        for (i in 1:ncol_imp) {
          CD_shift <- CDComputeTiming(CO, OD, COt, ParamList$MaxGap, order, np, 
            nc, nr, ParamList$nf_temp, COtsample, pastDistrib, futureDistrib, 
            ncot, k, col_to_imp[i], frame.radius)

          if (length(table(CD_shift$CD[, 1])) > 1) {
            log_CD <- list()
            log_CD[c("reglog", "CD")] <- ComputeModel(CD_shift$CD, regr, 
              ParamList$totV_temp, np, ParamList$nf_temp, k, ...)

            row_to_imp <- ParamList$REFORD_L[[order]][
              which(ParamList$REFORD_L[[order]][, 2] == col_to_imp[i]), 1]
            ODi <- CreatedModelImputationTiming(order, CO, log_CD$CD, COt, OD, 
              ODi, pastDistrib, futureDistrib, available, col_to_imp[i], 
              row_to_imp, ncot, nc, np, ParamList$nf_temp, k, 
              ParamList$totV_temp, regr, log_CD$reglog, noise, 
              CD_shift$shift, ParamList$MaxGap)
          } else {
            lev <- names(table(CD_shift$CD[, 1]))
            REFORDI <- as.matrix(ParamList$REFORD_L[[order]])
            if (ncol(REFORDI) == 1) {
              REFORDI <- t(REFORDI)
            }
            nr_REFORDI <- nrow(REFORDI)

            for (u in 1:nr_REFORDI) {
              i <- REFORDI[u, 1]
              # taking out the first coordinate (row
              # number in ORDER) from REFORDI
              j <- REFORDI[u, 2]
              ODi[i, j] <- lev
            }
          }
        }
      }
    }
  }
  return(ODi)
}
