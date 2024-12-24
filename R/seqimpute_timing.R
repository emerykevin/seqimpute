seqimpute_timing <- function(dataOD, var, imporder, regr = "multinom", np = 1, nf = 0, nfi = 1, 
  npt = 1, available = TRUE, covariates = matrix(NA, nrow = 1, ncol = 1),
  time.covariates = matrix(NA, nrow = 1, ncol = 1), pastDistrib = FALSE,
  futureDistrib = FALSE, noise=0, m = 1, SetRNGSeed = FALSE, 
  end.impute = TRUE, ParExec = TRUE, 
  ncores = NULL, frame.radius = 0, verbose = TRUE, ...)
{
  
  dataOD <- preliminaryChecks(dataOD, covariates, time.covariates,var=var)
  
  rownamesDataset <- rownames(dataOD$OD)
  nrowsDataset <- nrow(dataOD$OD)
  noise <- 0
  
  if (sum(is.na(dataOD$OD)) == 0) {
    if (verbose == TRUE) {
      message("This dataset has no missing values!")
    }
    return(dataOD$OD)
  }
  
  tmp <- check.predictors(np, nf, nfi, npt)
  np <- tmp$np
  nf <- tmp$nf
  nfi <- tmp$nfi
  npt <- tmp$npt
  
  regr <- check.regr(regr)
  dataOD$ncot <- check.ncot(dataOD$ncot,dataOD$nc)
  
  imporder <- OrderCreation(dataOD$OD, dataOD$nr, dataOD$nc, np, nf, npt, nfi, end.impute)
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
    
    if (max(imporder$ORDER)> 0) {
      if (verbose) {
        print("Imputation of the internal gaps...")
      }
      
      dataOD[["ODi"]] <- ModelImputationTiming(dataOD, dataOD$ODi,
        imporder$MaxGap, regr, np, nf, pastDistrib, futureDistrib, available,
        imporder$REFORD_L, frame.radius, verbose, ...)
    }
    if(imporder$MaxInitGapSize > 0){
      if(verbose){print("Imputation of the initial gaps...")}
      dataOD[["ODi"]] <- ModelImputationTiming(dataOD, dataOD$ODi,
                            imporder$MaxInitGapSize, regr, np=0, nfi, pastDistrib=FALSE, futureDistrib, available,
                            imporder$REFORDI_L, frame.radius, verbose, ...)
      
    }
    if (imporder$MaxTermGapSize >0) {
      if (verbose == TRUE) {
        print("Imputation of the terminal gaps...")
      }

      dataOD[["ODi"]] <- ModelImputationTiming(dataOD, dataOD$ODi,
                              imporder$MaxTermGapSize, regr, npt, nf=0, pastDistrib, 
                              futureDistrib=FALSE, available,
                              imporder$REFORDT_L, frame.radius, verbose, ...)
    }
    if (max(imporder$ORDERSLGLeft) != 0 & !is.null(imporder$ORDERSLGLeft)) {
      # Checking if we have to impute
      # left-hand side SLG
      if (verbose == TRUE) {
        print("Imputation of the left-hand side SLG...")
      }
      dataOD[["ODi"]] <- LSLGNAsImputeTiming(dataOD$OD, dataOD$ODi, dataOD$CO, 
        dataOD$COt, dataOD$COtsample, imporder$ORDERSLGLeft, pastDistrib, 
        futureDistrib, regr, np, dataOD$nr, nf, dataOD$nc, dataOD$ud, 
        dataOD$ncot, dataOD$nco, dataOD$k, dataOD$noise, available, 
        frame.radius, ...)
    }
    # right-hand side SLG
    if (max(imporder$ORDERSLGRight) != 0 & !is.null(imporder$ORDERSLGRight)) {
      # Checking if we have to impute right-hand
      # side SLG
      if (verbose == TRUE) {
        print("Imputation of the right-hand side SLG...")
      }

      dataOD[["ODi"]] <- RSLGNAsImputeTiming(dataOD$OD, dataOD$ODi, dataOD$CO, 
        dataOD$COt, dataOD$COtsample, imporder$ORDERSLGRight, pastDistrib, 
        futureDistrib, regr, np, dataOD$nr, nf, dataOD$nc, dataOD$ud, 
        dataOD$ncot, dataOD$nco, dataOD$k, dataOD$noise, 
        available, frame.radius, ...)
    }
    # Checking if we have to impute
    # Both-hand side SLG
    if (imporder$LongGap) {
      if (verbose == TRUE) {
        print("Imputation of the both-hand side SLG...")
      }
      for (h in 2:np) {
        if (sum(imporder$ORDERSLGBoth[, h - 1] == 0 & 
          imporder$ORDERSLGBoth[, h] != 0) > 0) {
          tt <- which(imporder$ORDERSLGBoth[, h - 1] == 0 & 
            imporder$ORDERSLGBoth[, h] != 0)
          tmpORDER <- matrix(0, nrow(imporder$ORDERSLGBoth), 
            ncol(imporder$ORDERSLGBoth))
          tmpORDER[tt, h:ncol(imporder$ORDERSLGBoth)] <- imporder$ORDERSLGBoth[tt, 
            h:ncol(imporder$ORDERSLGBoth)]
          dataOD[["ODi"]] <- RSLGNAsImputeTiming(dataOD$OD, dataOD$ODi, 
            dataOD$CO, dataOD$COt, dataOD$COtsample, tmpORDER, 
            pastDistrib, futureDistrib, regr, h - 1, dataOD$nr, 
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
  
  RESULT <- lapply(RESULT,FinalResultConvert, ODClass = dataOD$ODClass,
      ODlevels = dataOD$ODlevels, rownamesDataset = rownamesDataset, 
      nrowsDataset = nrowsDataset, nr = dataOD$nr, nc = dataOD$nc, 
      rowsNA = dataOD$rowsNA, mi = m)
  
  return(RESULT)
  
  
}

ModelImputationTiming <- function(data, ODi, MaxGap,regr, np, nf, 
                        pastDistrib, futureDistrib, available, REFORD_L, 
                        frame.radius, verbose,...){
  
  OD <- data$OD
  covariates <- data$CO
  time.covariates <- data$COt
  nc <- data$nc
  nr <- data$nr
  ncot <- data$ncot
  COtsample <- data$COtsample
  k <- data$k
  noise <- 0
  
  for (order in 1:MaxGap) {
    if(verbose == TRUE){
      print(paste0("Step ", order, "/", MaxGap))
    }

    if(is.vector(REFORD_L[[order]])){
      col_to_imp <- REFORD_L[[order]][2]
    }else{
      col_to_imp <- unique(sort(unique(REFORD_L[[order]])[, 2]))
    }
    ncol_imp <- length(col_to_imp)
    for (i in 1:ncol_imp){
      col <- col_to_imp[i]
      CD_shift <- CDComputeTiming(covariates, OD, time.covariates, MaxGap,
        order, min(np, col - 1), nc, nr, min(nf, nc - col), COtsample, pastDistrib, futureDistrib, ncot, k,
        col, frame.radius)
      if (length(table(CD_shift$CD[, 1])) > 1) {
        log_CD <- list()
        log_CD[c("reglog", "CD")] <- ComputeModel(CD_shift$CD, regr,
          min(np, col - 1), min(nf, nc - col), k, ...)
        if(is.vector(REFORD_L[[order]])){
          row_to_imp <- REFORD_L[[order]][1]
        }else{
          row_to_imp <- REFORD_L[[order]][
            which(REFORD_L[[order]][, 2] == col_to_imp[[i]]), 1]
        }


        ODi <- CreatedModelImputationTiming(order, covariates, log_CD$CD,
          time.covariates, OD, ODi, pastDistrib, futureDistrib, available,
          col, row_to_imp, ncot, nc, min(np, col - 1), min(nf, nc - col), k, regr,
          log_CD$reglog, CD_shift$shift, MaxGap)
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

  if ((np > 0 & nf > 0) & ((MaxGap %% 2 == 0 & order %% 2 == 0) |
                           (MaxGap %% 2 != 0 & order %% 2 != 0))) {
    shift <- MaxGap - order # jumping at the end of
    udp <- min(frame.radius, col - (MaxGap - order) - np - 1) #data in the past
    udf <- min(frame.radius, nc - col - nf) # number of usable data in the futur
    col_to_use <- (col - udp):(col + udf)
    ud <- udp + udf + 1
    # the gap
  } else {
    shift <- 0

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

  iter <- 1
  CD <- matrix(NA, nrow = nr * ud, ncol = 1)
  if(np>0){
    CDp <- matrix(NA, nrow = nr * ud, ncol = np)
  }
  if(nf>0){
    CDf <- matrix(NA, nrow = nr * ud, ncol = nf)
  }

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


    if(np>0&nf>0){
      if (shift == 0) {
        CDp[t1:t2, ] <- OD[, (j - np):(j - 1)]
        CDf[t1:t2, ] <- OD[, (j + MaxGap - order + 1):(j + MaxGap - order + nf)]
      } else {
        CDp[t1:t2, ] <- OD[, (j - shift - np):(j - shift - 1)]
        CDf[t1:t2, ] <- OD[, (j + 1):(j + nf)]
      }
    }else if(np==0){
      CDf[t1:t2,] <- OD[, (j + 1):(j + nf)]
    }else{
      CDp[t1:t2,] <- OD[, (j - np):(j - 1)]
    }

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
  if(np>0&nf>0){
    CD <- cbind(CD, CDp, CDf)
  }else if(np==0){
    CD <- cbind(CD, CDf)
  }else{
    CD <- cbind(CD, CDp)

  }


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

  CD_shift <- list()
  CD_shift[c("CD", "shift")] <- list(CD, shift)
  return(CD_shift)
}

CreatedModelImputationTiming <- function(order, CO, CD, COt, OD, ODi, 
  pastDistrib, futureDistrib, available, col, row_to_imp, ncot, nc, np, nf, k, 
  regr, reglog, shift, MaxGap) {
  
  if (available == TRUE) { 
    LOOKUP <- ODi
  } else { 
    LOOKUP <- OD
  }

  
  CDi <- matrix(NA,nrow=length(row_to_imp),ncol=1)
  
  if(np>0 & nf>0){
    if (shift == 0) {
      CDpi <- matrix(LOOKUP[row_to_imp, (col - np):(col - 1)],nrow=length(row_to_imp),ncol=np)
      CDfi <- matrix(LOOKUP[row_to_imp, (col + MaxGap - order + 1):(col + MaxGap - order + nf)],nrow=length(row_to_imp),ncol=nf)
      
    } else {
      CDpi <- matrix(LOOKUP[row_to_imp, (col - shift - np):(col - shift - 1)],nrow=length(row_to_imp),ncol=np)
      CDfi <- matrix(LOOKUP[row_to_imp, (col + 1):(col + nf)],nrow=length(row_to_imp),ncol=nf)
    }
    CDi <- cbind(CDi, CDpi, CDfi)
  }else if(np==0){
    CDfi <- matrix(LOOKUP[row_to_imp, (col + 1):(col + nf)],nrow=length(row_to_imp),ncol=nf)
    CDi <- cbind(CDi, CDfi)
    
  }else{
    CDpi <- matrix(LOOKUP[row_to_imp, (col - np):(col - 1)],nrow=length(row_to_imp),ncol=np)
    CDi <- cbind(CDi, CDpi)
  }
  
  
  if (pastDistrib) {
    CDdb <- matrix(NA,nrow=length(row_to_imp),ncol=k)
    ODt <- t(LOOKUP[row_to_imp,,drop=FALSE])
    ODt <- as.data.frame(ODt)
    tempOD <- lapply(ODt[(1:(col-1)), ,drop=FALSE], factor, levels = c(1:k, NA),
                     exclude = NULL)
    
    da_list <- lapply(tempOD, summary)
    da_matrix <- do.call(rbind, da_list)
    CDdb[,] <- da_matrix[, 1:k]/length(1:(col-1))
    CDi <- cbind(CDi,CDdb)
  }
  
  # Future distribution (i.e. After)
  if (futureDistrib) {
    CDda <- matrix(NA,nrow=length(row_to_imp),ncol=k)
    ODt <- t(LOOKUP[row_to_imp,,drop=FALSE])
    ODt <- as.data.frame(ODt)
    tempOD <- lapply(ODt[((col + 1):nc), ,drop=FALSE], factor, levels = c(1:k, NA),
                     exclude = NULL)
    
    da_list <- lapply(tempOD, summary)
    da_matrix <- do.call(rbind, da_list)
    CDda[,] <- da_matrix[, 1:k]/length((col + 1):nc)
    CDi <- cbind(CDi,CDda)
  }
  CDi <- as.data.frame(CDi)
  
  
  if (all(is.na(CO)) == FALSE) {
    
    if (is.null(dim(CO))) {
      COi <- CO[row_to_imp]
    } else {
      COi <- CO[row_to_imp,]
    }
    CDi <- cbind(CDi, COi)
    
  }
  
  if (ncot > 0) {
    COtselected_i <- as.data.frame(matrix(nrow = length(row_to_imp), ncol = 0))
    for (d in 1:(ncot / nc)) {
      COtselected_i <- cbind(COtselected_i, COt[row_to_imp, (col) + (d - 1) * nc])
    }
    CDi <- cbind(CDi, COtselected_i)
  }
  
  if (regr == "rf") {
    for (v in 2:(1 + np + nf)) {
      CDi[, v] <- factor(CDi[, v], levels = c(1:(k + 1)))
      CDi[, v][is.na(CDi[, v])] <- k + 1
    }
    CDi[, 1] <- factor(CDi[, 1], levels = levels(CD[, 1]))
  }else {
    # multinom
    CDi[, 1] <- factor(CDi[, 1], levels = c(1:k))
    for (v in 2:(1 + np + nf)) {
      CDi[, v] <- factor(CDi[, v], levels = levels(CD[, v]), exclude = NULL)
    }
  }
  
  colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
  
  if(regr=="multinom"){
    pred <- predict(reglog, CDi, type = "probs")
    names_saved <- reglog$lev
    if(length(reglog$lev)>2){
      for(u in 1:length(row_to_imp)){
        
        i <- row_to_imp[u]
        
        j <- col
        alea <- runif(1)
        if (alea > pred[u]) {
          sel <- as.numeric(reglog$lev[1])
        } else {
          sel <- as.numeric(reglog$lev[2])
        }
        ODi[i, j] <- sel[1]
      }
    }else{
      for(u in 1:length(row_to_imp)){
        i <- row_to_imp[u]
        
        j <- col
        alea <- runif(1)
        if (alea > pred[u]) {
          sel <- as.numeric(reglog$lev[1])
        } else {
          sel <- as.numeric(reglog$lev[2])
        }
        ODi[i, j] <- sel
      }
    }
  }else{
    pred <- predict(reglog, data = CDi, predict.all = TRUE)$predictions
    if(length(row_to_imp)==1){
      pred <- matrix(pred,nrow=1)
    }
    for(u in 1:length(row_to_imp)){
      i <- row_to_imp[u]
      
      j <- col
      post <- factor(pred[u,], levels = c(1:k))
      tab <- table(pred)
      tab <- tab / sum(tab)
      
      post <- cumsum(tab)
      
      alea <- runif(1)
      
      sel <- levels(CDi[, 1])[which(post >= alea)[1]]
      ODi[i, j] <- sel
    }
    
  }
  
  return(ODi)
}



LSLGNAsImputeTiming <- function(OD, ODi, CO, COt, COtsample, ORDERSLG, 
  pastDistrib, futureDistrib, regr, np, nr, nf, nc, ud, ncot, nco, k, noise, 
  available, frame.radius, ...)
  {
  ParamList <- list()
  
  for (h in 2:np) {
    if (max(ORDERSLG[, h]) > 0) {
      # Checking if a gap begins
      # somewhere on the current column
      # If that is not the case, there is
      # no need to perform
      # the entire following procedure
      # and we simply can go
      # to the next column of ORDERSLGLeft

      ParamList[c("ORDERSLG_temp","np_temp")] <- SLGMatrix_temp(nr, 
        nc, nf, h, ORDERSLG, nco, ncot, pastDistrib, futureDistrib, k)

      if (max(ParamList$ORDERSLG_temp) == 0) {
        next
      }

      ParamList[c("MaxGap", "REFORD_L", 
        "ORDERSLG_temp")] <- REFORDInit(ParamList$ORDERSLG_temp, nr, nc)

     
     
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
                ParamList$np_temp, nf, k, ...)

            row_to_imp <- ParamList$REFORD_L[[order]][
              which(ParamList$REFORD_L[[order]][, 2] == col_to_imp[i]), 1]
            # 3.3. Imputation using the just created model 
            ODi <- CreatedModelImputationTiming(order, CO, log_CD$CD, COt, OD, 
              ODi, pastDistrib, futureDistrib, available, col_to_imp[i], 
              row_to_imp, ncot, nc, ParamList$np_temp, nf, k, 
              regr, log_CD$reglog, 
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



RSLGNAsImputeTiming <- function(OD, ODi, CO, COt, COtsample, ORDERSLGRight, 
  pastDistrib, futureDistrib, regr, np, nr, nf, nc, ud, ncot, nco, k, 
  noise, available, frame.radius, ...)
{
 

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

      ParamList[c("ORDERSLGRight_temp",  
        "nf_temp")] <- SLGMatrixRight_temp(nr, nc, np, h, ORDERSLGRight, 
          nco, ncot, pastDistrib, futureDistrib, k)

      if (max(ParamList$ORDERSLGRight_temp) == 0) {
        next
      }

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
              np, ParamList$nf_temp, k, ...)

            row_to_imp <- ParamList$REFORD_L[[order]][
              which(ParamList$REFORD_L[[order]][, 2] == col_to_imp[i]), 1]
            
            ODi <- CreatedModelImputationTiming(order, CO, log_CD$CD, COt, OD, 
              ODi, pastDistrib, futureDistrib, available, col_to_imp[i], 
              row_to_imp, ncot, nc, np, ParamList$nf_temp, k, 
              regr, log_CD$reglog,
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
