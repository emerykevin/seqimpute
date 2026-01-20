mict.timing <- function(
    data, imporder, regr = "multinom", np = 1, nf = 0,
    nfi = 1, npt = 1, available = TRUE, pastDistrib = FALSE,
    futureDistrib = FALSE, m = 1, frame.radius = 0, verbose = TRUE, ...) {
  imp <- data$ODi
  if (imporder$maxInternal > 0) {
    if (verbose) {
      cat("  Imputation of the internal gaps...\n")
    }
    
    imp <- mict.timing.impute(
      data, imp,
      imporder$maxInternal, regr, np, nf, pastDistrib,
      futureDistrib, available,
      imporder$internal, frame.radius, verbose, ...
    )
  }
  if (imporder$maxInitial > 0) {
    if (verbose) {
      cat("  Imputation of the initial gaps...\n")
    }
    imp <- mict.timing.impute(data, imp,
                              imporder$maxInitial, regr,
                              np = 0, nfi,
                              pastDistrib = FALSE, futureDistrib, available,
                              imporder$initial, frame.radius, verbose, ...
    )
  }
  if (imporder$maxTerminal > 0) {
    if (verbose == TRUE) {
      cat("  Imputation of the terminal gaps...\n")
    }
    
    imp <- mict.timing.impute(data, imp,
                              imporder$maxTerminal, regr, npt,
                              nf = 0,
                              pastDistrib, futureDistrib = FALSE, available,
                              imporder$terminal, frame.radius, verbose, ...
    )
  }
  if (max(imporder$maxLeftSLG) > 0) {
    if (verbose == TRUE) {
      cat("  Imputation of the left-hand side SLG...\n")
    }
    
    for (h in 2:np) {
      if (imporder$maxLeftSLG[h] > 0) {
        np_temp <- h - 1
        
        imp <- mict.timing.impute(
          data, imp,
          imporder$maxLeftSLG[h], regr, np_temp,
          nf, pastDistrib, futureDistrib, available,
          imporder$SLGleft[[h]], frame.radius,
          verbose, ...
        )
      }
    }
  }
  if (max(imporder$maxRightSLG) > 0) {
    if (verbose == TRUE) {
      cat("  Imputation of the right-hand side SLG...\n")
    }
    
    for (h in (data$nc - 1):(data$nc - nf + 1)) {
      if (imporder$maxRightSLG[h] > 0) {
        nf_temp <- data$nc - h
        
        imp <- mict.timing.impute(
          data, imp,
          imporder$maxRightSLG[h], regr, np,
          nf_temp, pastDistrib, futureDistrib, available,
          imporder$SLGright[[h]], frame.radius,
          verbose, ...
        )
      }
    }
  }
  if (max(imporder$maxBothSLG) > 0) {
    if (verbose == TRUE) {
      cat("  Imputation of the both-hand side SLG...\n")
    }
    for (g in 2:np) {
      for (h in (data$nc - 1):(data$nc - nf + 1)) {
        if (imporder$maxBothSLG[g, h] > 0) {
          nf_temp <- data$nc - h
          np_temp <- g - 1
          
          imp <- mict.timing.impute(
            data, imp, imporder$maxBothSLG[g, h],
            regr, np_temp, nf_temp, pastDistrib, futureDistrib,
            available, imporder$SLGboth[[g]][[h]],
            frame.radius, verbose, ...
          )
        }
      }
    }
  }
  
  
  return(imp)
}

mict.timing.impute <- function(data, imp, MaxGap, regr, np, nf,
                               pastDistrib, futureDistrib, available,
                               REFORD_L, frame.radius, verbose, ...) {
  nc <- data$nc
  nr <- data$nr
  k <- data$k
  
  for (order in 1:MaxGap) {
    if (verbose == TRUE) {
      cat(paste0("    Step ", order, "/", MaxGap,"\n"))
    }
    
    if (is.vector(REFORD_L[[order]])) {
      col_to_imp <- REFORD_L[[order]][2]
    } else {
      col_to_imp <- unique(sort(unique(REFORD_L[[order]])[, 2]))
    }
    
    shift <- compute.shift(order, MaxGap, timing=TRUE, np, nf)
    
    for (i in 1:length(col_to_imp)) {
      col <- col_to_imp[i]
      
      if (is.vector(REFORD_L[[order]])) {
        row_to_imp <- REFORD_L[[order]][1]
      } else {
        row_to_imp <- REFORD_L[[order]][which(REFORD_L[[order]][, 2] == col), 1]
      }
      train <- compute.traindata(data, MaxGap, order, shift, min(np, col - 1),
                                 nc, nr, min(nf, nc - col), k, pastDistrib,
                                 futureDistrib, col, frame.radius, regr,
                                 timing = TRUE
      )
      
      if (length(table(train[, 1])) > 1) {
        reglog <- fitmodel(train, regr, ...)
        
        imp <- impute.timing(
          data, order, train,
          imp, pastDistrib, futureDistrib, available,
          col, row_to_imp, min(np, col - 1), min(nf, nc - col), regr,
          reglog, shift, MaxGap
        )
      } else {
        lev <- names(table(train[, 1]))
        imp[row_to_imp, col] <- lev
      }
    }
  }
  
  return(imp)
}

compute.traindata <- function(data, MaxGap, order, shift, np, nc, nr, nf, k,
                              pastDistrib, futureDistrib, col, frame.radius, regr,
                              timing = TRUE) {
  OD <- data$OD
  CO <- data$CO
  COt <- data$COt
  nc <- data$nc
  nr <- data$nr
  ncot <- data$ncot
  COtsample <- data$COtsample
  k <- data$k
  
  col_to_use <- compute.coluse(
    timing, MaxGap, np, nf, frame.radius, order,
    shift, nc, col
  )
  
  ud <- length(col_to_use)
  
  iter <- 1
  CD <- matrix(NA, nrow = nr * ud, ncol = 1)
  if (np > 0) {
    CDp <- matrix(NA, nrow = nr * ud, ncol = np)
  }
  if (nf > 0) {
    CDf <- matrix(NA, nrow = nr * ud, ncol = nf)
  }
  
  if (ncot > 0) {
    COtselected <- do.call(rbind, replicate(ud, COtsample, simplify = FALSE))
  }
  
  if (pastDistrib) {
    CDdb <- matrix(NA, nrow = nr * ud, ncol = k)
    db <- matrix(NA, nrow = nr, ncol = k)
  }
  
  if (futureDistrib) {
    CDda <- matrix(NA, nrow = nr * ud, ncol = k)
    da <- matrix(NA, nrow = nr, ncol = k)
  }
  
  for (j in col_to_use) {
    t1 <- (nr * (iter - 1) + 1)
    
    t2 <- nr * iter
    
    CD[t1:t2, 1] <- OD[, j]
    
    if (np > 0 & nf > 0) {
      if (shift == 0) {
        CDp[t1:t2, ] <- OD[, (j - np):(j - 1)]
        CDf[t1:t2, ] <- OD[, (j + MaxGap - order + 1):(j + MaxGap - order + nf)]
      } else {
        CDp[t1:t2, ] <- OD[, (j - shift - np):(j - shift - 1)]
        CDf[t1:t2, ] <- OD[, (j + 1):(j + nf)]
      }
    } else if (np == 0) {
      CDf[t1:t2, ] <- OD[, (j + 1):(j + nf)]
    } else {
      CDp[t1:t2, ] <- OD[, (j - np):(j - 1)]
    }
    
    if (ncot > 0) {
      COtselected <- COtselection(COtselected, COt, ncot, t1, t2, nr, nc, j,
                                  shifted = 0
      )
    }
    
    if (pastDistrib) {
      CDdb[t1:t2, ] <- compute.distrib(OD, nc, k, j, type = "past")
    }
    if (futureDistrib) {
      CDda[t1:t2, ] <- compute.distrib(OD, nc, k, j, type = "future")
    }
    
    iter <- iter + 1
  }
  if (np > 0 & nf > 0) {
    CD <- cbind(CD, CDp, CDf)
  } else if (np == 0) {
    CD <- cbind(CD, CDf)
  } else {
    CD <- cbind(CD, CDp)
  }
  
  if (pastDistrib) {
    CD <- cbind(CD, CDdb)
  }
  
  if (futureDistrib) {
    CD <- cbind(CD, CDda)
  }
  CD <- as.data.frame(CD)
  if (all(is.na(CO)) == FALSE) {
    if (is.null(dim(CO))) {
      CO <- matrix(CO, nrow = nrow(OD), ncol = 1)
    }
    COs <- do.call("rbind", rep(list(CO), ud))
    CD <- cbind(CD, COs)
  }
  if (ncot > 0) {
    CD <- cbind(CD, as.data.frame(COtselected))
  }
  
  npfi <- np + nf
  
  CD <- as.data.frame(CD)
  
  colnames(CD) <- paste("V", 1:ncol(CD), sep = "")
  
  if (regr == "rf") {
    CD[, (1:(1 + npfi))] <- lapply(CD[, (1:(1 + npfi))],
                                   factor,
                                   levels = c(1:k)
    )
    for (v in 2:(1 + npfi)) {
      CD[, v] <- factor(CD[, v], levels = c(1:(k + 1)))
      CD[, v][is.na(CD[, v])] <- k + 1
    }
    CD[, 1] <- factor(CD[, 1], levels = c(1:k))
    CD[, 1] <- droplevels(CD[, 1])
    CD <- CD[!is.na(CD[, 1]), ]
  } else if (regr == "multinom") {
    CD[, 1] <- factor(CD[, 1], levels = c(1:k))
    CD[, 1] <- droplevels(CD[, 1])
    
    if (npfi > 1) {
      CD[, (2:(1 + npfi))] <- lapply(CD[, (2:(1 + npfi))], factor,
                                     levels = c(1:k, NA), exclude = NULL
      )
    } else {
      CD[, 2] <- factor(CD[, 2], levels = c(1:k, NA), exclude = NULL)
    }
  } else {
  }
  
  return(CD)
}


compute.distrib <- function(data, nc, k, col, type) {
  if (type == "future") {
    col.keep <- (col + 1):nc
  } else {
    col.keep <- 1:(col - 1)
  }
  
  datat <- as.data.frame(t(data))
  tempdata <- lapply(datat[col.keep, , drop = FALSE], factor,
                     levels = c(1:k, NA),
                     exclude = NULL
  )
  
  distrib <- (do.call(rbind, lapply(tempdata, summary)) / length(col.keep))[, 1:k]
  
  return(distrib)
}

compute.shift <- function(order, MaxGap, timing, np, nf) {
  if ((np > 0 & nf > 0) & ((MaxGap %% 2 == 0 & order %% 2 == 0) |
                           (MaxGap %% 2 != 0 & order %% 2 != 0))) {
    shift <- MaxGap - order
  } else {
    shift <- 0
  }
  return(shift)
}

compute.coluse <- function(timing, MaxGap, np, nf, frame.radius, order, shift,
                           nc, col) {
  if (timing == TRUE) {
    if ((np > 0 & nf > 0) & ((MaxGap %% 2 == 0 & order %% 2 == 0) |
                             (MaxGap %% 2 != 0 & order %% 2 != 0))) {
      udp <- min(frame.radius, col - (MaxGap - order) - np - 1)
      udf <- min(frame.radius, nc - col - nf)
    } else {
      if (np > 0 & nf > 0) {
        udp <- min(frame.radius, col - np - 1)
        udf <- min(frame.radius, nc - col - (MaxGap - order) - nf)
      }
      if (np > 0 & nf == 0) {
        udp <- min(frame.radius, col - np - 1)
        udf <- min(frame.radius, nc - col)
      }
      if (nf > 0 & np == 0) {
        udf <- min(frame.radius, nc - col - nf)
        udp <- min(col - 1, frame.radius)
      }
    }
    col_to_use <- (col - udp):(col + udf)
  } else {
    if (np > 0 & nf > 0) {
      col_to_use <- (np + 1 + shift):(nc - MaxGap + order + shift - nf)
    } else {
      col_to_use <- (np + 1):(nc - nf)
    }
  }
  return(col_to_use)
}


fitmodel <- function(data, regr, ...) {
  factors_character <- paste("V", 2:ncol(data), sep = "")
  
  factors <- as.vector(factors_character)
  
  fmla <- as.formula(paste("V1~", paste(factors, collapse = "+")))
  
  if (regr == "rf") {
    if ("num.trees" %in% names(list(...))) {
      reglog <- ranger(fmla, data = data, ...)
    } else {
      reglog <- ranger(fmla, data = data, num.trees = 10, ...)
    }
  } else if (regr == "multinom") {
    reglog <- nnet::multinom(data,
                             maxit = 100, trace = FALSE,
                             MaxNwts = 1500, ...
    )
  } else {
  }
  
  return(reglog)
}


COtselection <- function(COtselected, COt, ncot, t1, t2, nr, nc, j, shifted) {
  COttemp <- as.data.frame(matrix(nrow = nr, ncol = 0))
  for (d in 1:(ncot / nc)) {
    COttemp <- cbind(COttemp, COt[, (j + shifted) + (d - 1) * nc])
  }
  COtselected[t1:t2, ] <- COttemp
  
  return(COtselected)
}


impute.timing <- function(data, order, CD, imp,
                          pastDistrib, futureDistrib, available, col,
                          row_to_imp, np, nf, regr, reglog, shift, MaxGap) {
  OD <- data$OD
  CO <- data$CO
  COt <- data$COt
  nc <- data$nc
  nr <- data$nr
  ncot <- data$ncot
  COtsample <- data$COtsample
  k <- data$k
  
  if (available == TRUE) {
    LOOKUP <- imp
  } else {
    LOOKUP <- OD
  }
  
  
  CDi <- matrix(NA, nrow = length(row_to_imp), ncol = 1)
  
  if (np > 0 & nf > 0) {
    if (shift == 0) {
      CDpi <- matrix(LOOKUP[row_to_imp, (col - np):(col - 1)],
                     nrow = length(row_to_imp), ncol = np
      )
      
      CDfi <- matrix(
        LOOKUP[
          row_to_imp,
          (col + MaxGap - order + 1):(col + MaxGap - order + nf)
        ],
        nrow = length(row_to_imp), ncol = nf
      )
    } else {
      CDpi <- matrix(LOOKUP[row_to_imp, (col - shift - np):(col - shift - 1)],
                     nrow = length(row_to_imp), ncol = np
      )
      
      CDfi <- matrix(LOOKUP[row_to_imp, (col + 1):(col + nf)],
                     nrow = length(row_to_imp), ncol = nf
      )
    }
    CDi <- cbind(CDi, CDpi, CDfi)
  } else if (np == 0) {
    CDfi <- matrix(LOOKUP[row_to_imp, (col + 1):(col + nf)],
                   nrow = length(row_to_imp), ncol = nf
    )
    CDi <- cbind(CDi, CDfi)
  } else {
    CDpi <- matrix(LOOKUP[row_to_imp, (col - np):(col - 1)],
                   nrow = length(row_to_imp), ncol = np
    )
    CDi <- cbind(CDi, CDpi)
  }
  
  
  if (pastDistrib) {
    CDdb <- matrix(NA, nrow = length(row_to_imp), ncol = k)
    CDdb[, ] <- compute.distrib(LOOKUP[row_to_imp, , drop = FALSE], nc, k, col,
                                type = "past"
    )
    
    CDi <- cbind(CDi, CDdb)
  }
  
  if (futureDistrib) {
    CDda <- matrix(NA, nrow = length(row_to_imp), ncol = k)
    CDda[, ] <- compute.distrib(LOOKUP[row_to_imp, , drop = FALSE], nc, k, col,
                                type = "future"
    )
    
    CDi <- cbind(CDi, CDda)
  }
  CDi <- as.data.frame(CDi)
  
  
  if (all(is.na(CO)) == FALSE) {
    if (is.null(dim(CO))) {
      COi <- CO[row_to_imp]
    } else {
      COi <- CO[row_to_imp, ]
    }
    CDi <- cbind(CDi, COi)
  }
  
  if (ncot > 0) {
    COtselected_i <- as.data.frame(matrix(nrow = length(row_to_imp), ncol = 0))
    
    indices <- seq(col, by = nc, length.out = ncot / nc)
    COtselected_i <- COt[row_to_imp, indices, drop = FALSE]
    CDi <- cbind(CDi, COtselected_i)
  }
  
  if (regr == "rf") {
    for (v in 2:(1 + np + nf)) {
      CDi[, v] <- factor(CDi[, v], levels = c(1:(k + 1)))
      CDi[, v][is.na(CDi[, v])] <- k + 1
    }
    CDi[, 1] <- factor(CDi[, 1], levels = levels(CD[, 1]))
  } else {
    CDi[, 1] <- factor(CDi[, 1], levels = c(1:k))
    for (v in 2:(1 + np + nf)) {
      CDi[, v] <- factor(CDi[, v], levels = levels(CD[, v]), exclude = NULL)
    }
  }
  
  colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
  
  if (regr == "multinom") {
    pred <- predict(reglog, CDi, type = "probs")
    names_saved <- reglog$lev
    if (length(row_to_imp) == 1) {
      pred <- matrix(pred, nrow = 1)
    }
    
    if (length(reglog$lev) > 2) {
      for (u in 1:length(row_to_imp)) {
        i <- row_to_imp[u]
        
        j <- col
        alea <- runif(1)
        post <- cumsum(pred[u, ])
        
        sel <- as.numeric(names_saved[which(post >= alea)])
        
        imp[i, j] <- sel[1]
      }
    } else {
      for (u in 1:length(row_to_imp)) {
        i <- row_to_imp[u]
        
        j <- col
        alea <- runif(1)
        if (alea > pred[u]) {
          sel <- as.numeric(reglog$lev[1])
        } else {
          sel <- as.numeric(reglog$lev[2])
        }
        imp[i, j] <- sel
      }
    }
  } else {
    for (u in 1:length(row_to_imp)) {
      i <- row_to_imp[u]
      
      j <- col
      
      pred <- predict(reglog,
                      data = CDi[u, ],
                      predict.all = TRUE
      )$predictions[1, ]
      pred <- factor(pred, levels = c(1:k))
      tab <- table(pred)
      tab <- tab / sum(tab)
      
      pred <- cumsum(tab)
      
      
      alea <- runif(1)
      
      sel <- levels(CDi[, 1])[which(pred >= alea)[1]]
      
      imp[i, j] <- sel[1]
    }
  }
  
  return(imp)
}