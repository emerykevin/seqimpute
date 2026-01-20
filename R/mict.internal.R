mict.internal <- function(
    data, imp, MaxGap,
    regr, nc, np, nf, nr, ncot, pastDistrib, futureDistrib, k,
    available, REFORD_L, noise, verbose, ...) {
  OD <- data$OD
  covariates <- data$CO
  time.covariates <- data$COt
  COtsample <- data$COtsample

  for (order in 1:MaxGap) {
    if (verbose == TRUE) {
      cat("    Step ", order, "/", MaxGap,"\n")
    }
    shift <- compute.shift(order, MaxGap, timing=FALSE, np, nf)

    train <- compute.traindata(data, MaxGap, order, shift, np, nc, nr, nf, k,
      pastDistrib, futureDistrib,
      col = 0,
      frame.radius = 0, regr, timing = FALSE
    )

    if (length(table(train[,1])) > 1) {
      reglog <- fitmodel(train, regr, ...)

      imp <- impute(
        order, covariates, train,
        time.covariates, COtsample, OD, imp, pastDistrib, futureDistrib,
        available, REFORD_L, ncot, nc, np, nf, k, regr, reglog, noise,
        shift, MaxGap
      )
    } else {
      lev <- names(table(train[,1]))
      REFORD <- as.matrix(REFORD_L[[order]])
      if (ncol(REFORD) == 1) {
        REFORD <- t(REFORD)
      }
      nr_REFORD <- nrow(REFORD)

      for (u in 1:nr_REFORD) {
        i <- REFORD[u, 1]
        j <- REFORD[u, 2]
        imp[i, j] <- lev
      }
    }
  }
  return(imp)
}


impute <- function(order, CO, train, COt, COtsample, OD, ODi,
                   pastDistrib, futureDistrib, available,
                   REFORD_L, ncot, nc, np, nf, k, regr, reglog,
                   noise, shift, MaxGap) {
  if (available == TRUE) {
    LOOKUP <- ODi
  } else {
    LOOKUP <- OD
  }

  REFORD <- as.matrix(REFORD_L[[order]])
  if (ncol(REFORD) == 1) {
    REFORD <- t(REFORD)
  }
  nr_REFORD <- nrow(REFORD)

  CDi <- matrix(NA, nrow = nr_REFORD, ncol = 1)

  if (np > 0) {
    CDpi <- matrix(NA, nrow = nr_REFORD, ncol = np)
  }
  if (nf > 0) {
    CDfi <- matrix(NA, nrow = nr_REFORD, ncol = nf)
  }

  if (pastDistrib) {
    CDdb <- matrix(NA, nrow = nr_REFORD, ncol = k)
  }

  if (futureDistrib) {
    CDda <- matrix(NA, nrow = nr_REFORD, ncol = k)
  }

  if (all(is.na(CO)) == FALSE) {
    if (is.null(dim(CO))) {
      COi <- as.data.frame(matrix(NA, nrow = nr_REFORD, ncol = 1))
    } else {
      COi <- as.data.frame(matrix(NA, nrow = nr_REFORD, ncol = ncol(CO)))
    }
  }
  if (ncot > 0) {
    COtselected <- do.call(rbind, replicate(nr_REFORD, COtsample,
      simplify = FALSE
    ))
  }

  for (u in 1:nr_REFORD) {
    i <- REFORD[u, 1]
    j <- REFORD[u, 2]

    shift <- as.numeric(shift)

    if (np > 0 & nf > 0) {
      CDpi[u, ] <- LOOKUP[i, (j - shift - np):(j - shift - 1)]

      CDfi[u, ] <- LOOKUP[i, (j - shift + MaxGap - order + 1):
      (j - shift + MaxGap - order + nf)]
    } else if (nf == 0) {
      CDpi[u, ] <- LOOKUP[i, (j - np):(j - 1)]
    } else {
      CDfi[u, ] <- LOOKUP[i, (j + 1):(j + nf)]
    }

    if (pastDistrib) {
      CDdb[u, ] <- compute.distrib(LOOKUP[i, , drop = FALSE],
        nc, k, j,
        type = "past"
      )
    }

    if (futureDistrib) {
      CDda[u, ] <- compute.distrib(LOOKUP[i, , drop = FALSE],
        nc, k, j,
        type = "future"
      )
    }

    if (all(is.na(CO)) == FALSE) {
      if (u == 1) {
        if (is.null(dim(CO))) {
          COi <- CO[i, ]
        } else {
          COi <- CO[i, , drop = FALSE]
        }
      } else {
        if (is.null(dim(CO))) {
          COi <- c(COi, CO[i, ])
        } else {
          COi <- rbind(COi, CO[i, , drop = FALSE])
        }
      }
    }
    if (ncot > 0) {
      COttemp <- as.data.frame(matrix(nrow = 1, ncol = 0))
      for (d in 1:(ncot / nc)) {
        COttemp <- cbind(COttemp, COt[i, j + (d - 1) * nc], row.names = NULL)
      }
      COtselected[u, ] <- COttemp
    }
  }

  if (np > 0 & nf > 0) {
    CDi <- cbind(CDi, CDpi, CDfi)
  } else if (np == 0) {
    CDi <- cbind(CDi, CDfi)
  } else {
    CDi <- cbind(CDi, CDpi)
  }
  if (pastDistrib) {
    CDi <- cbind(CDi, CDdb)
  }

  if (futureDistrib) {
    CDi <- cbind(CDi, CDda)
  }

  CDi <- as.data.frame(CDi)

  if (all(is.na(CO)) == FALSE) {
    CDi <- cbind(CDi, COi)
  }

  if (ncot > 0) {
    CDi <- cbind(CDi, COtselected, row.names = NULL)
  }

  colnames(CDi) <- paste("V", 1:ncol(CDi), sep = "")
  if (regr == "rf") {
    for (v in 2:(1 + np + nf)) {
      CDi[, v] <- factor(CDi[, v], levels = c(1:(k + 1)))
      CDi[, v][is.na(CDi[, v])] <- k + 1
    }
    CDi[, 1] <- factor(CDi[, 1], levels = levels(train[, 1]))
  } else {
    CDi[, 1] <- factor(CDi[, 1], levels = c(1:k))
    for (v in 2:(1 + np + nf)) {
      CDi[, v] <- factor(CDi[, v], levels = levels(train[, v]), exclude = NULL)
    }
  }

  if (regr == "multinom") {
    pred <- predict(reglog, CDi, type = "probs")
    names_saved <- reglog$lev
    if (nr_REFORD == 1) {
      pred <- matrix(pred, nrow = 1)
    }

    if (length(reglog$lev) > 2) {
      for (u in 1:nr_REFORD) {
        i <- REFORD[u, 1]
        j <- REFORD[u, 2]

        alea <- runif(1)
        post <- cumsum(pred[u, ])

        sel <- as.numeric(names_saved[which(post >= alea)])

        ODi[i, j] <- sel[1]
      }
    } else {
      for (u in 1:nr_REFORD) {
        i <- REFORD[u, 1]
        j <- REFORD[u, 2]

        alea <- runif(1)
        if (alea > pred[u]) {
          sel <- as.numeric(reglog$lev[1])
        } else {
          sel <- as.numeric(reglog$lev[2])
        }
        ODi[i, j] <- sel
      }
    }
  } else {
    for (u in 1:nr_REFORD) {
      i <- REFORD[u, 1]
      j <- REFORD[u, 2]

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

      ODi[i, j] <- sel[1]
    }
  }

  return(ODi)
}
