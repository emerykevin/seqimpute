mict.rightSLG <- function(data, imp,
                          pastDistrib, futureDistrib, regr, np, nr, nf, nc, ud,
                          ncot, nco, k, noise, available, REFORD_L, MaxGap, ...) {
  OD <- data$OD
  covariates <- data$CO
  time.covariates <- data$COt
  COtsample <- data$COtsample

  for (h in (nc - 1):(nc - nf + 1)) {
    if (MaxGap[h] > 0) {
      nf_temp <- nc - h
      for (order in 1:MaxGap[h]) {
        shift <- compute.shift(order, MaxGap[h], timing, np, nf_temp)

        train <- compute.traindata(data, MaxGap[h], order, shift, np, nc, nr,
          nf_temp, k, pastDistrib, futureDistrib,
          col = 0, frame.radius = 0, regr, timing = FALSE
        )


        if (length(table(train[, 1])) > 1) {
          reglog <- fitmodel(train, regr, ...)

          imp <- impute(
            order, covariates, train,
            time.covariates, COtsample, OD, imp,
            pastDistrib, futureDistrib, available,
            REFORD_L[[h]], ncot, nc, np, nf_temp,
            k, regr, reglog, noise,
            shift, MaxGap[h]
          )
        } else {
          lev <- names(table(train[, 1]))
          REFORD <- as.matrix(REFORD_L[[h]][[order]])
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
    }
  }
  return(imp)
}
