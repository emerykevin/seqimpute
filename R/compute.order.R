compute.order <- function(data, nr, nc, np, nf, npt, nfi, end.impute) {
  ORDER <- matrix(0, nr, nc)
  ORDER[is.na(data)] <- 1

  imporder <- list()


  imporder$gapInitial <- gap.size(ORDER, nr, 1, 2:nc)
  imporder$maxInitial <- max(imporder$gapInitial)

  imporder$gapTerminal <- gap.size(ORDER, nr, nc, (nc - 1):1)
  imporder$maxTerminal <- max(imporder$gapTerminal)

  for (i in 1:nr) {
    if (imporder$gapInitial[i] != 0) {
      ORDER[i, 1:imporder$gapInitial[i]] <- 0
    }

    if (imporder$gapTerminal[i] != 0) {
      ORDER[i, (nc - imporder$gapTerminal[i] + 1):nc] <- 0
    }
  }


  if (nfi == 0) {
    imporder$maxInitial <- 0
  }

  if (npt == 0 | end.impute == FALSE) {
    imporder$maxTerminal <- 0
  }

  if (imporder$maxInitial != 0) {
    imporder$initial <- order.initial(
      nr, nc, imporder$gapInitial,
      imporder$maxInitial
    )
  } else {
    imporder$initial <- list()
  }

  if (imporder$maxTerminal != 0) {
    imporder$terminal <- order.terminal(
      nr, nc, imporder$gapTerminal,
      imporder$maxTerminal
    )
  } else {
    imporder$terminal <- list()
  }




  ORDER2 <- ORDER
  ORDER3 <- ORDER

  for (i in 1:nr) {
    for (j in 2:nc) {
      if (ORDER[i, j - 1] == 1 & ORDER[i, j] == 1) {
        ORDER2[i, j] <- ORDER2[i, j - 1] + 1
        ORDER3[i, (j - (ORDER2[i, j] - 1)):j] <- ORDER2[i, j]
      }
    }
  }
  maxInternal <- max(max(ORDER2))

  if (max(ORDER) != 0) {
    if (np > 0 & nf > 0) {
      ORDER <- order.pastfuture(ORDER, ORDER3, np, nf, nr, nc, maxInternal)
    } else if (np > 0 & nf == 0) {
      ORDER <- order.past(ORDER, ORDER3, np, nf, nr, nc, maxInternal)
    } else {
      ORDER <- order.future(ORDER, ORDER3, np, nf, nr, nc, maxInternal)
    }

    tmp <- order.SLG(ORDER, nr, nc, np, nf)

    matOrderSLGleft <- tmp$ORDERSLGLeft
    matOrderSLGright <- tmp$ORDERSLGRight
    matOrderSLGboth <- tmp$ORDERSLGBoth
    LongGap <- tmp$LongGap

    ORDER <- ORDER - matOrderSLGleft - matOrderSLGright - matOrderSLGboth
    if (max(ORDER) != 0) {
      imporder[c("maxInternal", "internal")] <- order.internal(ORDER, nr, nc)
    } else {
      imporder[c("maxInternal", "internal")] <- list(0, list())
    }
  } else {
    matOrderSLGleft <- matrix(nrow = nr, ncol = nc, 0)
    matOrderSLGright <- matrix(nrow = nr, ncol = nc, 0)
    matOrderSLGboth <- matrix(nrow = nr, ncol = nc, 0)
    longGap <- FALSE
    imporder$maxInternal <- 0
  }

  if (max(matOrderSLGleft) > 0) {
    imporder$SLGleft <- vector("list", np)
    imporder$maxLeftSLG <- rep(0, np)
    for (h in 2:np) {
      if (max(matOrderSLGleft[, h]) > 0) {
        ORDERSLG_temp <- order.SLG.temp.left(
          nr, nc, h,
          matOrderSLGleft
        )$ORDERSLG_temp
        if (max(ORDERSLG_temp) == 0) {
          next
        }

        tmp <- order.internal(ORDERSLG_temp, nr, nc)

        imporder$SLGleft[[h]] <- tmp$REFORD_L
        imporder$maxLeftSLG[h] <- tmp$maxInternal
      }
    }
  } else {
    imporder$maxLeftSLG <- 0
  }

  if (max(matOrderSLGright) > 0) {
    imporder$SLGright <- vector("list", nc)
    imporder$maxRightSLG <- rep(0, nc - 1)
    for (h in (nc - 1):(nc - nf + 1)) {
      if (max(matOrderSLGright[, h]) > 0) {
        ORDERSLG_temp <- order.SLG.temp.right(
          nr,
          nc, h, matOrderSLGright
        )$ORDERSLGRight_temp
        if (max(ORDERSLG_temp) == 0) {
          next
        }
        tmp <- order.internal(ORDERSLG_temp, nr, nc)
        imporder$SLGright[[h]] <- tmp$REFORD_L
        imporder$maxRightSLG[h] <- tmp$maxInternal
      }
    }
  } else {
    imporder$maxRightSLG <- 0
  }


  if (max(matOrderSLGboth) > 0) {
    imporder$SLGboth <- list()
    imporder$maxBothSLG <- matrix(0, np, nc - 1)
    for (g in 2:np) {
      imporder$SLGboth[[g]] <- list()
      if (sum(matOrderSLGboth[, g - 1] == 0 &
        matOrderSLGboth[, g] != 0) > 0) {
        tt <- which(matOrderSLGboth[, g - 1] == 0 &
          matOrderSLGboth[, g] != 0)
        tmpORDER <- matrix(
          0, nrow(matOrderSLGboth),
          ncol(matOrderSLGboth)
        )
        tmpORDER[tt, g:ncol(matOrderSLGboth)] <- matOrderSLGboth[
          tt,
          g:ncol(matOrderSLGboth)
        ]

        for (h in (nc - 1):(nc - nf + 1)) {
          if (max(tmpORDER[, h]) > 0) {
            ORDERSLG_temp <- order.SLG.temp.right(
              nr, nc, h,
              tmpORDER
            )$ORDERSLGRight_temp
            if (max(ORDERSLG_temp) == 0) {
              next
            }
            tmp <- order.internal(ORDERSLG_temp, nr, nc)
            imporder$SLGboth[[g]][[h]] <- tmp$REFORD_L
            imporder$maxBothSLG[g, h] <- tmp$maxInternal
          }
        }
      }
    }
  } else {
    imporder$maxBothSLG <- 0
  }

  return(imporder)
}



gap.size <- function(ORDER, nr, OrderWidth, OrderList) {
  GapSize <- vector()
  for (i in 1:nr) {
    if (ORDER[i, OrderWidth] == 0) {
      GapSize[i] <- 0
    } else {
      GapSize[i] <- 1
      for (j in OrderList) {
        if (ORDER[i, j] == 1) {
          GapSize[i] <- GapSize[i] + 1
        } else {
          break
        }
      }
    }
  }
  return(GapSize)
}

order.initial <- function(nr, nc, gapInitial, maxInitial) {
  ORDERI <- matrix(0, nr, nc)
  for (i in 1:nr) {
    if (gapInitial[i] != 0) {
      ORDERI[i, 1:gapInitial[i]] <- c(maxInitial:
      (maxInitial + 1 - gapInitial[i]))
    } else {
      next
    }
  }
  initial <- compute.order.edges(maxInitial, nr, ORDERI, maxInitial:1)

  return(initial)
}


order.terminal <- function(nr, nc, gapTerminal, maxTerminal) {
  ORDERT <- matrix(0, nr, nc)
  for (i in 1:nr) {
    if (gapTerminal[i] != 0) {
      ORDERT[i, (nc - gapTerminal[i] + 1):nc] <- c((maxTerminal + 1 - gapTerminal[i]):
      maxTerminal)
    } else {
      next
    }
  }

  terminal <- compute.order.edges(maxTerminal, nr, ORDERT, (nc - maxTerminal + 1):nc)
  return(terminal)
}


order.internal <- function(ORDER, nr, nc) {
  maxInternal <- max(ORDER[ORDER != 0]) - (min(ORDER[ORDER != 0]) - 1)

  REFORD_L <- lapply(1:maxInternal, matrix, data = NA, nrow = 0, ncol = 2)

  non_zero <- which(ORDER > 0, arr.ind = TRUE)
  non_zero <- non_zero[(non_zero[, 1] <= nr) &
    (non_zero[, 2] <= nc), , drop = FALSE]

  ORDER[non_zero] <- ORDER[non_zero] - (min(ORDER[ORDER != 0]) - 1)

  non_zero <- non_zero[order(non_zero[, 1]), , drop = FALSE]
  ord_cord <- ORDER[non_zero]

  for (i in 1:maxInternal) {
    REFORD_L[[i]] <- non_zero[which(ord_cord == i), , drop = FALSE]
  }
  return(list(maxInternal = maxInternal, REFORD_L = REFORD_L))
}


compute.order.edges <- function(GapSize, nr, ORDER, GapSizelist) {
  REFORD_L <- lapply(1:GapSize, matrix, data = NA, nrow = 0, ncol = 2)

  non_zero <- which(ORDER > 0, arr.ind = TRUE)
  non_zero <- non_zero[(non_zero[, 1] <= nr) &
    (non_zero[, 2] %in% GapSizelist), , drop = FALSE]

  ord_cord <- ORDER[non_zero, drop = FALSE]

  for (i in 1:GapSize) {
    REFORD_L[[i]] <- non_zero[which(ord_cord == i), ]
  }

  return(REFORD_L)
}


order.pastfuture <- function(ORDER, ORDER3, np, nf, nr, nc, maxInternal) {
  ord <- integer(maxInternal)
  ord[1] <- 1
  iter_even <- 0
  iter_uneven <- 0
  for (i in 2:maxInternal) {
    if (i %% 2 == 0) {
      shift <- maxInternal - 2 - 3 * iter_even
      iter_even <- iter_even + 1
    } else {
      shift <- -1 - iter_uneven
      iter_uneven <- iter_uneven + 1
    }
    index <- i + shift
    ord[index] <- i
  }

  ifelse(maxInternal %% 2 == 0, ord <- ord, ord <- rev(ord))

  for (i in 1:nr) {
    j <- 1
    while (j <= nc) {
      if (ORDER3[i, j] != 0) {
        if (ORDER3[i, j] %% 2 == 0) {
          ORDER[i, j:(j + ORDER3[i, j] - 1)] <- ord[
            (floor(maxInternal / 2) - ORDER3[i, j] / 2 + 1):
            (floor(maxInternal / 2) + ORDER3[i, j] / 2)
          ]
        } else {
          ORDER[i, j:(j + ORDER3[i, j] - 1)] <- ord[
            (floor(maxInternal / 2) - floor(ORDER3[i, j] / 2) + 1):
            (floor(maxInternal / 2) + ceiling(ORDER3[i, j] / 2))
          ]
        }
        j <- j + ORDER3[i, j] + 1
      } else {
        j <- j + 1
      }
    }
  }

  return(ORDER)
}

order.past <- function(ORDER, ORDER3, np, nf, nr, nc, maxInternal) {
  for (i in 1:nr) {
    j <- 1
    while (j <= nc) {
      if (ORDER3[i, j] > 0) {
        numb <- ORDER3[i, j]
        ord <- c((maxInternal - numb + 1):maxInternal)
        ORDER[i, j:(j + numb - 1)] <- ord
        j <- j + numb + 1
      } else {
        j <- j + 1
      }
    }
  }

  return(ORDER)
}


order.future <- function(ORDER, ORDER3, np, nf, nr, nc, maxInternal) {
  for (i in 1:nr) {
    j <- nc
    while (j >= 1) {
      if (ORDER3[i, j] > 0) {
        numb <- ORDER3[i, j]
        ord <- c(maxInternal:(maxInternal - numb + 1))
        ORDER[i, (j - numb + 1):j] <- ord
        j <- j - numb - 1
      } else {
        j <- j - 1
      }
    }
  }
  return(ORDER)
}

order.SLG <- function(ORDER, nr, nc, np, nf) {
  ORDERSLG <- matrix(0, nrow = nr, ncol = nc)

  tempMinGapLeft <- matrix(0, nrow = nr, ncol = nc)
  tempmaxInternalLeft <- matrix(0, nrow = nr, ncol = nc)
  tempMinGapRight <- matrix(0, nrow = nr, ncol = nc)
  tempmaxInternalRight <- matrix(0, nrow = nr, ncol = nc)

  for (i in 1:nr) {
    if (np > 1) {
      j <- 2

      while (j <= np) {
        jump <- 1

        if (ORDER[i, j] > 0) {
          tempMinGapLeft[i, j] <- j

          while (ORDER[i, j] > 0) {
            ORDERSLG[i, j] <- ORDER[i, j]
            j <- j + 1
          }

          tempmaxInternalLeft[i, j] <- j - 1

          jump <- 1
        }

        j <- j + jump
      }
    }

    if (nf > 1) {
      j <- nc - 1

      while ((nc - j + 1) <= nf) {
        jump <- 1

        if (ORDER[i, j] > 0) {
          tempMinGapRight[i, j] <- j

          while (ORDER[i, j] > 0) {
            ORDERSLG[i, j] <- ORDER[i, j]
            j <- j - 1
          }

          tempmaxInternalRight[i, j] <- j + 1

          jump <- 1
        }

        j <- j - jump
      }
    }
  }

  ORDERSLGLeft <- matrix(nrow = nr, ncol = nc, 0)
  ORDERSLGRight <- matrix(nrow = nr, ncol = nc, 0)
  for (i in 1:nr) {
    if (sum(tempMinGapLeft[i, ]) > 0) {
      minGapLeft <- min(tempMinGapLeft[i, tempMinGapLeft[i, ] != 0])
      maxGapLeft <- max(tempmaxInternalLeft[i, ])
      ORDERSLGLeft[i, minGapLeft:maxGapLeft] <- ORDERSLG[i, minGapLeft:maxGapLeft]
    }

    if (sum(tempMinGapRight[i, ]) > 0) {
      minGapRight <- max(tempMinGapRight[i, tempMinGapRight[i, ] != 0])
      maxGapRight <- min(tempmaxInternalRight[
        i,
        tempmaxInternalRight[i, ] != 0
      ])
      ORDERSLGRight[i, maxGapRight:minGapRight] <- ORDERSLG[
        i,
        maxGapRight:minGapRight
      ]
    }
  }

  LongGap <- FALSE

  ORDERSLGBoth <- matrix(nrow = nr, ncol = nc, 0)

  if (sum(ORDERSLGLeft != 0 & ORDERSLGRight != 0) > 0) {
    LongGap <- TRUE
    ORDERSLGBoth[ORDERSLGLeft != 0 & ORDERSLGRight != 0] <- ORDERSLGLeft[
      ORDERSLGLeft != 0 & ORDERSLGRight != 0
    ]
    ORDERSLGRight[ORDERSLGBoth != 0] <- 0
    ORDERSLGLeft[ORDERSLGBoth != 0] <- 0
  }

  return(list(
    ORDERSLGLeft = ORDERSLGLeft, ORDERSLGRight = ORDERSLGRight,
    ORDERSLGBoth = ORDERSLGBoth, LongGap = LongGap
  ))
}




order.SLG.temp.left <- function(nr, nc, h, ORDERSLG) {
  ORDERSLG_temp <- matrix(0, nrow = nr, ncol = nc)

  for (i in 1:nr) {
    j <- h
    if (ORDERSLG[i, j] > 0 & ORDERSLG[i, j - 1] == 0) {
      while (ORDERSLG[i, j] > 0) {
        ORDERSLG_temp[i, j] <- ORDERSLG[i, j]
        j <- j + 1
      }
    }
  }

  np_temp <- h - 1

  return(list(ORDERSLG_temp = ORDERSLG_temp, np_temp = np_temp))
}


order.SLG.temp.right <- function(nr, nc, h, ORDERSLGRight) {
  ORDERSLGRight_temp <- matrix(0, nrow = nr, ncol = nc)

  for (i in 1:nr) {
    j <- h

    if (ORDERSLGRight[i, j] > 0 & ORDERSLGRight[i, j + 1] == 0) {
      while (ORDERSLGRight[i, j] > 0) {
        ORDERSLGRight_temp[i, j] <- ORDERSLGRight[i, j]
        j <- j - 1
      }
    }
  }

  nf_temp <- nc - h

  return(list(ORDERSLGRight_temp = ORDERSLGRight_temp, nf_temp = nf_temp))
}
