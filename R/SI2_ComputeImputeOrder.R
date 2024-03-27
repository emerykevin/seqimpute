# 2. Computation of the order of imputation of each MD 
#
#
# #############Compute the order of imputation
#
# In case of a factor dataset OD:
# RECODING of OD with numbers "1", "2", etc. instead of its "words"

ImputeOrderComputation <- function(ORDER, ORDER3, MaxGap, np, nf, nr, nc) {
  # 2.1. Model 1: use of previous and future observations ----------------------
  ORDER <- PrevAndFutCompute(ORDER, ORDER3, np, nf, nr, nc, MaxGap)

  # 2.2. Model 2: use of previous observations only ----------------------------
  ORDER <- PrevObsCompute(ORDER, ORDER3, np, nf, nr, nc, MaxGap)

  # 2.3. Model 3: use of future observations only ------------------------------
  ORDER <- FutObsCompute(ORDER, ORDER3, np, nf, nr, nc, MaxGap)

  # 6.1 Creation of ORDERSLG (ORDERSLGLeft and ORDERSLGRight)
  Ord_temp_L <- list()
  Ord_temp_L[c("ORDERSLG","tempMinGapLeft","tempMaxGapLeft","tempMinGapRight", 
    "tempMaxGapRight")] <- ORDERSLGCreation(ORDER, nr, nc, np, nf)
  ORDList <- list()
  ORDList[c("ORDERSLGLeft", "ORDERSLGRight", "ORDERSLGBoth", 
    "LongGap")] <- ORDERSLGLRCompute(nr, nc, Ord_temp_L$ORDERSLG, 
    Ord_temp_L$tempMinGapLeft, Ord_temp_L$tempMinGapRight, 
    Ord_temp_L$tempMaxGapLeft, Ord_temp_L$tempMaxGapRight)


  # Dummy that capture if there are gaps that are both too close to the
  # right edge of the sequence and the left edge of the sequence. It will
  # be imputed as a special case later on.


  # /!\ Final version of the matrix ORDER that we use through point 3.1 to
  # 3.3 of the program

  ORDER <- ORDER-ORDList$ORDERSLGLeft-ORDList$ORDERSLGRight-ORDList$ORDERSLGBoth

  # 2.4. Creation of matrices REFORD -------------------------------------------
  if (max(ORDER) != 0) {
    ORDList[c("MaxGap", "REFORD_L", "ORDER")] <- REFORDInit(ORDER, nr, nc)
  } else {
    ORDList[c("MaxGap", "REFORD_L", "ORDER")] <- list(MaxGap, list(), ORDER)
  }

  return(ORDList)
}



################################################################################
# Model 1: use of previous and future observations

PrevAndFutCompute <- function(ORDER, ORDER3, np, nf, nr, nc, MaxGap) {
  if (np > 0 & nf > 0) {
    ord <- integer(MaxGap) # initialization of the vector "ord"

    # Creation of the longest vector of the matrix
    ord[1] <- 1
    iter_even <- 0
    iter_uneven <- 0
    for (i in 2:MaxGap) {
      if (i %% 2 == 0) {
        shift <- MaxGap - 2 - 3 * iter_even
        iter_even <- iter_even + 1
      } else {
        shift <- -1 - iter_uneven
        iter_uneven <- iter_uneven + 1
      }
      index <- i + shift
      ord[index] <- i
    }

    ifelse(MaxGap %% 2 == 0, ord <- ord, ord <- rev(ord)) # reverse the order of
    # ord in case we are
    # in an even-first
    # case (that is in the
    # case of an uneven
    # MaxGap)


    # Creation of every shorter vector based on "ord" (taking some parts
    # of "ord")
    for (i in 1:nr) {
      j <- 1
      while (j <= nc) {
        if (ORDER3[i, j] != 0) { # meeting a value in ORDER3
          if (ORDER3[i, j] %% 2 == 0) {
            ORDER[i, j:(j + ORDER3[i, j] - 1)] <- ord[
              (floor(MaxGap / 2) - ORDER3[i, j] / 2 + 1):
              (floor(MaxGap / 2) + ORDER3[i, j] / 2)]
          } else {
            ORDER[i, j:(j + ORDER3[i, j] - 1)] <- ord[
              (floor(MaxGap / 2) - floor(ORDER3[i, j] / 2) + 1):
              (floor(MaxGap / 2) + ceiling(ORDER3[i, j] / 2))]
          }
          j <- j + ORDER3[i, j] + 1
        } else {
          j <- j + 1
        }
      }
    }
  }
  return(ORDER)
}



################################################################################
# Model 2: use of previous observations only

PrevObsCompute <- function(ORDER, ORDER3, np, nf, nr, nc, MaxGap) {
  if (np > 0 & nf == 0) { # Verifying that we are in case of model 2
    for (i in 1:nr) { # Beginning from row 1, we will go row by
      # row in the matrix
      j <- 1 # Setting the counter of the columns to 1
      while (j <= nc) { # We will cover the columns from
        # left to right (until we reach the
        # end of the row)
        if (ORDER3[i, j] > 0) { # If we meet a component
          # of ORDER3, it means
          # that we are entering a
          # gap of length equal to
          # this component
          numb <- ORDER3[i, j] # We then store this
          # size of gap in "numb"
          ord <- c((MaxGap - numb + 1):MaxGap) # We create a vector
          # "ord" of length numb
          # but going from where
          # we are in the
          # "residual part to
          # fill" to MaxGap
          ORDER[i, j:(j + numb - 1)] <- ord # We then insert the
          # vector "ord" at the
          # location where it has
          # to go on the current
          # row of the final
          # matrix ORDER
          j <- j + numb + 1 # We increment the counter of the
          # columns in order to skip the
          # entire gap ("+numb") and to skip
          # the directly next position as well
          # ("+1") (because we know that the
          # following location just after a
          # gap won't be another gap
          # (inevitably!)), before continuing
          # to analyze the rest of the current
          # row
        } else { # Otherwise...
          j <- j + 1 # ... we just simply increment the
          # counter of the columns to go to
          # the next column during the next
          # iteration
        }
      }
    }
  }
  return(ORDER)
}



################################################################################
# Model 3: use of future observations only

FutObsCompute <- function(ORDER, ORDER3, np, nf, nr, nc, MaxGap) {
  if (np == 0 & nf > 0) { # Verifying that we are effectively
    # in case of model 3
    for (i in 1:nr) { # Beginning from row 1, we will go
      # row by row in the matrix
      j <- nc # Setting the counter of
      # the columns to the end (i.e. the
      # total number of columns "nc")
      while (j >= 1) { # We will cover the columns from
        # right to left (until we reach
        # the beginning of the row)
        if (ORDER3[i, j] > 0) { # If we meet a
          # component of ORDER3,
          # it means that we are
          # entering a gap of
          # length equal to this
          # component
          numb <- ORDER3[i, j] # We then store this
          # size of gap in
          # "numb"
          ord <- c(MaxGap:(MaxGap - numb + 1)) # We create a vector
          # "ord" of length numb
          # but going from
          # MaxGap to where we
          # are in the "residual
          # part to fill"
          ORDER[i, (j - numb + 1):j] <- ord # We then insert the
          # vector "ord" at the
          # location where it
          # has to go (that is
          # "j-numb+1" elements
          # more to the left) on
          # the current row of
          # the final matrix
          # ORDER
          j <- j - numb - 1 # We decrement the counter of the
          # columns in order to skip the entire
          # gap going to the left ("-numb") and to
          # skip the directly next position as
          # well (still going to the left) ("-1")
          # (because we know that the previous
          # location just before a gap won't be
          # another gap (inevitably!)), before
          # continuing to analyze the "rest" of
          # the current row
        } else { # Otherwise...
          j <- j - 1 # ... we just simply decrement the
          # counter of the columns to go to the
          # previous column (i.e. one more column
          # on the left) during the next iteration
        }
      }
    }
  }
  return(ORDER)
}



################################################################################
# Creation of ORDERSLG (ORDERSLGLeft and ORDERSLGRight)
#
# Updating ORDER with "0" on every NAs belonging to a Specially Located Gap
# (SLG) (The purpose of this modification of ORDER is that we don't take into
# account SLG NAs at this moment of the program.
# We will first impute internal gaps, external gaps and consider SLG at the
# very end (as far as some SLG have been detected)

ORDERSLGCreation <- function(ORDER, nr, nc, np, nf) {
  # Initialization of matrix in which we will store the SLG
  ORDERSLG <- matrix(0, nrow = nr, ncol = nc)

  # Initialization of the range in which SLG could be found
  tempMinGapLeft <- matrix(0, nrow = nr, ncol = nc)
  tempMaxGapLeft <- matrix(0, nrow = nr, ncol = nc)
  tempMinGapRight <- matrix(0, nrow = nr, ncol = nc)
  tempMaxGapRight <- matrix(0, nrow = nr, ncol = nc)

  for (i in 1:nr) { # we will go through each line of ORDER

    if (np > 1) { # if np > 1, it may be possible that
      # SLG on the left-hand side of OD exist
      j <- 2

      while (j <= np) {
        jump <- 1

        if (ORDER[i, j] > 0) {
          tempMinGapLeft[i, j] <- j

          while (ORDER[i, j] > 0) {
            ORDERSLG[i, j] <- ORDER[i, j]
            j <- j + 1
          }

          tempMaxGapLeft[i, j] <- j - 1

          # jump <- max(tempMaxGapLeft[i,]) - max(tempMinGapLeft[i,])
          jump <- 1
        }

        j <- j + jump
      }
    }

    if (nf > 1) { # if nf > 1, it may be possible that SLG
      # on the right-hand side of OD exist
      j <- nc - 1

      while ((nc - j + 1) <= nf) {
        jump <- 1

        if (ORDER[i, j] > 0) {
          tempMinGapRight[i, j] <- j

          while (ORDER[i, j] > 0) {
            ORDERSLG[i, j] <- ORDER[i, j]
            j <- j - 1
          }

          tempMaxGapRight[i, j] <- j + 1

          # jump <- max(tempMinGapRight[i,]) - max(tempMaxGapRight[i,])
          jump <- 1
        }

        j <- j - jump
      }
    }
  }
  return(list(ORDERSLG, tempMinGapLeft, tempMaxGapLeft, 
    tempMinGapRight, tempMaxGapRight))
}



################################################################################
# Computation of ORDERSLG (ORDERSLGLeft and ORDERSLGRight)
#
# Extracting extrema from tempMinGapLeft, tempMaxGapLeft,
# tempMinGapRight and tempMaxGapRight
# And creation of ORDERSLGLeft and ORDERSLGRight (matrices for both
# groups of SLG (i.e. one on the left- and the other one on the right-
# hand side of the matrix ORDERSLG)

ORDERSLGLRCompute <- function(nr, nc, ORDERSLG, tempMinGapLeft, 
  tempMinGapRight, tempMaxGapLeft, tempMaxGapRight)
{
  ORDERSLGLeft <- matrix(nrow = nr, ncol = nc, 0)
  ORDERSLGRight <- matrix(nrow = nr, ncol = nc, 0)
  for (i in 1:nr) {
    if (sum(tempMinGapLeft[i, ]) > 0) {
      minGapLeft <- min(tempMinGapLeft[i, tempMinGapLeft[i, ] != 0])
      maxGapLeft <- max(tempMaxGapLeft[i, ])
      ORDERSLGLeft[i,minGapLeft:maxGapLeft] <- ORDERSLG[i,minGapLeft:maxGapLeft]
    }

    if (sum(tempMinGapRight[i, ]) > 0) {
      minGapRight <- max(tempMinGapRight[i, tempMinGapRight[i, ] != 0])
      maxGapRight <- min(tempMaxGapRight[i, tempMaxGapRight[i, ] != 0])
      ORDERSLGRight[i, maxGapRight:minGapRight] <- ORDERSLG[i,
        maxGapRight:minGapRight]
    }
  }

  LongGap <- FALSE

  # We create the matrix ORDERSLGBoth
  ORDERSLGBoth <- matrix(nrow = nr, ncol = nc, 0)

  if (sum(ORDERSLGLeft != 0 & ORDERSLGRight != 0) > 0) {
    LongGap <- TRUE
    ORDERSLGBoth[ORDERSLGLeft != 0 & ORDERSLGRight != 0] <- ORDERSLGLeft[
      ORDERSLGLeft != 0 & ORDERSLGRight != 0]
    ORDERSLGRight[ORDERSLGBoth != 0] <- 0
    ORDERSLGLeft[ORDERSLGBoth != 0] <- 0
  }

  return(list(ORDERSLGLeft, ORDERSLGRight, ORDERSLGBoth, LongGap))
}
