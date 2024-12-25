
#
#
# ############################################################################
#
# Coding of missing data in function of the length of the gap
# ORDER: matrix of the same size of OD giving the imputation order of each
# MD (0 for observed data and 1 for MD) --> there will be 1 everywhere there
# is MD and 0 everywhere else
# ORDER2: matrix of the same size of OD numbering MD into each gap (for
# example from 1 to 3 for a gap of length 3)
# ORDER3: matrix of the same size of OD replacing each MD by the length of
# the gap it belongs to.

OrderCreation <- function(OD, nr, nc, np, nf, npt, nfi, end.impute) {
  # Creation of matrix ORDER
  ORDER <- matrix(0, nr, nc) # initialization of matrix ORDER with 0 everywhere
  SEL <- is.na(OD) == TRUE # creation of matrix SEL, constituted of TRUE where
  # there is MD in OD and of FALSE everywhere else
  ORDER[SEL] <- 1 # setting some 1 in ORDER at the location where in
  # SEL we have some TRUE


  # Creation of vector InitGapSize (i.e. a vector containing the size of the
  # initial gaps of each line)
  ORDList <- list()
  

  ORDList$InitGapSize <- GapSizeCreation(ORDER, 
                                         nr, 1, 2:nc)
  ORDList$MaxInitGapSize <- max(ORDList$InitGapSize)
  
  # Creation of vector TermGapSize (i.e. a vector containing the size of the
  # terminal gaps of each line)
  
  ORDList$TermGapSize <- GapSizeCreation(ORDER, 
                                         nr, nc, (nc - 1):1)
  ORDList$MaxTermGapSize <- max(ORDList$TermGapSize)
  
    

  # Updating of ORDER with "0" on every external NAs
  # (The purpose of this modification of ORDER is that we don't take into
  # account external NAs at this moment of the program.
  # We will first impute internal gaps and consider external gaps further
  # (as far as nfi and npt are greater than 0))
  for (i in 1:nr) {
    if (ORDList$InitGapSize[i] != 0) {
      ORDER[i, 1:ORDList$InitGapSize[i]] <- vector("numeric", 
        ORDList$InitGapSize[i])
    }

    if (ORDList$TermGapSize[i] != 0) {
      ORDER[i, (nc - ORDList$TermGapSize[i] + 1):nc] <- vector("numeric", 
        ORDList$TermGapSize[i])
    } else {
      next
    }
  }

  if(nfi==0){
    ORDList$MaxInitGapSize <- 0
  }
  
  if(npt==0 | end.impute==FALSE){
    ORDList$MaxTermGapSize <- 0
  }
  # Creation of matrices ORDER2 and ORDER3
  ORDER2 <- ORDER # initially both ORDER2 and
  ORDER3 <- ORDER # ORDER3 are equal to ORDER

  for (i in 1:nr) { # in matrix ORDER, we go line
    # by line...
    for (j in 2:nc) { # ... from column to column
      # (actually exactly as we
      # read a book) (and /!\
      # beginning from column 2!)
      if (ORDER[i, j - 1] == 1 & ORDER[i, j] == 1) { # if the previous value
        # of ORDER is equal to 1
        # and its current value
        # is also equal to 1,
        # then...
        ORDER2[i, j] <- ORDER2[i, j - 1] + 1 # ... construction of
        # ORDER2 for this
        # iteration: the current
        # (j) corresponding
        # value of ORDER2 is
        # assigned by its
        # previous (j-1) value
        # incremented by 1
        # ... construction of ORDER3 for this iteration: all the values
        # in ORDER3 from the beginning of the gap (j-(ORDER2[i,j]-1)) up
        # to the current (j) location are assigned by the current (j)
        # corresponding value in ORDER2
        ORDER3[i, (j - (ORDER2[i, j] - 1)):j] <- ORDER2[i, j]
      }
    }
  }
  MaxGap <- max(max(ORDER2)) 
  
  if (max(ORDER) != 0) {
    ORDER <- PrevAndFutCompute(ORDER, ORDER3, np, nf, nr, nc, MaxGap)
    
    # 2.2. Model 2: use of previous observations only ----------------------------
    ORDER <- PrevObsCompute(ORDER, ORDER3, np, nf, nr, nc, MaxGap)
    
    # 2.3. Model 3: use of future observations only ------------------------------
    ORDER <- FutObsCompute(ORDER, ORDER3, np, nf, nr, nc, MaxGap)
    
    # 6.1 Creation of ORDERSLG (ORDERSLGLeft and ORDERSLGRight)
    Ord_temp_L <- list()
    Ord_temp_L[c("ORDERSLG","tempMinGapLeft","tempMaxGapLeft","tempMinGapRight", 
                 "tempMaxGapRight")] <- ORDERSLGCreation(ORDER, nr, nc, np, nf)
    ORDList[c("ORDERSLGLeft", "ORDERSLGRight", "ORDERSLGBoth", 
              "LongGap")] <- ORDERSLGLRCompute(nr, nc, Ord_temp_L$ORDERSLG, 
                                               Ord_temp_L$tempMinGapLeft, Ord_temp_L$tempMinGapRight, 
                                               Ord_temp_L$tempMaxGapLeft, Ord_temp_L$tempMaxGapRight)
    
    ORDER <- ORDER-ORDList$ORDERSLGLeft-ORDList$ORDERSLGRight-ORDList$ORDERSLGBoth
    if (max(ORDER) != 0) {
      ORDList[c("MaxGap", "REFORD_L", "ORDER")] <- REFORDInit(ORDER, nr, nc)
    } else {
      ORDList[c("MaxGap", "REFORD_L", "ORDER")] <- list(MaxGap, list(), ORDER)
    }
  }else{
    ORDList$ORDERSLGLeft <- matrix(nrow = nr, ncol = nc, 0)
    ORDList$ORDERSLGRight <- matrix(nrow = nr, ncol = nc, 0)
    ORDList$ORDERSLGBoth <- matrix(nrow = nr, ncol = nc, 0)
    ORDList$LongGap <- FALSE
  }
  
  if(ORDList$MaxInitGapSize != 0){
    ORDList$REFORDI_L <- REFORDICreation(nr, nc, ORDList$InitGapSize,
                               ORDList$MaxInitGapSize)
  }else{
    ORDList$REFORDI_L <- list()
  }

  if(ORDList$MaxTermGapSize != 0){
    ORDList$REFORDT_L <- REFORDTCreation(nr, nc, ORDList$TermGapSize,
                                           ORDList$MaxTermGapSize)
  }else{
    ORDList$REFORDT_L <- list()
  }
  return(ORDList)
}


#
# #####Creation of Gapsize
#
# Creation of vector InitGapSize (i.e. a vector containing the size of the
# initial gaps of each line)

GapSizeCreation <- function(ORDER, nr, OrderWidth, OrderList) {
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
  #MaxGapSize <- max(GapSize)
  return(GapSize)
}
