# 6. Imputation of SLG NAs
#
#
################################################################################
# Left-hand side SLG imputation

LSLGNAsImpute <- function(OD, ODi, covariates, time.covariates, COtsample, 
  pastDistrib, futureDistrib, regr, np, nr, nf, nc, ud, ncot, nco, 
  k, noise, available, REFORD_L, MaxGap,...) 
{ 
  # # Creation of the temporary SLG matrices for the left-hand
  # # side of OD
  for (h in 2:np) {
    if (MaxGap[h] > 0) {
      # # Checking if a gap begins
      # # somewhere on the current column
      # # If that is not the case, there is
      # # no need to perform
      # # the entire following procedure
      # # and we simply can go
      # # to the next column of ORDERSLGLeft
      # 
      # ParamList[c("ORDERSLG_temp","np_temp")] <- SLGMatrix_temp(nr, 
      #   nc, h, ORDERSLG)
      # 
      # if (max(ParamList$ORDERSLG_temp) == 0) {
      #   next
      # }
      # 
      # # In a similar manner to part 2.4., we go here one
      # # single time through the reduced version
      # # ORDERSLGLeft_temp of ORDERSLG and we create
      # # MaxGapSLGLeft "REFORDSLGRLeft_" matrices
      # # collecting the coordinates of each corresponding
      # # values in ORDERSLGLeft_temp which are greater
      # # than 0
      # 
      # 
      # # REFORDSLGLeft matrices
      # # Initialization of the REFORDSLGLeft matrices
      # ParamList[c("MaxGap", "REFORD_L", "ORDERSLG_temp")] <- REFORDInit(
      #   ParamList$ORDERSLG_temp, nr, nc)
      # 
      # # MaxGapSLGLeft <- REFORDOD[[1]]
      # # REFORDSLG_L <- REFORDOD[[2]]
      # 
      # # 6.3.1.LEFT Building of the different data matrices CD ----
      # # for the computation of the model for every SLG
      # # on the left-hand side of OD
      # 
      np_temp <- h-1
      for (order in 1:MaxGap[h]) {
        ParamList <- list()
        # ParamList[c("CD", "shift")] <- SLGCDMatBuild(covariates, 
        #   time.covariates, OD, order, MaxGap[h], np_temp, 
        #   ncot, nr, nc, nf, COtsample, pastDistrib, futureDistrib, k)
        
        ParamList[c("CD", "shift")] <- CDCompute(covariates, OD, time.covariates, 
                                      MaxGap[h], order, np_temp, nc, nr, nf, COtsample, 
                                      pastDistrib, futureDistrib, ncot, k) 
        
        # 6.3.2.LEFT Computation of the model 
        if (length(table(ParamList$CD[, 1])) > 1) {
          log_CD <- list()
          log_CD[c("reglog", "CD")] <- ComputeModel(ParamList$CD, regr, 
              np_temp, nf, k, ...)
          
          # 6.3.3.LEFT Imputation using the just created model 
          # ODi <- SLGCreatedModelImpute(covariates, OD, log_CD$CD, ODi, 
          #     time.covariates, ncot, nf, nc, regr, k, 
          #     log_CD$reglog, noise, available, REFORD_L[[h]], 
          #     np_temp, pastDistrib, futureDistrib, 
          #     order, MaxGap[h], ParamList$shift)
          
          ODi <- CreatedModelImputation(order, covariates, log_CD$CD, 
                                        time.covariates, OD, ODi, pastDistrib, futureDistrib, available, 
                                        REFORD_L[[h]], ncot, nc, np_temp, nf, k, regr, log_CD$reglog, noise, 
                                        ParamList$shift, MaxGap[h])
        }else{
          lev <- names(table(ParamList$CD[, 1]))
          REFORD <- as.matrix(REFORD_L[[h]][[order]])
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
  }
  return(ODi)
}


##############################################################################
# Right-hand side SLG imputation

RSLGNAsImpute <- function(OD, ODi, covariates, time.covariates, COtsample, 
  pastDistrib, futureDistrib, regr, np, nr, nf, nc, ud, ncot, 
  nco, k, noise, available, REFORD_L, MaxGap,...)
{
  # Checking if we have to impute right-hand
  # side SLG


  # 6.2.RIGHT Computation of the order of imputation of each MD ----------------

  # Initialization of a list to take all the variable returned by the functions

  # Creation of the temporary SLG matrices for the right-hand
  # side of OD
  for (h in (nc - 1):(nc - nf + 1)) {
    if (MaxGap[h] > 0) {
      # 6.3.RIGHT Imputation of the missing data 
      nf_temp <- nc - h
      for (order in 1:MaxGap[h]) {
        ParamList <- list()
        ParamList[c("CD","shift")] <- CDCompute(covariates, OD, time.covariates, 
                                                MaxGap[h], order, np, nc, nr, nf_temp, COtsample, 
                                                pastDistrib, futureDistrib, ncot, k) 
        if (length(table(ParamList$CD[, 1])) > 1) {
          log_CD <- list()
          log_CD[c("reglog", "CD")] <- ComputeModel(ParamList$CD, regr, 
              np, nf_temp, k, ...)
          
          ODi <- CreatedModelImputation(order, covariates, log_CD$CD, 
                                        time.covariates, OD, ODi, pastDistrib, futureDistrib, available, 
                                        REFORD_L[[h]], ncot, nc, np, nf_temp, k, regr, log_CD$reglog, noise, 
                                        ParamList$shift, MaxGap[h])
        }else{
          lev <- names(table(ParamList$CD[, 1]))
          REFORD <- as.matrix(REFORD_L[[h]][[order]])
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
  }
  return(ODi)
}

SLGMatrix_temp <- function(nr, nc, h, ORDERSLG) 
{
  # Creation of the temporary SLG matrices for the left-hand
  # side of OD

  # initialization of a new temporary
  # ORDERSLG_temp matrix
  ORDERSLG_temp <- matrix(0, nrow = nr, ncol = nc)

  for (i in 1:nr) { # going from top to bottom
    j <- h # storing the current column
    # we are checking (and
    # reinitializing j)

    if (ORDERSLG[i, j] > 0 & ORDERSLG[i, j - 1] == 0) {
      while (ORDERSLG[i, j] > 0) {
        # going from left to right
        ORDERSLG_temp[i, j] <- ORDERSLG[i, j]
        j <- j + 1
      }
    }
  }



  # Update of np and totV
  np_temp <- h - 1


  # (i.e. updating matrix ORDERSLG_temp with
  # coherent value of "order"
  # (going from 1 to MaxGap))

  # Adjusting the matrix ORDERSLG_temp and
  # storing the coordinates of the NA to impute
  return(list(ORDERSLG_temp=ORDERSLG_temp, np_temp=np_temp))
}




SLGMatrixRight_temp <- function(nr, nc, h, ORDERSLGRight) 
{
  # initialization of a new temporary
  # ORDERSLGRight_temp matrix
  ORDERSLGRight_temp <- matrix(0, nrow = nr, ncol = nc)

  for (i in 1:nr) {
    j <- h # storing the current column we are
    # checking (and reinitializing j)

    if (ORDERSLGRight[i, j] > 0 & ORDERSLGRight[i, j + 1] == 0) {
      while (ORDERSLGRight[i, j] > 0) {
        ORDERSLGRight_temp[i, j] <- ORDERSLGRight[i, j]
        j <- j - 1
      }
    }
  }



  # Update of nf and totV
  nf_temp <- nc - h
  
  # (i.e. updating matrix ORDERSLGRight_temp with
  # coherent value of "order" (going from 1 to
  # MaxGapSLGRight))

  # Adjusting the matrix ORDERSLGRight_temp and
  # storing the coordinates of the NA to impute
  return(list(ORDERSLGRight_temp=ORDERSLGRight_temp, nf_temp=nf_temp))
}

# 
# ################################################################################
# # Compute the CD matrix for SLG
# 
# SLGCDMatBuild <- function(CO, COt, OD, order, MaxGapSLGLeft, np, ncot, nr, nc, 
#   nf, COtsample, pastDistrib, futureDistrib, k) 
# {
#   ud <- nc - (MaxGapSLGLeft - order + np + nf)
#   
#   frameSize <- MaxGapSLGLeft - order + np + nf + 1
# 
#   CD <- matrix(NA, nrow = nr * ud, ncol = 1)
#   
#   if((np > 0 & nf > 0) & ((MaxGapSLGLeft %% 2 == 0 & order %% 2 == 0) | 
#     (MaxGapSLGLeft %% 2 != 0 & order %% 2 != 0))) {
#     shift <- MaxGapSLGLeft - order
#   } else {
#     shift <- 0
#    
#   }
# 
#   iter <- 1
#  
#   if (np > 0 & nf == 0) {
#     CD <- PastVICompute(CD, CO, OD, ncot, frameSize, iter, nr, nc, ud, np, 
#       COtsample, COt, pastDistrib, futureDistrib, k)
#   }else if(nf >0 & np==0){
#     CD <- FutureVICompute(CD, CO, OD, ncot, frameSize, iter, nr, nc, ud, np, 
#                             COtsample, COt, pastDistrib, futureDistrib, k, nf)
#   }else {
#     CD <- PastFutureVICompute(CD, CO, OD, ncot, frameSize, iter, nr, nc, ud, 
#       np, COtsample, COt, pastDistrib, futureDistrib, k, nf, shift)
#   }
#   return(list(CD, shift))
# }
# 
# 
# 
# ################################################################################
# # Impute SLG using created model
# 
# SLGCreatedModelImpute <- function(CO, OD, CD, ODi, COt, ncot, nf, nc, regr, k, 
#   reglog_6, noise, available, REFORDSLG_L, np, pastDistrib, 
#   futureDistrib, order, MaxGap, shift)
# {
#   # Structure and building of the data matrix CDi
#   # The first column of CDi is the dependent
#   # variable (VD, response variable) that we have
#   # to implement during the current iteration
#   # (i.e. automatically an NA).
#   # The following columns are the corresponding
#   # independent variables
#   # (VIs, explanatory variables)
#   # coming from the past (np>0) or the future
#   # (nf>0) (ordered by time and corresponding to
#   # the current predictive pattern) and the
#   # distribution of the possible values (i.e. all
#   # possible categorical variables numbered from 1
#   # to k and of course the value NA) respectively
#   # Before and After the NA to impute
#   # (The fact that every lines of CDi
#   # are identical is related to the working of the
#   # function mlogit that has to have as much lines
#   # in CDi as there are categories of the variable
#   # (see the parameter "k")
#   # --> so, CDi is composed of k identical lines)
#   # according to the value of np and nf
# 
# 
# 
# 
# 
#   # Analysing the value of parameter available
#   if (available == TRUE) {
#     # we take the previously imputed data
#     # into account
#     LOOKUP <- ODi
#   } else {
#     # that is available == FALSE and thus we
#     # don't take the previously imputed data
#     # into account
#     LOOKUP <- OD
#   }
# 
# 
# 
#   # Assigning the current "REFORDSLGLeft_order"
#   # matrix to the variable matrix REFORDSLGLeft
#   # (according to the current value of "order")
# 
#   REFORDSLGLeft <- as.matrix(REFORDSLG_L[[order]])
#   if (ncol(REFORDSLGLeft) == 1) {
#     REFORDSLGLeft <- t(REFORDSLGLeft)
#   }
#   nr_REFORD <- nrow(REFORDSLGLeft)
# 
# 
# 
#   if (np > 0 & nf == 0) {
#     # only PAST VIs do existe)
#     ODi <- ODiImputePAST(CO, ODi, CD, COt, REFORDSLGLeft, nr_REFORD, 
#       pastDistrib, futureDistrib, k, np, nf, nc, ncot, reglog_6, 
#       LOOKUP, regr, noise)
#   }else if(np==0 & nf>0){
#     ODi <- ODiImputeFUTURE(CO, ODi, CD, COt, REFORDSLGLeft, nr_REFORD, 
#       pastDistrib, futureDistrib, k, np, nf, nc, ncot, reglog_6, 
#       LOOKUP, regr, noise)
#   }else {
#     # meaning np>0 and nf>0 and that, thus,
#     # PAST as well as FUTURE VIs do exist
#     ODi <- ODiImputePF(CO, ODi, CD, COt, REFORDSLGLeft, nr_REFORD, pastDistrib, 
#       futureDistrib, k, np, nf, nc, ncot, reglog_6, LOOKUP, regr, 
#       noise, shift, MaxGap, order)
#   }
#   return(ODi)
# }
