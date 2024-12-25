# 3. Imputation using a specific model -----------------------------------------


################################################################################
# Impute data using a specific model

ModelImputation <- function(OD, covariates, time.covariates, ODi, MaxGap, 
  regr, nc, np, nf, nr, ncot, COtsample, pastDistrib, futureDistrib, k, 
  available, REFORD_L, noise, verbose, ...)
{
  for (order in 1:MaxGap) {
    if(verbose==TRUE){
      print(paste0("Step ", order, "/", MaxGap))
    }
    # /!\ "order" corresponds to the
    # values of the components of ORDER (i.e. the number
    # of the iteration, the order in which the
    # values are going to be imputed)
    # 3.1. Building of the data matrix CD for the computation of the model ----
    # CD_shift <- CDComputeTiming(CO=covariates, OD=OD, COt=time.covariates, 
    #                     MaxGap=MaxGap, order=order, np=np, nc=nc, nr=nr, nf=nf, COtsample, 
    #                     pastDistrib=pastDistrib, futureDistrib=futureDistrib, 
    #                     ncot=ncot, k=k, col=np+1, frame.radius=nc)
    
    CD_shift <- CDCompute(covariates, OD, time.covariates, MaxGap, order, np, 
                          nc, nr, nf, COtsample, pastDistrib, futureDistrib, ncot, k)
    if (length(table(CD_shift$CD[, 1])) > 1) {
      # 3.2. Computation of the model (Dealing with the LOCATIONS of imputation
      log_CD <- list()
      log_CD[c("reglog", "CD")] <- ComputeModel(CD_shift$CD, regr, np, 
                                                nf, k, ...)
      # 3.3. Imputation using the just created model 
      ODi <- CreatedModelImputation(order, covariates, log_CD$CD, 
        time.covariates, OD, ODi, pastDistrib, futureDistrib, available, 
        REFORD_L, ncot, nc, np, nf, k, regr, log_CD$reglog, noise, 
        CD_shift$shift, MaxGap)
    }else{
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
  return(ODi)
}

################################################################################
# Building of the data matrix CD for the computation of the model

CDCompute <- function(CO, OD, COt, MaxGap, order, np, nc, nr, nf, COtsample, 
  pastDistrib, futureDistrib, ncot, k) {
  # Building of a data matrix for the computation of the
  # model
  ud <- nc - (MaxGap - order + np + nf) # number of usable data


  # for each row of OD
  frameSize <- MaxGap - order + np + nf + 1 # size of the current
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
  
  # the value of np and nf
  CD <- matrix(NA, nrow = nr * ud, ncol = 1)
  # initialization of
  # the current very left part of
  # the predictive model matrix
  # ("matrice de modele de
  # prediction") with NA
  # everywhere (/!\ for each
  # "order", we are going to
  # build such a CD)

  # Dealing with the change of shape of the prediction
  # frame (according to whether the imputed data is
  # located at the beginning (left) of a gap or at the end
  # (right)).
  # The purpose of this if statement is to detect if we
  # have to insert a shift (to jump at the end of the gap)
  # or not
  if ((np > 0 & nf > 0) & ((MaxGap %% 2 == 0 & order %% 2 == 0) | 
    (MaxGap %% 2 != 0 & order %% 2 != 0))) {
    shift <- MaxGap - order # jumping at the end of
    # the gap
  } else {
    shift <- 0 # no shift is needed (note that
    # no shift is needed for the
    # case of model 2 (only past)
    # and model 3 (only future))
  }

  iter <- 1 # initialisation of the number
  # of iterations of the
  # following for loops
  # Only PAST
  if (np > 0 & nf == 0) {
    CD <- PastVICompute(CD, CO, OD, ncot, frameSize, iter, nr, nc, ud, np, 
      COtsample, COt, pastDistrib, futureDistrib, k)
    # Only FUTURE
  } else if (np == 0 & nf > 0) {
    # only FUTURE VIs do exist
    CD <- FutureVICompute(CD, CO, OD, ncot, frameSize, iter, nr, nc, ud, np, 
      COtsample, COt, pastDistrib, futureDistrib, k, nf)

    # PAST and FUTURE
  } else {
    # meaning np>0 and nf>0 and that, thus,
    # PAST as well as FUTURE VIs do exist
    CD <- PastFutureVICompute(CD, CO, OD, ncot, frameSize, iter, nr, nc, ud, 
      np, COtsample, COt, pastDistrib, futureDistrib, k, nf, shift)
  }

  CD_shift <- list()
  CD_shift[c("CD", "shift")] <- list(CD, shift)
  return(CD_shift)
}



################################################################################
# Imputation using the just created model 

CreatedModelImputation <- function(order, CO, CD, COt, OD, ODi, pastDistrib, 
  futureDistrib, available, REFORD_L, ncot, nc, np, nf, k, regr, reglog, 
  noise, shift, MaxGap)
{
  
  if (available == TRUE) { # we take the previously imputed
    # data into account
    LOOKUP <- ODi
  } else { # that is available == FALSE and thus we
    # don't take the previously imputed
    # data into account
    LOOKUP <- OD
  }


  REFORD <- as.matrix(REFORD_L[[order]])
  if (ncol(REFORD) == 1) {
    REFORD <- t(REFORD)
  }
  nr_REFORD <- nrow(REFORD)
  if (np > 0 & nf == 0) { # only PAST VIs do exist
    ODi <- ODiImputePAST(CO, ODi, CD, COt, REFORD, nr_REFORD, pastDistrib, 
      futureDistrib, k, np, nf, nc, ncot, reglog, LOOKUP, regr, noise)

    #-------------------------------------------------------------------------
  } else if (np == 0 & nf > 0) { # only FUTURE VIs do exist
    ODi <- ODiImputeFUTURE(CO, ODi, CD, COt, REFORD, nr_REFORD, pastDistrib, 
      futureDistrib, k, np, nf, nc, ncot, reglog, LOOKUP, regr, noise)
  } else { # meaning np>0 and nf>0 and that,
    # thus, PAST as well as FUTURE VIs
    # do exist
    ODi <- ODiImputePF(CO, ODi, CD, COt, REFORD, nr_REFORD, pastDistrib, 
      futureDistrib, k, np, nf, nc, ncot, reglog, LOOKUP, regr, noise, 
      shift, MaxGap, order)
  }
  return(ODi)
}



################################################################################
# Imputation where only FUTURE VIs exist

ODiImputeFUTURE <- function(CO, ODi, CD, COt, REFORD, nr_REFORD, pastDistrib, 
  futureDistrib, k, np, nf, nc, ncot, reglog, LOOKUP, regr, noise)
{
  for (u in 1:nr_REFORD) {
    i <- REFORD[u, 1]
    # taking out the first coordinate
    # (row number in ORDER) from REFORD
    j <- REFORD[u, 2]
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
    
    
    # if (regr == "lm" | regr == "lrm") {
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
    #   # multinom
    #   CDi[, 1] <- factor(CDi[, 1], levels = c(1:k))
    #   for (v in 2:(1 + np + nf)) {
    #     CDi[, v] <- factor(CDi[, v], levels = levels(CD[, v]), exclude = NULL)
    #   }
    # }
    # The last values of CDi must be of type numeric
    # (distributions)
    if (pastDistrib | futureDistrib) {
      CDi[,(1+np+nf+1):ncol(CDi)] <- lapply(CDi[,(1+np+nf+1):ncol(CDi)], 
        as.numeric)
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
    if (max(is.na(CDi[1, 2:ncol(CDi)])) == 0) {
      ODi <- RegrImpute(ODi, CDi, regr, reglog, noise, i, j, k)
    }
  }
  return(ODi)
}
