InitCorectControl <- function(regr, ODClass, OD, nr, nc, k, np, nf, nco, ncot,
                              nfi, npt, pastDistrib, futureDistrib, totV, totVi, totVt, noise) {
  # 0.1 Testing "regr" effectively either "mlogit", "lm" or "lrm" ------------
  if((regr != "rf") & (regr != "multinom")){
    stop("/!\\ regr defines the type of regression model you want to use.
               It has to be either assigned to character 'multinom'
               (for multinomialregression) or'rf' (for random forests)")
  }



  # 0.4 Eventually discarding the consideration of pastDistrib and futureDistrib
  # Making sure that pastDistrib and futureDistrib are set to FALSE by default
  # in case OD is made of "numeric" variables

  ####### numeric ########
  # if (ODClass == "numeric") {
  #   pastDistrib <- FALSE
  #   futureDistrib <- FALSE
  #   # Update of the totV variables
  #   # (which become then smaller since
  #   # we don't consider pastDistrib and futureDistrib)
  #   totV <- 1 + np + nf + nco + (ncot / nc)
  #   totVi <- 1 + nfi + nco + (ncot / nc)
  #   totVt <- 1 + npt + nco + (ncot / nc)
  # }









  # 0.5 Testing not np<0, nor nf<0 nor both ==0 ------------------------------
  # No VIs is not possible (we should have at least np>0 or nf>0)
  if (np == 0 & nf == 0) {
    stop("/!\\ We can't have np as well as nf equal to '0' at the same
               time")
  }
  # Negative value for np as well as nf raises an error
  if (np < 0 | nf < 0) {
    stop("/!\\ np and nf can't be negative numbers")
  }



  # 0.7 Testing not nfi<0, nor npt<0 ---------------------------------------
  # Negative value for nfi as well as npt raises an error
  if (nfi < 0 | npt < 0) {
    stop("/!\\ nfi and npt can't be negative numbers")
  }


  # 0.8 Testing the right construction of COt ---------------------------------
  # Since COt contains the time-dependent covariates, each of them is
  # represented by a submatrix that has the same number of columns as OD.
  # The total number of columns of COt is necessarily a multiple of the number
  # of column of OD
  if (ncot %% nc != 0) {
    stop("/!\\ Each time-dependent covariates contained in COt has to have the
      same number of columns as the dataset.")
  }


  # 0.9 Taking the absolute value of the parameter "noise" --------------------
  # Since "noise" is the variance of the elements of the final vector pred, it
  # can't be negative
  noise <- abs(noise)

  return(list(pastDistrib, futureDistrib, totV, totVi, totVt, noise))
}


check.deprecated <- function(...) {
  nms <- names(list(...))
  replace.args <- list(CO = "covariates", COt = "time.covariates", OD = "data")
  wrn <- names(replace.args) %in% nms
  if (any(wrn)) {
    for (i in which(wrn)) {
      msg <- paste0(
        "The '", names(replace.args)[i],
        "' argument is no longer supported. Please use '",
        replace.args[i], "' instead.")
      warning(msg)
    }
  }

  deleted.args <- c("mice.return", "include")
  wrn <- deleted.args %in% nms
  if (any(wrn)) {
    for (i in which(wrn)) {
      msg <- paste0(
        "The '", deleted.args[i], "' argument is no longer supported. An object
        of class seqimp is returned by the function since version 2.0"
      )
      warning(msg)
    }
  }

  invisible(NULL)
}


dataxtract <- function(data, var)
{
  if (inherits(data, "tbl_df"))
    data <- as.data.frame(data)
  if (missing(var) || is.null(var) || is.na(var[1]))
    seqdata <- data
  else seqdata <- subset(data, , var)
  return(seqdata)
}

dataputback <- function(data, var, data.traj){
  if (missing(var) || is.null(var) || is.na(var[1]))
    seqdata <- data
  else{
    data[,var] <- data.traj
  }
  return(data)
}

covxtract <- function(data, covariates)
{
  if(!inherits(covariates,"data.frame")){
    if(missing(covariates) || is.null(covariates) || is.na(covariates[1])){
      data.cov <- matrix(NA, nrow = 1, ncol = 1)

    }else if(length(covariates)==nrow(data) & !covariates[1]%in%colnames(data)){
      data.cov <- covariates
    }else{
        data.cov <- subset(data, , covariates)
      }
  }else{
    data.cov <- covariates
  }
  return(data.cov)
}

preliminaryChecks <- function(OD, CO, COt, var){

  # if (inherits(OD, "stslist")) {
  #   valuesNA <- c(attr(OD, "nr"), attr(OD, "void"))
  #   OD <- data.frame(OD)
  #   OD[OD == valuesNA[1] | OD == valuesNA[2]] <- NA
  # }else{
  #   CO <- covxtract(OD, CO)
  #   COt <- covxtract(OD, COt)
  # 
  #   OD <- dataxtract(OD, var)
  # }

  data <- list()
  # Updating the number of columns of CO
  data["nco"] <- compute.ncol(CO)
  # Updating the number of columns of COt
  data["ncot"] <- compute.ncol(COt)

  data$ncot <- check.ncot(data$ncot,ncol(OD))


  # Deleting entire rows of OD filled only with NAs
  # (and deleting corresponding lines in CO and COt)
  data$rowsNA <- rows.allmiss(OD)
  if(length(data$rowsNA)>0){
    data$OD <- OD[-data$rowsNA,]
    if(data$nco>0){
      data$CO <- CO[-data$rowsNA,]
    }
    if(data$ncot>0){
      data$COt <- COt[-data$rowsNA,]
    }
  }else{
    data$OD <- OD
    data$CO <- CO
    data$COt <- COt
  }
  data[c("OD", "ODi", "ODClass", "ODlevels", "k",
    "nr", "nc")] <- factorAndMatrixChecks(data$OD)

  data$COtsample <- compute.COtsample(data$COt, data$ncot, data$nr, data$nc)

  return(data)
}


check.predictors <- function(np, nf, npt, nfi){
  if (np == 0 & nf == 0) {
    stop("/!\\ We can't have np as well as nf equal to '0' at the same
               time")
  }
  if (np < 0 | nf < 0) {
    stop("/!\\ np and nf can't be negative numbers")
  }
  if (nfi < 0 | npt < 0) {
    stop("/!\\ nfi and npt can't be negative numbers")
  }
  return(list(np=np, nf=nf, npt=npt, nfi=nfi))
}

check.regr <- function(regr){
  if((regr != "rf") & (regr != "multinom")){
    stop("/!\\ regr defines the type of regression model you want to use.
               It has to be either assigned to character 'multinom'
               (for multinomialregression) or'rf' (for random forests)")
  }
  return(regr)
}

check.ncot <- function(ncot,nc){
  if (ncot %% nc != 0) {
    stop("/!\\ Each time-dependent covariates contained in COt has to have the
      same number of columns as the dataset.")
  }
  return(ncot)
}

compute.ncol <- function(x){
  if (all(is.na(x)) == FALSE){
    if (is.null(dim(x))) {
      return(1)
    } else {
      return(ncol(x))
    }
  }
  return(0)
}

rows.allmiss <- function(OD){
  rowsNA <- c()
  for (i in 1:nrow(OD)) {
    if (all(is.na(OD[i, ]))) {
      rowsNA <- c(rowsNA, i)
    }
  }
  return(rowsNA)
}


deleteNaRows <- function(OD, CO, COt) {
  rowsNA <- c()
  for (i in 1:nrow(OD)) {
    if (all(is.na(OD[i, ]))) {
      rowsNA <- c(rowsNA, i)
    }
  }

  if (length(rowsNA) > 0) {
    OD <- OD[-rowsNA, ]
    if (all(is.na(CO)) == FALSE) {
      if (is.null(dim(CO))) {
        CO <- CO[-rowsNA]
      } else {
        CO <- CO[-rowsNA, ]
      }
    }
    if (all(is.na(COt)) == FALSE) {
      COt <- COt[-rowsNA, ]
    }
  }
  return(list(OD, CO, COt, rowsNA))
}

compute.COtsample <- function(COt, ncot, nr, nc)
{
  COtsample <- vector()
  if (ncot > 0) {

    COtsample <- as.data.frame(matrix(nrow = nr, ncol = 0))

    for (d in 1:(ncot / nc)) {
      COtsample <- cbind(COtsample, COt[, 1 + (d - 1) * nc])
    }
  }

  return(COtsample)
}
factorAndMatrixChecks <- function(OD) {
  nc <- ncol(OD)
  nr <- nrow(OD)
  ODClass <- class(OD[1, 1])

  if ((ODClass != "factor") & (ODClass != "character")) {
    stop("/!\\ The class of the variables contained in your original dataset
           should be either 'factor' or 'character'")
  }

  #*************************************
  ODlevels <- vector()
  
  ODlevels <- sort(unique(as.vector(as.matrix(OD))))
  k <- length(ODlevels)
  OD <- as.data.frame(sapply(OD, mapvalues,
                             from = ODlevels,
                             to = as.character(as.vector(1:length(ODlevels)))
  ))
 
  OD <- apply(as.matrix(OD), 2, as.numeric)
  

  ODi <- OD

  if (ODClass == "factor") {
    for (j in 1:nc) {
      ODi[, j] <- factor(x = OD[, j], levels = c(1:k))
    }
  }

  return(list(OD, ODi, ODClass, ODlevels, k, nr, nc))
}


check.cores <- function(ncores, available, m){
  if(is.null(ncores)) {
    ncores <- min(available - 1, m)
  }
  else {
    if(ncores > available){
      warning(paste("'ncores' exceeds the maximum number of available cores on
                    your machine, and is set to",
                    min(available - 1, m)))
    }

    if(ncores > m){
      warning(paste("'ncores' exceeds the number of imputations, and is set to",
                    min(available - 1, m)))
    }

    ncores <- min(available - 1, m, ncores)
  }
  ncores
}
