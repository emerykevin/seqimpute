#' seqimpute: Imputation of missing data in longitudinal categorical data
#'
#' @description The seqimpute package implements the MICT and MICT-timing 
#' methods. These are multiple imputation methods for longitudinal data. 
#' The core idea of the algorithms is to fills gaps of missing data, which is 
#' the typical form of missing data in a longitudinal setting, recursively from 
#' their edges. The prediction is based on either a multinomial or a 
#' random forest regression model. Covariates and time-dependent covariates 
#' can be included in the model. 
#' 
#' The MICT-timing algorithm is an extension of the MICT algorithm designed 
#' to address a key limitation of the latter: its assumption that position in 
#' the trajectory is irrelevant. 
#'
#' @details The imputation process is divided into several steps, depending on
#' the type of gaps of missing data. The order of imputation of the gaps are:
#' \describe{
#'  \item{\code{Internal gap: }}{there is at least \code{np} observations 
#'  before an internal gap and \code{nf} after the gap}
#'  
#'  \item{\code{Initial gap: }}{gaps situated at the very beginning 
#'  of a trajectory}
#'  
#'  \item{\code{Terminal gap: }}{gaps situated at the very end
#'  of a trajectory}
#'  \item{\code{Left-hand side specifically located gap (SLG): }}{gaps 
#'  that have at least \code{nf} observations after the gap, but less than
#'  \code{np} observation before it}
#'  \item{\code{Right-hand side SLG: }}{gaps 
#'  that have at least \code{np} observations before the gap, but less than
#'  \code{nf} observation after it}
#'  \item{\code{Both-hand side SLG: }}{gaps 
#'  that have less than \code{np} observations before the gap, and less than
#'  \code{nf} observations after it}
#' }
#' 
#' 
#' The primary difference between the MICT and MICT-timing 
#' algorithms lies in their approach to selecting patterns from other 
#' sequences for fitting the multinomial model. While the MICT algorithm 
#' considers all similar patterns regardless of their temporal placement, 
#' MICT-timing restricts pattern selection to those that are temporally 
#' closest to the missing value. This refinement ensures that the 
#' imputation process adequately accounts for temporal dynamics, resulting 
#' in more accurate imputed values.
#'
#'
#' @param data a data frame containing sequences of a categorical
#' variable with missing data (coded as \code{NA})
#' @param var the list of columns containing the trajectories. 
#' Default is NULL, i.e. all the columns. 
#' @param np number of previous observations in the imputation model 
#' of the internal gaps.
#' @param nf number of future observations in the imputation model 
#' of the internal gaps.
#' @param m number of multiple imputations  (default: \code{5}).
#' @param timing a logical value that specifies if the MICT algorithm 
#' (timing=FALSE) or the MICT-timing algorithm (timing=TRUE) should be used.
#' @param frame.radius parameter relative to the MICT-timing algorithm 
#' specifying the radius of the timeframe.
#' @param covariates the list of columns containing the covariates to include
#' in the imputation process
#' @param time.covariates the list of columns containing the time-varying 
#' covariates to include in the imputation process
#' @param regr a character specifying the imputation method. If 
#' \code{regr="multinom"}, multinomial models are used, while 
#' if \code{regr="rf"}, random forest models are used.
#' @param npt number of previous observations in the imputation model 
#' of the terminal gaps.
#' @param nfi number of future observations in the imputation model 
#' of the initial gaps.
#' @param ParExec logical. If \code{TRUE}, the multiple imputations are run 
#' in parallel. This allows faster run time depending of how many cores 
#' the processor has.
#' @param ncores integer. Number of cores to be used for the parallel 
#' computation. If no value is set for this parameter, the number of cores 
#' will be set to the maximum number of CPU cores minus 1.
#' @param SetRNGSeed an integer that is used to set the seed in the case of 
#' parallel computation. Note that setting \code{set.seed()} alone before the 
#' seqimpute function won't work in case of parallel computation.
#' @param verbose logical. If \code{TRUE}, seqimpute will print history and 
#' warnings on console. Use \code{verbose=FALSE} for silent computation.
#' @param available a logical value allowing the user to choose whether 
#' to consider the already imputed data in the predictive model 
#' (\code{available = TRUE}) or not (\code{available = FALSE}).
#' @param pastDistrib a logical indicating if the past distribution should be 
#' used as predictor in the imputation model.
#' @param futureDistrib a logical indicating if the future distribution 
#' should be used as predictor in the imputation model.
#' @param ... Named arguments that are passed down to the imputation functions.
#'
#' @author Kevin Emery <kevin.emery@@unige.ch>, Andre Berchtold,  
#' Anthony Guinchard, and Kamyar Taher
#'
#' @return Returns an S3 object of class \code{seqimp}.
#'
#' @examples
#'
#' # Default multiple imputation of the trajectories of game addiction with the
#' # MICT algorithm
#' 
#' \dontrun{
#' set.seed(5)
#' imp1 <- seqimpute(data = gameadd, var = 1:4)
#' 
#' 
#' # Default multiple imputation with the MICT-timing algorithm
#' set.seed(3)
#' imp2 <- seqimpute(data = gameadd, var = 1:4, timing = TRUE)
#' 
#'
#' # Inclusion in the MICt-timing imputation process of the three background 
#' # characteristics (Gender, Age and Track), and the time-varying covariate 
#' # about gambling
#' 
#' 
#' set.seed(4)
#' imp3 <- seqimpute(data = gameadd, var = 1:4, covariates = 5:7, 
#'   time.covariates = 8:11)
#'
#'   
#' # Parallel computation
#' 
#' 
#' imp4 <- seqimpute(data = gameadd, var = 1:4, covariates = 5:7, 
#'   time.covariates = 8:11, ParExec = TRUE, ncores=5, SetRNGSeed = 2)
#' }
#'
#' @references HALPIN, Brendan (2012). Multiple imputation for life-course 
#' sequence data. Working Paper WP2012-01, Department of Sociology, 
#' University of Limerick. http://hdl.handle.net/10344/3639.
#' @references HALPIN, Brendan (2013). Imputing sequence data: Extensions to 
#' initial and terminal gaps, Stata's. Working Paper WP2013-01, 
#' Department of Sociology, 
#' University of Limerick. http://hdl.handle.net/10344/3620
#'
#'
#' @export
seqimpute <- function(data, var = NULL, np = 1, nf = 1, m = 5, timing = FALSE, 
  frame.radius = 0, covariates = NULL, 
  time.covariates = NULL, regr = "multinom", 
  npt = 1, nfi = 1, ParExec = FALSE, ncores = NULL, 
  SetRNGSeed = FALSE, verbose = TRUE, available = TRUE, pastDistrib = FALSE,
  futureDistrib = FALSE,...)
{
  
  call <- match.call()
  check.deprecated(...)
  covariates <- covxtract(data, covariates)
  time.covariates <- covxtract(data, time.covariates)
  
  data <- dataxtract(data, var)
  
  if (timing == FALSE) {
    imputed <- seqimpute_standard(data,
      np = np, nf = nf, m = m, covariates = covariates,
      time.covariates = time.covariates, regr = regr, nfi = nfi, npt = npt,
      available = available, pastDistrib = pastDistrib,
      futureDistrib = futureDistrib, noise = 0, ParExec = ParExec, 
      ncores = ncores, SetRNGSeed = SetRNGSeed, verbose = verbose, ...)
    method <- "MICT"
  }else {
    imputed <- seqimpute_timing(data,
      np = np, nf = nf, m = m, covariates = covariates,
      time.covariates = time.covariates, regr = regr, nfi = nfi, npt = npt,
      available = available, pastDistrib = pastDistrib,
      futureDistrib = futureDistrib, noise = 0, ParExec = ParExec, 
      ncores = ncores, SetRNGSeed = SetRNGSeed, verbose = verbose, ...)
    method <- "MICT-timing"
  }
  
  seqimpobj <- list(data = data, imp = imputed, m = m, method = method, 
      np = np, nf = nf, regr=regr, call = call)
  
  oldClass(seqimpobj) <- "seqimp"
  
  
  seqimpobj
}


seqimpute_standard <- function(data, 
  covariates = matrix(NA, nrow = 1, ncol = 1), 
  time.covariates = matrix(NA, nrow = 1, ncol = 1), np = 1, nf = 1, m = 1, 
  regr = "multinom", nfi = 1, npt = 1, available = TRUE, pastDistrib = FALSE,
  futureDistrib = FALSE, noise = 0, ParExec = FALSE, ncores = NULL,
  SetRNGSeed = FALSE, verbose = TRUE, ...)
{
  if (inherits(data, "stslist")) {
    valuesNA <- c(attr(data, "nr"), attr(data, "void"))
    data <- data.frame(data)
    data[data == valuesNA[1] | data == valuesNA[2]] <- NA
  }

  if (sum(is.na(data)) == 0) {
    if (verbose == TRUE) {
      message("This dataset has no missing values!")
    }
    return(data)
  }

  rownamesDataset <- rownames(data)
  nrowsDataset <- nrow(data)

  # 0. Initial tests and manipulations on parameters --------------------------
  dataOD <- preliminaryChecks(OD = data, CO = covariates, 
    COt = time.covariates, np = np, nf = nf, nfi = nfi, npt = npt, 
    pastDistrib = pastDistrib, futureDistrib = futureDistrib)
  
  dataOD[c("pastDistrib", "futureDistrib", "totV", "totVi", "totVt", 
    "noise")] <- InitCorectControl(regr, dataOD$ODClass, dataOD$OD, dataOD$nr, 
    dataOD$nc, dataOD$k, np, nf, dataOD$nco, dataOD$ncot, nfi, npt, 
    pastDistrib, futureDistrib, dataOD$totV, dataOD$totVi, dataOD$totVt, noise)
  
  # 1. Analysis of OD and creation of matrices ORDER, ORDER2 and ORDER3 
  dataOD[c("MaxInitGapSize", "InitGapSize", "MaxTermGapSize", "TermGapSize", 
    "MaxGap", "ORDER", "ORDER2", "ORDER3")] <- OrderCreation(dataOD$OD, 
    dataOD$nr, dataOD$nc)
  # 2. Computation of the order of imputation of each MD 
  if (max(dataOD$ORDER) != 0) {
    dataOD[c("ORDERSLGLeft", "ORDERSLGRight", "ORDERSLGBoth", "LongGap", 
      "MaxGap", "REFORD_L", "ORDER")] <- ImputeOrderComputation(dataOD$ORDER, 
      dataOD$ORDER3, dataOD$MaxGap, np, nf, dataOD$nr, dataOD$nc)
  } else {
    dataOD$ORDERSLGLeft <- matrix(nrow = dataOD$nr, ncol = dataOD$nc, 0)
    dataOD$ORDERSLGRight <- matrix(nrow = dataOD$nr, ncol = dataOD$nc, 0)
    dataOD$ORDERSLGBoth <- matrix(nrow = dataOD$nr, ncol = dataOD$nc, 0)
    dataOD$LongGap <- FALSE
  }

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
        message("/!\\ The number of multiple imputation is 1, parallel 
            processing is only available for m > 1.")
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
    
    # 3. Imputation using a specific model --------------------------------
    if (max(dataOD$ORDER) != 0) {
      # Otherwise if there is only 0 in ORDER,
      # there is no need to impute internal gaps
      # and we directly jump to the imputation of
      # external gaps (i.e. points 4. and 5.)
      if (verbose == TRUE) {
        print("Imputation of the internal gaps...")
      }
      dataOD[["ODi"]] <- ModelImputation(OD = dataOD$OD, 
          covariates = dataOD$CO, time.covariates = dataOD$COt, 
          ODi = dataOD$ODi, MaxGap = dataOD$MaxGap, totV = dataOD$totV, 
          totVi = dataOD$totVi, regr = regr, nc = dataOD$nc, np = np, nf = nf, 
          nr = dataOD$nr, ncot = dataOD$ncot, COtsample = dataOD$COtsample, 
          pastDistrib = dataOD$pastDistrib, 
          futureDistrib = dataOD$futureDistrib, k = dataOD$k, 
          available = available, REFORD_L = dataOD$REFORD_L, 
          noise = dataOD$noise, verbose,...)
    }
    # 4. Imputing initial NAs ---------------------------------------------
    if ((nfi != 0) & (dataOD$MaxInitGapSize != 0)){ 
      if (verbose == TRUE) {
        print("Imputation of the initial gaps...")
      }
      # # we only impute the initial gaps if nfi > 0
      dataOD[["ODi"]] <- ImputingInitialNAs(OD = dataOD$OD, 
          covariates = dataOD$CO, time.covariates = dataOD$COt, 
          ODi = dataOD$ODi, totVi = dataOD$totVi, COtsample = dataOD$COtsample,
          futureDistrib = dataOD$futureDistrib, 
          InitGapSize = dataOD$InitGapSize, 
          MaxInitGapSize = dataOD$MaxInitGapSize, nr = dataOD$nr, 
          nc = dataOD$nc, ud = dataOD$ud, nco = dataOD$nco, 
          ncot = dataOD$ncot, nfi = nfi, regr = regr, k = dataOD$k, 
          available = available, noise = dataOD$noise, ...)
    }
    # 5. Imputing terminal NAs ----------------------------------------------
    if ((npt != 0) & (dataOD$MaxTermGapSize != 0)) {
      # we only impute the terminal
      # gaps if npt > 0
      if (verbose == TRUE) {
        print("Imputation of the terminal gaps...")
      }
      dataOD[["ODi"]] <- ImputingTerminalNAs(OD = dataOD$OD, 
          covariates = dataOD$CO, time.covariates = dataOD$COt, 
          ODi = dataOD$ODi, COtsample = dataOD$COtsample, 
          MaxTermGapSize = dataOD$MaxTermGapSize, 
          TermGapSize = dataOD$TermGapSize, pastDistrib = dataOD$pastDistrib, 
          regr = regr, npt = npt, nco = dataOD$nco, ncot = dataOD$ncot, 
          totVt = dataOD$totVt, nr = dataOD$nr, nc = dataOD$nc, ud = dataOD$ud, 
          available = available, k = dataOD$k, noise = dataOD$noise, ...)
    }
    # 6. Imputing SLG NAs --------------------------------------------------
      # Checking if we have to impute
      # left-hand side SLG
    if (max(dataOD$ORDERSLGLeft) != 0) { 
      if (verbose == TRUE) {
        print("Imputation of the left-hand side SLG...")
      }
      dataOD[["ODi"]] <- LSLGNAsImpute(OD = dataOD$OD, ODi = dataOD$ODi, 
          covariates = dataOD$CO, time.covariates = dataOD$COt, 
          COtsample = dataOD$COtsample, ORDERSLG = dataOD$ORDERSLGLeft,
          pastDistrib = dataOD$pastDistrib, 
          futureDistrib = dataOD$futureDistrib, regr = regr, np = np, 
          nr = dataOD$nr, nf = nf, nc = dataOD$nc, 
          ud = dataOD$ud, ncot = dataOD$ncot,nco = dataOD$nco, k = dataOD$k, 
          noise = dataOD$noise, available = available, ...)
    }
    # right-hand side SLG
    if (max(dataOD$ORDERSLGRight) != 0) {
      # Checking if we have to impute right-hand
      # side SLG
      if (verbose == TRUE) {
        print("Imputation of the right-hand side SLG...")
      }
      dataOD[["ODi"]] <- RSLGNAsImpute(OD = dataOD$OD, ODi = dataOD$ODi, 
          covariates = dataOD$CO, time.covariates = dataOD$COt, 
          COtsample = dataOD$COtsample, ORDERSLGRight = dataOD$ORDERSLGRight,
          pastDistrib = dataOD$pastDistrib, 
          futureDistrib = dataOD$futureDistrib, regr = regr, np = np, 
          nr = dataOD$nr, nf = nf, nc = dataOD$nc, 
          ud = dataOD$ud, ncot = dataOD$ncot,nco = dataOD$nco, k = dataOD$k, 
          noise = dataOD$noise, available = available, ...)
    }
    # Checking if we have to impute
    # Both-hand side SLG
    if (dataOD$LongGap){
      if (verbose == TRUE) {
        print("Imputation of the both-hand side SLG...")
      }
      for (h in 2:np) {
        if (sum(dataOD$ORDERSLGBoth[, h - 1] == 0 & 
                dataOD$ORDERSLGBoth[, h] != 0) > 0) {
          tt <- which(dataOD$ORDERSLGBoth[, h - 1] == 0 & 
                        dataOD$ORDERSLGBoth[, h] != 0)
          tmpORDER <- matrix(0, nrow(dataOD$ORDERSLGBoth),
            ncol(dataOD$ORDERSLGBoth))
          tmpORDER[tt, h:ncol(dataOD$ORDERSLGBoth)] <- dataOD$ORDERSLGBoth[tt, 
              h:ncol(dataOD$ORDERSLGBoth)]

          dataOD[["ODi"]] <- RSLGNAsImpute(OD = dataOD$OD, ODi = dataOD$ODi, 
              covariates = dataOD$CO, time.covariates = dataOD$COt, 
              COtsample = dataOD$COtsample, ORDERSLGRight = tmpORDER,
              pastDistrib = dataOD$pastDistrib, 
              futureDistrib = dataOD$futureDistrib, regr = regr, np = h - 1, 
              nr = dataOD$nr, nf = nf, nc = dataOD$nc, ud = dataOD$ud, 
              ncot = dataOD$ncot, nco = dataOD$nco, k = dataOD$k, 
              noise = dataOD$noise, available = available, ...)
        }
      }
    }

    return(dataOD$ODi)
  }
  if (ParParams) {
    parallel::stopCluster(cl)
  }
  
  names(RESULT) <- paste0("imp",1:m)

  # RESULT <- rbind(cbind(replicate(dataOD$nr, 0), dataOD$OD), RESULT)

  # X. Final conversions -----------------------------------------------------
  RESULT <- lapply(RESULT,FinalResultConvert, ODClass = dataOD$ODClass,
      ODlevels = dataOD$ODlevels, rownamesDataset = rownamesDataset, 
      nrowsDataset = nrowsDataset, nr = dataOD$nr, nc = dataOD$nc, 
      rowsNA = dataOD$rowsNA, mi = m)
  return(RESULT)

}
