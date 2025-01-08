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
#' imputation process adequately accounts for temporal dynamics, imping 
#' in more accurate imputed values.
#'
#'
#' @param data Either a data frame containing sequences of a categorical 
#' variable, where missing data are coded as \code{NA}, or a state sequence 
#' object created using the \link[TraMineR]{seqdef} function. If using a 
#' state sequence object, any "void" elements will also be treated as missing. 
#' See the \code{end.impute} argument if you wish to skip imputing values 
#' at the end of the sequences.
#' @param var A specifying the columns of the dataset 
#' that contain the trajectories. Default is \code{NULL}, meaning all columns 
#' are used.
#' @param np Number of prior states to include in the imputation model 
#' for internal gaps.
#' @param nf Number of subsequent states to include in the imputation model 
#' for internal gaps.
#' @param m Number of multiple imputations to perform (default: \code{5}).
#' @param timing Logical, specifies the imputation algorithm to use. 
#' If \code{FALSE}, the MICT algorithm is applied; if \code{TRUE}, the 
#' MICT-timing algorithm is used.
#' @param frame.radius Integer, relevant only for the MICT-timing algorithm, 
#' specifying the radius of the timeframe.
#' 
#' @param covariates List of the columns of the dataset
#' containing covariates to be included in the imputation model.
#' 
#' @param time.covariates List of the columns of the dataset
#'  with time-varying covariates to include in the imputation model.
#' 
#' @param regr Character specifying the imputation method. Options include 
#' \code{"multinom"} for multinomial models and \code{"rf"} for random forest 
#' models.
#' 
#' @param npt Number of prior observations in the imputation model for 
#' terminal gaps (i.e., gaps at the end of sequences).
#' 
#' @param nfi Number of future observations in the imputation model for 
#' initial gaps (i.e., gaps at the beginning of sequences).
#' 
#' @param ParExec Logical, indicating whether to run multiple imputations 
#' in parallel. Setting to \code{TRUE} can improve computation time depending 
#' on available cores.
#' 
#' @param ncores Integer, specifying the number of cores to use for parallel 
#' computation. If unset, defaults to the maximum number of CPU cores minus one.
#' 
#' @param SetRNGSeed Integer, to set the random seed for reproducibility in 
#' parallel computations. Note that setting \code{set.seed()} alone does not 
#' ensure reproducibility in parallel mode.
#' 
#' 
#' @param end.impute Logical. If \code{FALSE}, missing data at the end of 
#' sequences will not be imputed.
#' 
#' @param verbose Logical, if \code{TRUE}, displays progress and warnings 
#' in the console. Use \code{FALSE} for silent computation.
#' 
#' @param available Logical, specifies whether to consider already imputed 
#' data in the predictive model. If \code{TRUE}, previous imputations are 
#' used; if \code{FALSE}, only original data are considered.
#' 
#' @param pastDistrib Logical, if \code{TRUE}, includes the past distribution 
#' as a predictor in the imputation model.
#' 
#' @param futureDistrib Logical, if \code{TRUE}, includes the future 
#' distribution as a predictor in the imputation model.
#' @param ... Named arguments that are passed down to the imputation functions.
#'
#' @author Kevin Emery <kevin.emery@@unige.ch>, Andre Berchtold,  
#' Anthony Guinchard, and Kamyar Taher
#'
#' @return An object of class \code{seqimp}, which is a list with the following 
#' elements:
#' \describe{
#'   \item{\code{data}}{A \code{data.frame} containing the original 
#'   (incomplete) data.}
#'   \item{\code{imp}}{A list of \code{m} \code{data.frame} corresponding to 
#'   the imputed datasets.}
#'   \item{\code{m}}{The number of imputations.}
#'   \item{\code{method}}{A character vector specifying whether MICT or 
#'   MICT-timing was used.}
#'   \item{\code{np}}{Number of prior states included in the imputation model.}
#'   \item{\code{nf}}{Number of subsequent states included in the imputation 
#'   model.}
#'   \item{\code{regr}}{A character vector specifying whether multinomial or
#'   random forest imputation models were applied.}
#'   \item{\code{call}}{The call that created the object.}
#' }
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
#' @references Halpin, B. (2012). Multiple imputation for life-course 
#' sequence data. Working Paper WP2012-01, Department of Sociology, 
#' University of Limerick. http://hdl.handle.net/10344/3639.
#' @references Halpin, B. (2013). Imputing sequence data: Extensions to 
#' initial and terminal gaps, Stata's. Working Paper WP2013-01, 
#' Department of Sociology, 
#' University of Limerick. http://hdl.handle.net/10344/3620
#' @references Emery, K., Studer, M., & Berchtold, A. (2024). Comparison of 
#' imputation methods for univariate categorical longitudinal data. 
#' Quality & Quantity, 1-25. 
#' https://link.springer.com/article/10.1007/s11135-024-02028-z
#'
#' @export
seqimpute <- function(data, var = NULL, np = 1, nf = 1, m = 5, timing = FALSE, 
                      frame.radius = 0, covariates = NULL, 
                      time.covariates = NULL, regr = "multinom",npt = 1, 
                      nfi = 1, ParExec = FALSE, ncores = NULL, 
                      SetRNGSeed = FALSE, end.impute = TRUE, verbose = TRUE,
                      available = TRUE, pastDistrib = FALSE, 
                      futureDistrib = FALSE,...)
{
  
  call <- match.call()
  check.deprecated(...)
  
  if (inherits(data, "stslist")) {
    valuesNA <- c(attr(data, "nr"), attr(data, "void"))
    data <- data.frame(data)
    data[data == valuesNA[1] | data == valuesNA[2]] <- NA
  }else{
    covariates <- covxtract(data, covariates)
    time.covariates <- covxtract(data, time.covariates)
    
    data <- dataxtract(data, var)
  }
  dataOD <- check.data(data, covariates, time.covariates, var)
  if (sum(is.na(dataOD$OD)) == 0) {
    if (verbose == TRUE) {
      message("This dataset has no missing values!")
    }
    return(dataOD$OD)
  }
  
  tmp <- check.predictors(np, nf, nfi, npt)
  np <- tmp$np
  nf <- tmp$nf
  nfi <- tmp$nfi
  npt <- tmp$npt
  
  regr <- check.regr(regr)
  
  imporder <- compute.order(dataOD$OD, dataOD$nr, dataOD$nc, np, nf, npt, nfi, 
                            end.impute)
  
  if(ParExec){
    available.cores <- parallelly::availableCores(logical = TRUE)
    ncores <- check.cores(ncores, available.cores, m)
  }else{
    ncores <- 1
  }
  
  if(ncores > 1){
    cl <- parallel::makeCluster(ncores)
    doSNOW::registerDoSNOW(cl) 
    if (SetRNGSeed) {
      doRNG::registerDoRNG(SetRNGSeed)
    }
    pb <- txtProgressBar(max = m, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    ParParams <- TRUE
  }else{
    if (SetRNGSeed) {
      set.seed(SetRNGSeed)
    }
    
    foreach::registerDoSEQ()
    opts <- NULL
    
    ParParams <- FALSE
  }
  
  o <- NULL
  imp <- foreach(o = 1:m, .inorder = TRUE, 
                    .options.snow = opts) %dopar% {
                      
                      if (!ParParams) {
                        if (verbose == TRUE) {
                          cat("iteration :", o, "/", m, "\n")
                        }
                      }
                      
                      
                      if (timing == FALSE) {
                        imp <- mict(dataOD, imporder = imporder,
                                      np = np, nf = nf, m = m, regr = regr, 
                                      nfi = nfi, npt = npt,
                                      available = available, 
                                      pastDistrib = pastDistrib,
                                      futureDistrib = futureDistrib,
                                      verbose = verbose,...)
                      }else {
                        imp <- mict.timing(dataOD, imporder = imporder,
                                          np = np, nf = nf, m = m, regr = regr, 
                                          nfi = nfi, npt = npt,
                                          available = available, 
                                          pastDistrib = pastDistrib,
                                          futureDistrib = futureDistrib,
                                          verbose = verbose, 
                                          frame.radius=frame.radius,...)
                      }
                    }
  if (ParParams) {
    parallel::stopCluster(cl)
  }
  names(imp) <- paste0("imp",1:m)
  
 
  imp <- lapply(imp, final.transform, ODClass = dataOD$ODClass,
                  ODlevels = dataOD$ODlevels, 
                  rownamesDataset = rownames(dataOD$OD), 
                  nrowsDataset = nrow(dataOD$OD), nr = dataOD$nr, 
                  nc = dataOD$nc, rowsNA = dataOD$rowsNA, mi = m)
  
  if(timing==TRUE){
    method <- "MICT-timing"
  }else{
    method <- "MICT"
    
  }
  seqimpobj <- list(data = data, imp = imp, m = m, method = method, 
                    np = np, nf = nf, regr = regr, call = call)
  
  oldClass(seqimpobj) <- "seqimp"
  
  
  seqimpobj
}


mict <- function(dataOD, imporder=NULL, np = 1, nf = 1, m = 1, 
  regr = "multinom", nfi = 1, npt = 1, available = TRUE, pastDistrib = FALSE,
  futureDistrib = FALSE, verbose = TRUE, ...)
{
  imp <- dataOD$ODi
  noise <- 0
    
    if (imporder$maxInternal != 0) {
      
      if (verbose == TRUE) {
        print("Imputation of the internal gaps...")
      }
      

     imp <- mict.internal(data=dataOD, imp, MaxGap = imporder$maxInternal, 
                            regr = regr, nc = dataOD$nc, np = np, nf = nf, 
                            nr = dataOD$nr, ncot = dataOD$ncot,
                            pastDistrib = pastDistrib, 
                            futureDistrib = futureDistrib, k = dataOD$k, 
                            available = available, 
                            REFORD_L = imporder$REFORD_L, noise = dataOD$noise, 
                            verbose=verbose,...)
    }
    if (imporder$maxInitial != 0){ 
      if (verbose == TRUE) {
        print("Imputation of the initial gaps...")
      }
      # # we only impute the initial gaps if nfi > 0
     imp <- mict.initial(dataOD, imp, futureDistrib = futureDistrib, 
          REFORDI_L = imporder$initial, 
          MaxInitGapSize = imporder$maxInitial, nr = dataOD$nr, 
          nc = dataOD$nc, ud = dataOD$ud, nco = dataOD$nco, 
          ncot = dataOD$ncot, nfi = nfi, regr = regr, k = dataOD$k, 
          available = available, noise = dataOD$noise, ...)
    }
    if (imporder$maxTerminal != 0) {
      # we only impute the terminal
      # gaps if npt > 0
      if (verbose == TRUE) {
        print("Imputation of the terminal gaps...")
      }
     imp <- mict.terminal(dataOD, imp, MaxTermGapSize = imporder$maxTerminal, 
          REFORDT_L = imporder$terminal, pastDistrib = pastDistrib, 
          regr = regr, npt = npt, nco = dataOD$nco, ncot = dataOD$ncot, 
          nr = dataOD$nr, nc = dataOD$nc, ud = dataOD$ud, 
          available = available, k = dataOD$k, noise = dataOD$noise, ...)
    }
      # Checking if we have to impute
      # left-hand side SLG
  if (max(imporder$maxLeftSLG)>0) {
    if (verbose == TRUE) {
        print("Imputation of the left-hand side SLG...")
      }
     imp <- mict.leftSLG(dataOD, imp, pastDistrib = pastDistrib, 
          futureDistrib = futureDistrib, regr = regr, np = np, 
          nr = dataOD$nr, nf = nf, nc = dataOD$nc, 
          ud = dataOD$ud, ncot = dataOD$ncot,nco = dataOD$nco, k = dataOD$k, 
          noise = dataOD$noise, available = available, REFORD_L = imporder$REFORDSLGLeft,
          MaxGap=imporder$maxLeftSLG,...)
    }
    # right-hand side SLG
  if (max(imporder$maxRightSLG)>0) {
    if (verbose == TRUE) {
        print("Imputation of the right-hand side SLG...")
      }
     imp <- mict.rightSLG(dataOD, imp, pastDistrib = pastDistrib, 
          futureDistrib = futureDistrib, regr = regr, np = np, 
          nr = dataOD$nr, nf = nf, nc = dataOD$nc, 
          ud = dataOD$ud, ncot = dataOD$ncot,nco = dataOD$nco, k = dataOD$k, 
          noise = dataOD$noise, available = available, REFORD_L = imporder$REFORDSLGRight,
          MaxGap=imporder$maxRightSLG,...)
    }
  if (max(imporder$maxBothSLG)>0) {
    if (verbose == TRUE) {
        print("Imputation of the both-hand side SLG...")
      }
      for (h in 2:np) {
        if (sum(imporder$ORDERSLGBoth[, h - 1] == 0 & 
                imporder$ORDERSLGBoth[, h] != 0) > 0) {
          
         imp <- mict.rightSLG(dataOD, imp, pastDistrib = pastDistrib, 
              futureDistrib = futureDistrib, regr = regr, np = h - 1, 
              nr = dataOD$nr, nf = nf, nc = dataOD$nc, ud = dataOD$ud, 
              ncot = dataOD$ncot, nco = dataOD$nco, k = dataOD$k, 
              noise = dataOD$noise, available = available, REFORD_L = imporder$REFORDSLGBoth[[h]],
              MaxGap=imporder$maxBothSLG[h,],...)
        }
      }
    }

    return(imp)

}
