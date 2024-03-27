#' Plot the most common patterns of missing data.
#' 
#' @description
#' Plot function that renders the most frequent patterns of missing data. This
#' function is based on the \link[TraMineR]{seqfplot} function. 
#' 
#'
#' @param data a data.frame where missing data are coded as \code{NA} or 
#' a state sequence object built with \link[TraMineR]{seqdef} function
#' @param var the list of columns containing the trajectories. 
#' Default is NULL, i.e. all the columns. 
#' @param with.complete a logical stating if complete trajectories 
#' should be included or not in the plot
#'
#' @param ... parameters to be passed to the \link[TraMineR]{seqfplot} function
#' 
#' @details
#' This plot function is based on the \link[TraMineR]{seqfplot} function. 
#' To see which arguments can be changed, see the \link[TraMineR]{seqfplot} 
#' help. In particular, the number of most frequent 
#' patterns to be plotted can be changed with the argument \code{idxs}. By 
#' default, the 10 most frequent patterns are plotted. 
#' 
#'
#' @author Kevin Emery
#'
#' @examples
#' # Plot the 10 most common patterns of missing data
#' 
#' seqmissfplot(gameadd, var=1:4)
#' 
#' # Plot the 10 most common patterns of missing data discarding 
#' # complete trajectories
#' 
#' seqmissfplot(gameadd, var=1:4, with.missing = FALSE)
#' 
#' # Plot only the 5 most common patterns of missing data discarding 
#' # complete trajectories
#' 
#' seqmissfplot(gameadd, var=1:4, with.missing = FALSE, idxs = 1:5)
#' 
#' @export
seqmissfplot <- function(data, var = NULL, with.complete = TRUE, ...) {
  data <- dataxtract(data, var)
  
  if(inherits(data, "stslist")) {
    seqdata <- data
  }else{
    seqdata <- suppressMessages(seqdef(data, right=NA))
  }
  if (with.complete == TRUE) {
    seqmiss <- seqdata
  } else {
    seqmiss <- seqwithmiss(seqdata)
  }
  misspatterns <- matrix(NA, nrow(seqmiss), ncol(seqmiss))
  misspatterns <- as.data.frame(misspatterns)
  colnames(misspatterns) <- colnames(seqmiss)

  misspatterns <- matrix(NA, nrow(seqmiss), ncol(seqmiss))
  misspatterns <- as.data.frame(misspatterns)
  colnames(misspatterns) <- colnames(seqmiss)
  if (!is.na(attr(seqmiss, "nr"))) {
    misspatterns[seqmiss == attr(seqmiss, "nr")] <- "missing"
    misspatterns[seqmiss != attr(seqmiss, "nr")] <- "observed"
    seqtest <- suppressMessages(TraMineR::seqdef(misspatterns, 
        alphabet = c("observed", "missing"), cpal = c("blue", "red"), 
        xtstep = attr(seqmiss, "xtstep")))
    
    TraMineR::seqfplot(seqtest, ...)
  } else {
    misspatterns[is.na(seqmiss)] <- "missing"
    misspatterns[!is.na(seqmiss)] <- "observed"
    seqtest <- suppressMessages(TraMineR::seqdef(misspatterns, 
        alphabet = c("observed", "missing"), cpal = c("blue", "red"), 
        xtstep = attr(seqmiss, "xtstep")))
    
    TraMineR::seqfplot(seqtest, ...)
  }
}


#' Plot all the patterns of missing data.
#' 
#' #' @description
#' Plot function that renders all the patterns of missing data. This
#' function is based on the \link[TraMineR]{seqIplot} function. 
#'
#' @param data a data.frame where missing data are coded as \code{NA} or 
#' a state sequence object built with \link[TraMineR]{seqdef} function
#' @param var the list of columns containing the trajectories. 
#' Default is NULL, i.e. all the columns. 
#' @param with.complete a logical stating if complete trajectories 
#' should be included or not in the plot
#' @param ... parameters to be passed to the \link[TraMineR]{seqIplot}function
#' 
#' @examples
#' # Plot all the patterns of missing data
#' 
#' seqmissIplot(gameadd, var=1:4)
#' 
#' # Plot all the patterns of missing data discarding 
#' # complete trajectories
#' 
#' seqmissIplot(gameadd, var=1:4, with.missing = FALSE)
#' @author Kevin Emery
#'
#' @export
seqmissIplot <- function(data, var = NULL, with.complete = TRUE, ...) {
  data <- dataxtract(data, var)
  
  if(inherits(data, "stslist")) {
    seqdata <- data
  }else{
    seqdata <-suppressMessages(seqdef(data, right=NA))
  }
  if (with.complete == TRUE) {
    seqmiss <- seqdata
  } else {
    seqmiss <- seqwithmiss(seqdata)
  }
  
  misspatterns <- matrix(NA, nrow(seqmiss), ncol(seqmiss))
  misspatterns <- as.data.frame(misspatterns)
  colnames(misspatterns) <- colnames(seqmiss)
  if (!is.na(attr(seqmiss, "nr"))) {
    misspatterns[seqmiss == attr(seqmiss, "nr")] <- "missing"
    misspatterns[seqmiss != attr(seqmiss, "nr")] <- "observed"
    seqtest <- suppressMessages(TraMineR::seqdef(misspatterns, 
      alphabet = c("observed", "missing"), cpal = c("blue", "red"), 
      xtstep = attr(seqmiss, "xtstep")))

    TraMineR::seqIplot(seqtest, ...)
  } else {
    misspatterns[is.na(seqmiss)] <- "missing"
    misspatterns[!is.na(seqmiss)] <- "observed"
    seqtest <- suppressMessages(TraMineR::seqdef(misspatterns, 
      alphabet = c("observed", "missing"), cpal = c("blue", "red"), 
      xtstep = attr(seqmiss, "xtstep")))

    TraMineR::seqIplot(seqtest, ...)
  }
}

#' Identification and visualization of states that best characterize sequences
#' with missing data
#' 
#' @description Function based on the \link[TraMineRextras]{seqimplic}. 
#' Identification and visualization of the states that best characterize the 
#' sequence with missing data vs. the sequences without missing data at each 
#' position (time point). See the \link[TraMineRextras]{seqimplic} help 
#' for more details on how it works.
#'
#' @param data a data frame where missing data are coded as \code{NA} or 
#' a state sequence object built with \link[TraMineR]{seqdef} function
#' @param var the list of columns containing the trajectories. 
#' Default is NULL, i.e. all the columns.
#' @param ... parameters to be passed to the \link[TraMineRextras]{seqimplic}
#' function
#' 
#' @return returns a \code{seqimplic} object that can be plotted and printed. 
#' 
#' @examples
#' 
#' # For illustration purpose, we simulate missing data on the mvad dataset,
#' # available in the TraMineR package. The state "joblessness" state has a 
#' # higher probability of triggering a missing gap
#' 
#' \dontrun{
#' data(mvad, package = "TraMineR")
#' mvad.miss <- seqaddNA(mvad, var = 17:86, states.high = "joblessness")
#' 
#' # The states that best characterize sequences with missing data
#' implic <- seqmissimplic(mvad.miss, var = 17:86)
#' 
#' # Visualization of the results
#' plot(implic)
#' }
#' 
#' @author Kevin Emery
#'
#' @export
seqmissimplic <- function(data, var = NULL, ...){
  
  data <- dataxtract(data, var)
  
  if(inherits(data, "stslist")) {
    seqdata <- data
  }else{
    seqdata <-suppressMessages(seqdef(data, right=NA))
  }
  tt <- rep("missing", nrow(seqdata))
  tt[rowSums(seqdata == attr(seqdata, "nr")) == 0] <- "observed"
  imp <- suppressWarnings(TraMineRextras::seqimplic(seqdata, tt, na.rm=FALSE))
  return(imp)
}