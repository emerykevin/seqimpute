#' Plot the most common patterns of missing data.
#' 
#' @description
#' This function plots the most frequent patterns of missing data, based on the 
#' \link[TraMineR]{seqfplot} function.
#' 
#' @param data Either a data frame containing sequences of a categorical 
#' variable, where missing data are coded as \code{NA}, or a state sequence 
#' object created using the \link[TraMineR]{seqdef} function. 
#' 
#' @param var A vector specifying the columns of the dataset 
#' that contain the trajectories. Default is \code{NULL}, meaning all columns 
#' are used.
#' 
#' @param with.complete Logical, if \code{TRUE}, complete trajectories 
#' will be included in the plot.
#' 
#' @param void.miss Logical, if \code{TRUE}, treats void elements as 
#' missing values. Applies only to state sequence objects created with 
#' \link[TraMineR]{seqdef}. Note that the default behavior of \code{seqdef} 
#' is to treat missing data at the end of sequences as void elements.
#' 
#' @param ... Additional parameters passed to the \link[TraMineR]{seqfplot} 
#' function.
#' @details
#' This plot function is based on the \link[TraMineR]{seqfplot} function, 
#' allowing users to visualize patterns of missing data within sequences. 
#' For details on additional customizable arguments, see the 
#' \link[TraMineR]{seqfplot} documentation. 
#' 
#' By default, this function plots the 10 most frequent patterns. The number 
#' of patterns to be plotted can be adjusted using the \code{idxs} argument 
#' in \code{seqfplot}.
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
#' 
#' @export
seqmissfplot <- function(data, var = NULL, with.complete = TRUE, 
                         void.miss = TRUE, ...) {
  if(inherits(data, "stslist")) {
    seqmiss <- data
  }else{
    data <- dataxtract(data, var)
    seqmiss <- suppressMessages(seqdef(data, right=NA))
  }
  
  misspatterns <- matrix(NA, nrow(seqmiss), ncol(seqmiss))
  misspatterns <- as.data.frame(misspatterns)
  colnames(misspatterns) <- colnames(seqmiss)
  
  if(void.miss==FALSE){
    misspatterns[seqmiss == attr(seqmiss, "nr")] <- "missing"
    misspatterns[seqmiss != attr(seqmiss, "nr")] <- "not missing"
    
    if (with.complete == FALSE) {
      misspatterns <- misspatterns[rowSums(misspatterns=="missing")!=0,]
    }
    seqtest <- suppressMessages(TraMineR::seqdef(misspatterns, 
                                                 alphabet = c("not missing", "missing"), cpal = c("blue", "red"), 
                                                 xtstep = attr(seqmiss, "xtstep")))
    
    TraMineR::seqfplot(seqtest, ...)
    
  }else{
    misspatterns[seqmiss == attr(seqmiss, "nr") 
                 | seqmiss == attr(seqmiss, "void")] <- "missing"
    misspatterns[seqmiss != attr(seqmiss, "nr") 
                 & seqmiss != attr(seqmiss, "void")] <- "not missing"
    
    if (with.complete == FALSE) {
      misspatterns <- misspatterns[rowSums(misspatterns=="missing")!=0,]
    }
    seqtest <- suppressMessages(TraMineR::seqdef(misspatterns, 
                                                 alphabet = c("not missing", "missing"), cpal = c("blue", "red"), 
                                                 xtstep = attr(seqmiss, "xtstep")))
    
    TraMineR::seqfplot(seqtest, ...)
  }
}


#' Plot all the patterns of missing data.
#' 
#' @description
#' This function plots all patterns of missing data within sequences, based on 
#' the \link[TraMineR]{seqIplot} function.
#' 
#' @param data Either a data frame containing sequences of a categorical 
#' variable, where missing data are coded as \code{NA}, or a state sequence 
#' object created using the \link[TraMineR]{seqdef} function. 
#' 
#' @param var A vector specifying the columns of the dataset 
#' that contain the trajectories. Default is \code{NULL}, meaning all columns 
#' are used.
#' 
#' @param with.complete Logical, if \code{TRUE}, complete trajectories 
#' will be included in the plot.
#' 
#' @param void.miss Logical, if \code{TRUE}, treats void elements as 
#' missing values. Applies only to state sequence objects created with 
#' \link[TraMineR]{seqdef}. Note that the default behavior of \code{seqdef} 
#' is to treat missing data at the end of sequences as void elements.
#' 
#' @param ... Additional parameters passed to the \link[TraMineR]{seqIplot} function.
#' 
#' @details
#' This function uses \link[TraMineR]{seqIplot} to visualize all patterns of missing 
#' data within sequences. For further customization options, refer to the 
#' \link[TraMineR]{seqIplot} documentation.
#' 
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
seqmissIplot <- function(data, var = NULL, with.complete = TRUE, 
                         void.miss = TRUE, ...) {
  if(inherits(data, "stslist")) {
    seqmiss <- data
  }else{
    data <- dataxtract(data, var)
    seqmiss <- suppressMessages(seqdef(data, right=NA))
  }
  
  misspatterns <- matrix(NA, nrow(seqmiss), ncol(seqmiss))
  misspatterns <- as.data.frame(misspatterns)
  colnames(misspatterns) <- colnames(seqmiss)
  
  if(void.miss==FALSE){
    misspatterns[seqmiss == attr(seqmiss, "nr")] <- "missing"
    misspatterns[seqmiss != attr(seqmiss, "nr")] <- "not missing"
    
    if (with.complete == FALSE) {
      misspatterns <- misspatterns[rowSums(misspatterns=="missing")!=0,]
    }
    seqtest <- suppressMessages(TraMineR::seqdef(misspatterns, 
                                                 alphabet = c("not missing", "missing"), cpal = c("blue", "red"), 
                                                 xtstep = attr(seqmiss, "xtstep")))
    
    TraMineR::seqIplot(seqtest, ...)
    
  }else{
    misspatterns[seqmiss == attr(seqmiss, "nr") 
                 | seqmiss == attr(seqmiss, "void")] <- "missing"
    misspatterns[seqmiss != attr(seqmiss, "nr") 
                 & seqmiss != attr(seqmiss, "void")] <- "not missing"
    
    if (with.complete == FALSE) {
      misspatterns <- misspatterns[rowSums(misspatterns=="missing")!=0,]
    }
    seqtest <- suppressMessages(TraMineR::seqdef(misspatterns, 
                                                 alphabet = c("not missing", "missing"), cpal = c("blue", "red"), 
                                                 xtstep = attr(seqmiss, "xtstep")))
    
    TraMineR::seqIplot(seqtest, ...)
  }
}

#' Identification and visualization of states that best characterize sequences
#' with missing data
#' 
#' @description
#' This function identifies and visualizes states that best characterize 
#' sequences with missing data at each position (time point), comparing them to 
#' sequences without missing data at each position (time point). It is based on 
#' the \link[TraMineRextras]{seqimplic} function. For more information on the 
#' methodology, see the \code{seqimplic} documentation.
#' 
#' @param data Either a data frame containing sequences of a categorical 
#' variable, where missing data are coded as \code{NA}, or a state sequence 
#' object created using the \link[TraMineR]{seqdef} function. 
#' 
#' @param var A vector specifying the columns of the dataset 
#' that contain the trajectories. Default is \code{NULL}, meaning all columns 
#' are used.
#' 
#' @param void.miss Logical, if \code{TRUE}, treats void elements as 
#' missing values. This argument applies only to state sequence objects created 
#' with \link[TraMineR]{seqdef}. Note that the default behavior of \code{seqdef} 
#' is to treat missing data at the end of sequences as void elements.
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
seqmissimplic <- function(data, var = NULL, void.miss = TRUE,...){
  
  if(inherits(data, "stslist")) {
    seqdata <- data
  }else{
    data <- dataxtract(data, var)
    seqdata <-suppressMessages(seqdef(data, right=NA))
  }
  tt <- rep("missing", nrow(seqdata))
  
  if(void.miss==FALSE){
    tt[rowSums(seqdata == attr(seqdata, "nr")) == 0] <- "not missing"
  }else{
    tt[rowSums(seqdata == attr(seqdata, "nr")) == 0 & 
         rowSums(seqdata == attr(seqdata, "void")) == 0] <- "not missing"
  }
  imp <- suppressWarnings(TraMineRextras::seqimplic(seqdata, tt, na.rm=FALSE))
  return(imp)
}