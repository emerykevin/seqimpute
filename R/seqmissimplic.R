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
#' with \link[TraMineR]{seqdef}. Note that the default behavior of
#' \link[TraMineR]{seqdef} is to treat missing data at the end of sequences as
#' void elements.
#' @param ... parameters to be passed to the \link[TraMineRextras]{seqimplic}
#' function
#'
#' @return returns a \link[TraMineRextras]{seqimplic} object that can be 
#' plotted and printed.
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
seqmissimplic <- function(data, var = NULL, void.miss = TRUE, ...) {
  if (inherits(data, "stslist")) {
    seqdata <- data
  } else {
    data <- dataxtract(data, var)
    seqdata <- suppressMessages(seqdef(data, right = NA))
  }
  tt <- rep("missing", nrow(seqdata))
  
  if (void.miss == FALSE) {
    tt[rowSums(seqdata == attr(seqdata, "nr")) == 0] <- "not missing"
  } else {
    tt[rowSums(seqdata == attr(seqdata, "nr")) == 0 &
         rowSums(seqdata == attr(seqdata, "void")) == 0] <- "not missing"
  }
  imp <- suppressWarnings(TraMineRextras::seqimplic(seqdata, tt, na.rm = FALSE))
  return(imp)
}

