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
#' seqmissfplot(gameadd, var = 1:4)
#'
#' # Plot the 10 most common patterns of missing data discarding
#' # complete trajectories
#'
#' seqmissfplot(gameadd, var = 1:4, with.missing = FALSE)
#'
#' # Plot only the 5 most common patterns of missing data discarding
#' # complete trajectories
#'
#' seqmissfplot(gameadd, var = 1:4, with.missing = FALSE, idxs = 1:5)
#'
#' @export
seqmissfplot <- function(data, var = NULL, with.complete = TRUE,
                         void.miss = TRUE, ...) {
  if (inherits(data, "stslist")) {
    seqmiss <- data
  } else {
    data <- dataxtract(data, var)
    seqmiss <- suppressMessages(seqdef(data, right = NA))
  }
  
  misspatterns <- matrix(NA, nrow(seqmiss), ncol(seqmiss))
  misspatterns <- as.data.frame(misspatterns)
  colnames(misspatterns) <- colnames(seqmiss)
  
  if (void.miss == FALSE) {
    misspatterns[seqmiss == attr(seqmiss, "nr")] <- "missing"
    misspatterns[seqmiss != attr(seqmiss, "nr")] <- "not missing"
    
    if (with.complete == FALSE) {
      misspatterns <- misspatterns[rowSums(misspatterns == "missing") != 0, ]
    }
    seqtest <- suppressMessages(TraMineR::seqdef(misspatterns,
                                                 alphabet = c("not missing", "missing"), cpal = c("#4F91CD", "#E34A33"),
                                                 xtstep = attr(seqmiss, "xtstep")
    ))
    
    TraMineR::seqfplot(seqtest, ...)
  } else {
    misspatterns[seqmiss == attr(seqmiss, "nr") |
                   seqmiss == attr(seqmiss, "void")] <- "missing"
    misspatterns[seqmiss != attr(seqmiss, "nr") &
                   seqmiss != attr(seqmiss, "void")] <- "not missing"
    
    if (with.complete == FALSE) {
      misspatterns <- misspatterns[rowSums(misspatterns == "missing") != 0, ]
    }
    seqtest <- suppressMessages(TraMineR::seqdef(misspatterns,
                                                 alphabet = c("not missing", "missing"), cpal = c("#4F91CD", "#E34A33"),
                                                 xtstep = attr(seqmiss, "xtstep")
    ))
    
    TraMineR::seqfplot(seqtest, ...)
  }
}
