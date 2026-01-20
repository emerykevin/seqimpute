
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
#' @param ... Additional parameters passed to the
#' \link[TraMineR]{seqIplot} function.
#'
#' @details
#' This function uses \link[TraMineR]{seqIplot} to visualize all patterns of
#' missing data within sequences. For further customization options, refer to
#' the \link[TraMineR]{seqIplot} documentation.
#'
#'
#' @examples
#' # Plot all the patterns of missing data
#'
#' seqmissIplot(gameadd, var = 1:4)
#'
#' # Plot all the patterns of missing data discarding
#' # complete trajectories
#'
#' seqmissIplot(gameadd, var = 1:4, with.missing = FALSE)
#' @author Kevin Emery
#'
#' @export
seqmissIplot <- function(data, var = NULL, with.complete = TRUE,
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
    
    TraMineR::seqIplot(seqtest, ...)
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
    
    TraMineR::seqIplot(seqtest, ...)
  }
}
