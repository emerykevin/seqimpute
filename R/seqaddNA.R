#' Generation of missing on longitudinal categorical data.
#'
#' @description Generation of missing data in sequence based on a Markovian
#' approach.
#'
#' @details
#'
#' The first time point of a trajectory has a \code{pstart.low} probability to
#' be missing. For the next time points, the probability to be missing depends
#' on the previous time point. There are four cases:
#'
#' 1. If the previous time point is missing and the maximum length of a
#' missing gap, which is specified by the argument \code{maxgap}, is reached,
#' the time point is set as observed.
#'
#' 2. If the previous time point is missing, but the maximum length of a gap is
#' not reached, there is a \code{pcont} probability that this time point is missing.
#'
#' 3. If the previous time point is observed and the previous time point belongs
#' to the list of states specified by \code{pstart.high}, the probability to
#' be missing is \code{pstart.high}.
#'
#' 4. If the previous time point is observed but the previous time point does not
#' belong to the list of states specified by \code{pstart.high}, the
#' probability to be missing is \code{pstart.low}.
#'
#' If the proportion of missing data in a given trajectory exceeds the
#' proportion specified by \code{maxprop}, the missing data simulation is
#' repeated for the sequence.
#'
#' @param data A data frame containing sequences of a categorical (multinomial)
#' variable, where missing data are coded as \code{NA}.
#'
#' @param var A vector specifying the columns of the dataset
#' that contain the trajectories. Default is \code{NULL}, meaning all columns
#' are used.
#'
#' @param states.high A list of states with a higher probability of
#' initiating a subsequent missing data gap.
#'
#' @param pstart.high Probability of starting a missing data gap for the
#' states specified in the \code{states.high} argument.
#'
#' @param pstart.low Probability of starting a missing data gap for all
#' other states.
#'
#' @param pcont Probability of a missing data gap to continue.
#'
#' @param propdata Proportion of trajectories for which missing data
#' is simulated, as a decimal between 0 and 1.
#'
#' @param maxgap Maximum length of a missing data gap.
#'
#' @param maxprop Maximum proportion of missing data allowed in a sequence,
#' as a decimal between 0 and 1.
#'
#' @param only.traj Logical, if \code{TRUE}, only the
#' trajectories (specified in \code{var}) are returned. If \code{FALSE},
#' the entire data frame is returned.
#'
#' @return A data frame with simulated missing data.
#'
#'
#' @author Kevin Emery
#'
#' @examples
#' # Generate MCAR missing data on the mvad dataset
#' # from the TraMineR package
#'
#' \donttest{
#' data(mvad, package = "TraMineR")
#' mvad.miss <- seqaddNA(mvad, var = 17:86)
#'
#'
#' # Generate missing data on mvad where joblessness is more likely to trigger
#' # a missing data gap
#' mvad.miss2 <- seqaddNA(mvad, var = 17:86, states.high = "joblessness")
#' }
#'
#' @export
seqaddNA <- function(data, var = NULL, states.high = NULL, pstart.high = 0.1, 
                        pstart.low = 0.005, pcont=0.66, propdata=1, maxgap = 3, 
                        maxprop=0.75, only.traj = FALSE)
{
  data.traj <- data
  
  sizehalf <- round(propdata * nrow(data.traj))
  rowsmiss <- sample(1:nrow(data.traj), size = sizehalf, replace = FALSE)
  matrix.missing <- matrix(1, nrow(data.traj), ncol(data.traj))
  for (i in 1:length(rowsmiss)) {
    nmis <- ncol(data.traj)
    while (nmis > floor(maxprop * ncol(data.traj))) {
      length.gap <- 0
      for (j in 1:ncol(data.traj)) {
        if (length.gap == maxgap) {
          matrix.missing[rowsmiss[i], j] <- 1
          length.gap <- 0
        } else {
          if (j == 1) {
            matrix.missing[rowsmiss[i], j] <- sample(x = c(0, 1), size = 1, 
                                                     prob = c(pstart.low, 1 - pstart.low))
          } else {
            if (matrix.missing[rowsmiss[i], j - 1] == 1) {
              if (data.traj[rowsmiss[i], j - 1] %in% states.high) {
                matrix.missing[rowsmiss[i], j] <- sample(x = c(0, 1), size = 1, 
                                                         prob = c(pstart.high, 1 - pstart.high))
              } else {
                matrix.missing[rowsmiss[i], j] <- sample(x = c(0, 1), size = 1, 
                                                         prob = c(pstart.low, 1 - pstart.low))
              }
            } else {
              matrix.missing[rowsmiss[i], j] <- sample(x = c(0, 1), size = 1, 
                                                       prob = c(pcont, 1-pcont))
            }
          }
          if (matrix.missing[rowsmiss[i], j] == 0) {
            length.gap <- length.gap + 1
          } else {
            length.gap <- 0
          }
        }
      }
      nmis <- sum(matrix.missing[rowsmiss[i], ] == 0)
    }
  }
  
  data.traj[matrix.missing == 0] <- NA
  
  if(is.null(var) | only.traj == TRUE){
    data <- data.traj
  }else{
    data <- dataputback(data, var=var, data.traj=data.traj)
  }
  return(data)
}
