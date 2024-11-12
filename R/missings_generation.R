#' Generation of missing on longitudinal categorical data.
#' 
#' @description Generation of missing data under the form of gaps, which
#' is the typical form of missing data with longitudinal data.
#' It simulates MCAR or MAR missing data.
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
#' @param propdata Proportion of observations for which missing data 
#' is simulated, as a decimal between 0 and 1.
#' 
#' @param maxgap Maximum length of a missing data gap.
#' 
#' @param maxprop Maximum proportion of missing data allowed in a sequence, 
#' as a decimal between 0 and 1. If the proportion exceeds this value, the 
#' simulation is rerun for the sequence.
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
#' \dontrun{
#' data(mvad, package = "TraMineR")
#' mvad.miss <- seqaddNA(mvad, var = 17:86)
#' 
#' 
#' # Generate missing data on mvad where joblessness is more likely to trigger 
#' # a missing data gap
#' mvad.miss2 <- seqaddNA(mvad, var = 17:86,  states.high = "joblessness")
#' }
#'
#' @export
seqaddNA <- function(data, var = NULL, states.high = NULL, propdata = 1, 
  pstart.high = 0.1, pstart.low = 0.005, maxgap = 3, maxprop=0.75, 
  only.traj = FALSE)
{
  data.traj <- dataxtract(data, var)
  
  sizehalf <- round(propdata * nrow(data.traj))
  rowsmiss <- sample(1:nrow(data.traj), size = sizehalf, replace = FALSE)
  matrix.missing <- matrix(1, nrow(data.traj), ncol(data.traj))
  for (i in 1:length(rowsmiss)) {
    nmis <- ncol(data.traj)
    length.gap <- 0
    while (nmis > floor(maxprop * ncol(data.traj))) {
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
                prob = c(66, 34))
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
