#' Generation of missing on longitudinal categorical data.
#' 
#' @description Generation of missing data under the form of gaps, which
#' is the typical form of missing data with longitudinal data.
#' It simulates MCAR or MAR missing data.
#'
#' @param data a data frame containing sequences of a multinomial 
#' variable with missing data (coded as \code{NA})
#' @param var the list of columns containing the trajectories. 
#' Default is \code{NULL}, i.e. all the columns. 
#' @param states.high list of states that have a larger probability of 
#' triggering a subsequent missing data gap
#' @param pstart.high probability to start a missing data for the 
#' states specified with the \code{states.high} argument
#' @param pstart.low probability to start a missing data for the 
#' other states
#' @param propdata proportion  observations for which missing data is simulated
#' @param maxgap maximum length of a missing data gap
#' 
#' @param only.traj logical that specifies whether only the trajectories should 
#' be returned (\code{only.traj=TRUE}), or 
#' the whole data (\code{only.traj=FALSE})
#'
#' @return Returns a data frame on which missing data were simulated
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
  pstart.high = 0.1, pstart.low = 0.005, maxgap = 3, only.traj = FALSE)
{
  data.traj <- dataxtract(data, var)
  
  sizehalf <- round(propdata * nrow(data.traj))
  rowsmiss <- sample(1:nrow(data.traj), size = sizehalf, replace = FALSE)
  matrix.missing <- matrix(1, nrow(data.traj), ncol(data.traj))
  for (i in 1:length(rowsmiss)) {
    nmis <- ncol(data.traj)
    length.gap <- 0
    while (nmis > floor(0.75 * ncol(data.traj))) {
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
