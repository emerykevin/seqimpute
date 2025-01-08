#' Extract all the trajectories with at least one missing value
#'
#' @param data either a data frame containing sequences of a multinomial 
#' variable with missing data (coded as \code{NA}) or a state sequence 
#' object built with the TraMineR package
#' @param var the list of columns containing the trajectories. 
#' Default is NULL, i.e. all the columns.
#' @return Returns either a data frame or a state sequence object, 
#' depending the type of data that was provided to the function
#' 
#' @examples
#' 
#' # Game addiction dataset
#' data(gameadd)
#' # Extract the trajectories without any missing data
#' gameadd.withmiss <- seqwithmiss(gameadd, var = 1:4)
#' 
#' @author Kevin Emery
#'
#' @export
seqwithmiss <- function(data, var = NULL) {
  data <- dataxtract(data, var)
  if (inherits(data, "stslist")){
    # deal with the case where the user used NA as the internal code 
    # for missing data
    if (is.na(attr(data, "nr"))){
      tmp <- data
      for (i in 1:ncol(data)) {
        tmp[, i] <- as.character(data[, i])
      }
      rowsNA <- rowSums(is.na(tmp))
    } else {
      rowsNA <- rowSums(data == attr(data, "nr"))
    }
  }else{
    rowsNA <- rowSums(is.na(data))
  }
  
  return(data[rowsNA != 0, ])
}