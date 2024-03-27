#' Summary of the types of gaps among a dataset
#'
#' The \code{seqQuickLook()} function aimed at providing an overview of the
#' number and size of the different types of gaps 
#' spread in the original dataset.
#'
#' @param data a data.frame where missing data are coded as NA or 
#' a state sequence object built with \link[TraMineR]{seqdef} function
#' @param var the list of columns containing the trajectories. 
#' Default is NULL, i.e. all the columns. 
#' @param np number of previous observations in the imputation model of the 
#' internal gaps.
#' @param nf number of future observations in the imputation model of the 
#' internal gaps.
#' 
#' @details The distinction between internal and SLG gaps depends on the 
#' number of previous (\code{np}) and future (\code{nf}) observations that are 
#' set for the \code{MICT} and \code{MICT-timing} algorithms.
#' @author Andre Berchtold and Kevin Emery
#'
#' @return Returns a  \code{data.frame} object that summarizes, for each 
#' type of gaps (Internal Gaps, Initial Gaps, Terminal Gaps, 
#' LEFT-hand side SLG, RIGHT-hand side SLG, Both-hand side SLG),
#' the minimum length, the maximum length, the total number of gaps and
#' the total number of missing they contain.
#'
#' @examples
#' data(gameadd)
#'
#' seqQuickLook(data = gameadd, var = 1:4, np = 1, nf = 1)
#'
#' @export


seqQuickLook <- function(data, var=NULL, np = 1, nf = 1) {
  
  data <- dataxtract(data, var)

  if (inherits(data, "stslist")) {
    valuesNA <- attr(data, "nr")
    data <- data.frame(data)
    data[data == valuesNA] <- NA
  }

  
  MDGapsChart <- matrix(0, 6, 4)
  MDGapsChart <- as.data.frame(MDGapsChart)
  rownames(MDGapsChart) <- c("Internal Gaps", "Initial Gaps", "Terminal Gaps", 
    "LEFT-hand side SLG", "RIGHT-hand side SLG", "BOTH-hand side SLG")
  colnames(MDGapsChart) <- c("MinGapSize","MaxGapSize","numbOfGaps","sumNAGaps")
  
  if(sum(is.na(data))==0){
    return(MDGapsChart)
  }
  
  # 1. Initial tests on parameters --------------------------------------------


  # 1.1 Testing the class of the variables of the original dataset data ------
  dataClass <- class(data[1, 1])
  if ((dataClass != "factor") & (dataClass != "character") & 
      (dataClass != "numeric")) {
    stop("/!\\ The class of the variables contained in your original dataset
         should be either 'factor','character' or 'numeric'")
  }


  if (dataClass == "factor" | dataClass == "character") {
    k <- length(sort(unique(as.vector(as.matrix(data)))))
  } else {
    k <- max(data)
  }

  datadata <- list()
  datadata[c("data", "CO", "COt", "rowsNA")] <- deleteNaRows(data, 
    matrix(NA, nrow = 1, ncol = 1), matrix(NA, nrow = 1, ncol = 1))



  # Naming the number of rows and columns of the original dataset data
  nr <- nrow(datadata$data)
  nc <- ncol(datadata$data)





  # 1. Analysis of data and creation of matrices ORDER, ORDER2 and ORDER3 -----
  datadata[c("MaxInitGapSize", "InitGapSize", "MaxTermGapSize", "TermGapSize", 
    "MaxGap", "ORDER", "ORDER2", "ORDER3")] <- OrderCreation(datadata$data, 
          nr, nc)


  # 2. Computation of the order of imputation of each MD 
  # if (max(datadata$ORDER) != 0) {
  #   datadata[c("ORDERSLGLeft", "ORDERSLGRight", "ORDERSLGBoth", "LongGap", 
  #     "MaxGap","REFORD_L","ORDER")] <- ImputeOrderComputation(datadata$ORDER, 
  #         datadata$ORDER3, datadata$MaxGap, np, nf, nr, nc)
  # } else {
  #   datadata$ORDERSLGLeft <- matrix(nrow = datadata$nr, ncol = datadata$nc,0)
  #   datadata$ORDERSLGRight <- matrix(nrow =datadata$nr, ncol = datadata$nc,0)
  #   datadata$ORDERSLGBoth <- matrix(nrow = datadata$nr, ncol = datadata$nc,0)
  #   datadata$LongGap <- FALSE
  # }
  
  
  datadata[c("ORDERSLGLeft", "ORDERSLGRight", "ORDERSLGBoth", "LongGap", 
      "MaxGap", "REFORD_L", "ORDER")] <- ImputeOrderComputation(datadata$ORDER, 
      datadata$ORDER3, datadata$MaxGap, np, nf, nr, nc)

  if (max(datadata$ORDER) > 0) {
    MDGapsChart[1, 1] <- min(datadata$ORDER[datadata$ORDER > 0])
    MDGapsChart[1, 2] <- max(datadata$ORDER)
  }
  # Number of Internal Gaps
  # Transforming datadata$ORDER in a single row vector
  datadata$ORDERVect <- as.vector(t(datadata$ORDER))
  # Transforming this single row vector into class character
  datadata$ORDERVectChar <- paste(datadata$ORDERVect, collapse = "")
  # Identifying the patterns "0 1" (this is the signature  of an internal gap
  # (it always indicates the beginning of an internal gap!))
  MDGapsChart[1, 3] <- str_count(datadata$ORDERVectChar, pattern = "01")
  MDGapsChart[1, 4] <- sum(datadata$ORDER > 0)


  if (max(datadata$InitGapSize) > 0) {
    MDGapsChart[2,1] <- min(datadata$InitGapSize[datadata$InitGapSize > 0])
    MDGapsChart[2,2] <- max(datadata$InitGapSize)
    MDGapsChart[2,3]<-sum(table(datadata$InitGapSize[datadata$InitGapSize>0]))
    MDGapsChart[2,4] <- sum(datadata$InitGapSize)
  }

  if (max(datadata$TermGapSize) > 0) {
    MDGapsChart[3,1] <- min(datadata$TermGapSize[datadata$TermGapSize > 0])
    MDGapsChart[3,2] <- max(datadata$TermGapSize)
    MDGapsChart[3,3]<-sum(table(datadata$TermGapSize[datadata$TermGapSize>0]))
    MDGapsChart[3,4] <- sum(datadata$TermGapSize)
  }

  LeftGap <- rowSums(datadata$ORDERSLGLeft > 0)
  if (max(LeftGap) > 0) {
    MDGapsChart[4, 1] <- min(LeftGap[LeftGap > 0])
    MDGapsChart[4, 2] <- max(LeftGap)
    MDGapsChart[4, 3] <- sum(table(LeftGap[LeftGap > 0]))
    MDGapsChart[4, 4] <- sum(LeftGap)
  }

  RightGap <- rowSums(datadata$ORDERSLGRight > 0)
  if (max(RightGap) > 0) {
    MDGapsChart[5, 1] <- min(RightGap[RightGap > 0])
    MDGapsChart[5, 2] <- max(RightGap)
    MDGapsChart[5, 3] <- sum(table(RightGap[RightGap > 0]))
    MDGapsChart[5, 4] <- sum(RightGap)
  }

  BothGap <- rowSums(datadata$ORDERSLGBoth > 0)
  if (max(BothGap) > 0) {
    MDGapsChart[6, 1] <- min(BothGap[BothGap > 0])
    MDGapsChart[6, 2] <- max(BothGap)
    MDGapsChart[6, 3] <- sum(table(BothGap[BothGap > 0]))
    MDGapsChart[6, 4] <- sum(BothGap)
  }

  return(MDGapsChart)
}
