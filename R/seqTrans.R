#' Spotting impossible transitions in longitudinal categorical data
#'
#' The purpose of \code{seqTrans} is to spot impossible transitions
#' in longitudinal categorical data.
#'
#' @param data a data frame containing sequences of a multinomial 
#' variable with missing data (coded as \code{NA})
#' @param var the list of columns containing the trajectories. 
#' Default is NULL, i.e. all the columns. 
#' @param trans \code{character} vector gathering the impossible transitions. 
#' For example: trans <- c("1->3","1->4","2->1","4->1","4->3")
#'
#' @author Andre Berchtold and Kevin Emery
#'
#' @return It returns a matrix where each row is the position of an 
#' impossible transition.
#'
#' @examples
#' data(gameadd)
#'
#' seqTransList <- seqTrans(data = gameadd, var = 1:4, trans = c("yes->no"))
#'
#' @importFrom stringr str_count
#' @importFrom stringr str_detect
#' @importFrom stringr str_locate
#' @importFrom stringr str_locate_all
#'
#' @importFrom graphics plot
#'
#' @importFrom stats as.formula
#' @importFrom stats cutree
#' @importFrom stats lm
#' @importFrom stats predict
#' @importFrom stats rnorm
#' @importFrom stats runif
#'
#' @importFrom utils capture.output
#'
#' @importFrom Amelia missmap
#'
#' @importFrom TraMineR seqdef
#' @importFrom TraMineR seqfplot
#' @importFrom TraMineR seqdplot
#' @importFrom TraMineR seqsubm
#' @importFrom TraMineR seqdist
#'
#' @importFrom cluster agnes
#'
#' @importFrom plyr mapvalues
#'
#' @importFrom dfidx dfidx
#'
#' @importFrom rms lrm
#'
#' @importFrom mice as.mids
#'
#' @importFrom mlr makeClassifTask
#'
#' @importFrom ranger ranger
#'
#' @importFrom stats model.matrix
#'
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#'
#' @importFrom dplyr n_distinct
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate_if
#'
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' 
#' @importFrom parallelly availableCores
#'
#' @export
seqTrans <- function(data, var = NULL, trans) {
  
  data <- dataxtract(data, var)
  
  
  impTrans <- trans
  nr <- nrow(data)
  nc <- ncol(data)


 
  dataClass <- class(data[1, 1])
  if ((dataClass != "factor") & (dataClass != "character") 
      & (dataClass != "numeric")) {
    stop("/!\\ The class of the variables contained in your original dataset
         should be either 'factor', 'character', or 'numeric'")
  }


  if (dataClass == "factor"|dataClass=="character") {
    k <- length(sort(unique(as.vector(as.matrix(data)))))
  } else {
    k <- max(data)
  }

  # rows with NA
  i <- 1
  numbOfNAFilledLines <- 0
  while (i <= nrow(data)) {
    if (all(is.na(data[i, ]))) {
      data <- data[-i, ]
      numbOfNAFilledLines <- numbOfNAFilledLines + 1
    }
    i <- i + 1
  }
  # Updating the number of rows in data
  nr <- nrow(data)


  for (i in 1:length(impTrans)) {
    if (!str_detect(impTrans[i], "->")) {
      stop("/!\\ Warning, you should construct your transition(s) vector trans 
        with little arrows as 
        follow: trans <- c('...->...', '...->...', etc.).")
    }
    # Testing if we are effectively analyzing a transition or not
    locDash <- str_locate(impTrans[i], "-")
    firstState <- substr(impTrans[i], 1, locDash - 1)
    locSpike <- str_locate(impTrans[i], ">")
    secondState <- substr(impTrans[i], locSpike + 1, nchar(impTrans[1]))
 
  }









  # 3. Spotting impossible transitions ----------------------------


  ## Setup
  #
  # Transforming every line of data into class character and adding a dash 
  # inbetween every state
  dashes <- replicate(nc, "->")
  CharAndDashes <- function(data) {
    mytestChar <- paste(as.vector(rbind(data, dashes)), collapse = "")
  }
  dataCharAndDashes <- apply(data, 1, CharAndDashes)






  # 3.1 Number of listed impossible transitions -------------------------------
  # Identifying the patterns of the impossible transitions in each line of data
  countImpTrans <- function(dataCharAndDashes) {
    str_count(dataCharAndDashes, pattern = paste0("(?=", impTrans, ")"))
  }
  numbOfImpTransByRow <- sapply(dataCharAndDashes, countImpTrans)
  if (length(impTrans) > 1) {
    numbOfImpTrans <- rowSums(numbOfImpTransByRow)
  } else {
    numbOfImpTrans <- sum(numbOfImpTransByRow)
  }
  # Interrupting the program in case no impossible transitions among impTrans 
  # have been found
  if (sum(numbOfImpTrans) == 0) {
    message("Your dataset has no impossible transitions!")

    seqTransList <- matrix(0, 0, 2)
    colnames(seqTransList) <- c("row", "col")

    return(seqTransList)
  }

  if (length(impTrans) == 1) {
    numbOfImpTransByRow <- t(as.matrix(numbOfImpTransByRow))
  }





  # we will have "3 14" in place i,j
  startLocMat <- matrix(NA, nrow = length(impTrans), ncol = nrow(data))
  for (i in 1:length(impTrans)) {
    # We need this roundabout with "(?=)" because in the case of a transition
    # yes->yes->yes for example, if the user is interested in the transition
    # yes->yes, the simple version without "(?=)" won't spot
    # the second transtion because a part of it is used in the first one.
    Tmplist <- str_locate_all(dataCharAndDashes, paste0("(?=",impTrans[i],")"))
    for (j in 1:length(dataCharAndDashes)) {
      tempMat <- Tmplist[[j]]
      if (nrow(tempMat) > 0) {
        if (nrow(tempMat) > 1) {
          tempMat <- paste(tempMat[, 1], collapse = " ")
        }
        startLocMat[i, j] <- as.matrix(tempMat[1])
      }
    }
  }



  # Function to extract the different digits when a string is composed
  # of multiple digits
  Numextract <- function(string) {
    unlist(regmatches(string, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", string)))
  }
  # Function that translates the value corresponding to where the transition
  # occures in the string to where it occures in a observation
  # To do so, we extract the substring until the desired a transition occures
  # in the string and we compute the number +1 of ">" in this substring.
  # For example, if the observation j is of the form "yes->yes->no->yes->yes"
  # and the value 15 is stored in StartLocMat, we extract the substring
  # "yes->yes->no->", we compute 3 + 1. The target transition therefore occurs
  # in position 4.
  GetRealColPosition <- function(dataCharAndDashes, RowPos, ColPos) {
    dataCharAndDashesSubstring <- substr(dataCharAndDashes[RowPos], 1, ColPos)
    realColPosition <- str_count(dataCharAndDashesSubstring, pattern = ">") + 1
    return(realColPosition)
  }
  list_rows <- list()
  list_cols <- list()

  for (i in 1:length(impTrans)) {
    if (numbOfImpTrans[i] > 0) {
      TmpRows <- c()
      TmpCols <- c()
      for (j in 1:ncol(startLocMat)) {
        if (numbOfImpTransByRow[i, j] > 0) {
          detectTemp <- str_detect(startLocMat[i, j], " ")
          if (detectTemp == FALSE) {
            TmpRows <- c(TmpRows, j)
            TmpCols <- c(TmpCols, 
              GetRealColPosition(dataCharAndDashes = dataCharAndDashes, 
              RowPos = j, ColPos = as.numeric(startLocMat[i, j])))
          } else {
            hiddenCol <- str_count(startLocMat[i, j], " ")
            TmpRows <- c(TmpRows, rep(j, hiddenCol + 1))
            extractedNumbers <- Numextract(startLocMat[i, j])
            TmpVec <- c()
            for (k in 1:length(extractedNumbers)) {
              TmpVec[k] <- GetRealColPosition(
                dataCharAndDashes=dataCharAndDashes, 
                RowPos = j, ColPos = as.numeric(extractedNumbers[k]))
            }
            TmpCols <- c(TmpCols, TmpVec)
          }
        }
      }
      list_rows[[i]] <- TmpRows
      list_cols[[i]] <- TmpCols
    }
  }


  MaxNumTransitions <- max(sapply(list_cols, function(x) length(x)))

  rowMat <- matrix(NA, length(impTrans), MaxNumTransitions)
  colMat <- matrix(NA, length(impTrans), MaxNumTransitions)
  #
  # # Converting into a dataframe
  # rowMat <- as.data.frame(rowMat)
  # # Renaming the columns of rowMat
  # colnames_rowMat <- paste(1:ncol(rowMat),")",sep='')
  # rownames(rowMat) <- impTrans
  # colnames(rowMat) <- colnames_rowMat
  #
  # # Converting into a dataframe
  # colMat <- as.data.frame(colMat)
  # # Renaming the columns of colMat
  # colnames_colMat <- paste(1:ncol(colMat),")",sep='')
  # rownames(colMat) <- impTrans
  # colnames(colMat) <- colnames_colMat
  #
  for (i in 1:length(impTrans)) {
    if (numbOfImpTrans[i] > 0) {
      rowMat[i, 1:length(list_rows[[i]])] <- list_rows[[i]]
      colMat[i, 1:length(list_cols[[i]])] <- list_cols[[i]]
    }
  }

  if (length(impTrans) > 1) {
    rowMat <- rowMat[rowSums(!is.na(rowMat)) > 0, ]
    colMat <- colMat[rowSums(!is.na(colMat)) > 0, ]
  }



  # 3.3 Summary data frame --------------------------------------------------
  #
  # Update of impTrans with an arrow
  impTransOverview <- data.frame(c(impTrans, "", "Total:"), 
    c(numbOfImpTrans, "", sum(numbOfImpTrans)))
  colnames(impTransOverview) <- c("Transitions", "Occurence")

  if (length(rowMat) > 0) {
    seqTransList <- matrix(0, length(rowMat), 2)
    seqTransList[, 1] <- rowMat
    seqTransList[, 2] <- colMat
  }
  colnames(seqTransList) <- c("row", "col")


  return(seqTransList)
}
