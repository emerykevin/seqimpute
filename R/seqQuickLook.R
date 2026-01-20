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


seqQuickLook <- function(data, var = NULL, np = 1, nf = 1) {
  data <- dataxtract(data, var)
  
  if (inherits(data, "stslist")) {
    valuesNA <- attr(data, "nr")
    data <- data.frame(data)
    data[data == valuesNA] <- NA
  }
  
  
  MDGapsChart <- matrix(0, 6, 4)
  MDGapsChart <- as.data.frame(MDGapsChart)
  rownames(MDGapsChart) <- c(
    "Internal Gaps", "Initial Gaps", "Terminal Gaps",
    "LEFT-hand side SLG", "RIGHT-hand side SLG", "BOTH-hand side SLG"
  )
  colnames(MDGapsChart) <- c("MinGapSize", "MaxGapSize", "numbOfGaps", "sumNAGaps")
  
  if (sum(is.na(data)) == 0) {
    return(MDGapsChart)
  }
  
  if (!inherits(data[1, 1], c("factor", "character", "numeric"))){
    stop("/!\\ The class of the variables contained in your original dataset
         should be either 'factor', 'character', or 'numeric'")
  }
  
  if (inherits(data[1, 1], c("factor", "character"))) {
    k <- length(sort(unique(as.vector(as.matrix(data)))))
  } else {
    k <- max(data)
  }
  
  rowsNA <- rows.allmiss(data)
  
  nr <- nrow(data)
  nc <- ncol(data)
  
  imporder <- compute.order(data, nr, nc, np, nf, 1, 1,
                            end.impute = TRUE
  )
  
  if (imporder$maxInternal > 0) {
    numGapsInternal <- compute.num.gaps(imporder$internal, imporder$maxInternal)
    MDGapsChart[1, 1] <- min(which(numGapsInternal > 0))
    MDGapsChart[1, 2] <- imporder$maxInternal
    MDGapsChart[1, 3] <- sum(numGapsInternal)
    MDGapsChart[1, 4] <- sum(numGapsInternal * (1:imporder$maxInternal))
  } else {
    MDGapsChart[1, ] <- 0
  }
  
  if (imporder$maxInitial > 0) {
    numGapsInitial <- compute.num.gaps(imporder$initial, imporder$maxInitial)
    MDGapsChart[2, 1] <- min(which(numGapsInitial > 0))
    MDGapsChart[2, 2] <- imporder$maxInitial
    MDGapsChart[2, 3] <- sum(numGapsInitial)
    MDGapsChart[2, 4] <- sum(numGapsInitial * (1:imporder$maxInitial))
  } else {
    MDGapsChart[2, ] <- 0
  }
  
  if (imporder$maxTerminal > 0) {
    numGapsTerminal <- compute.num.gaps(imporder$terminal, imporder$maxTerminal)
    MDGapsChart[3, 1] <- min(which(numGapsTerminal > 0))
    MDGapsChart[3, 2] <- imporder$maxTerminal
    MDGapsChart[3, 3] <- sum(numGapsTerminal)
    MDGapsChart[3, 4] <- sum(numGapsTerminal * (1:imporder$maxTerminal))
  } else {
    MDGapsChart[3, ] <- 0
  }
  
  if (max(imporder$maxLeftSLG) > 0) {
    numGapsSLG <- matrix(0, np, max(imporder$maxLeftSLG))
    for (h in 2:np) {
      if (imporder$maxLeftSLG[h] > 0) {
        numGapsSLG[h, 1:imporder$maxLeftSLG[h]] <- compute.num.gaps(imporder$SLGleft[[h]], imporder$maxLeftSLG[h])
      }
    }
    MDGapsChart[4, 1] <- min(which(numGapsSLG > 0))
    MDGapsChart[4, 2] <- max(imporder$maxLeftSLG)
    MDGapsChart[4, 3] <- sum(numGapsSLG)
    MDGapsChart[4, 4] <- sum(t(numGapsSLG) * (1:max(imporder$maxLeftSLG)))
  } else {
    MDGapsChart[4, ] <- 0
  }
  
  if (max(imporder$maxRightSLG) > 0) {
    numGapsSLG <- matrix(0, nc - 1, max(imporder$maxRightSLG))
    for (h in (nc - 1):(nc - nf + 1)) {
      if (imporder$maxRightSLG[h] > 0) {
        numGapsSLG[h, 1:imporder$maxRightSLG[h]] <- compute.num.gaps(imporder$SLGright[[h]], imporder$maxRightSLG[h])
      }
    }
    
    MDGapsChart[5, 1] <- min(which(colSums(numGapsSLG) > 0))
    MDGapsChart[5, 2] <- max(imporder$maxRightSLG)
    MDGapsChart[5, 3] <- sum(numGapsSLG)
    MDGapsChart[5, 4] <- sum(t(numGapsSLG) * (1:max(imporder$maxRightSLG)))
  } else {
    MDGapsChart[5, ] <- 0
  }
  
  
  if (max(imporder$maxBothSLG) > 0) {
    minGap <- max(imporder$maxBothSLG)
    sumGaps <- 0
    sumNumGaps <- 0
    for (g in 2:np) {
      for (h in (nc - 1):(nc - nf + 1)) {
        if (imporder$maxBothSLG[g, h] > 0) {
          tmpGaps <- compute.num.gaps(imporder$SLGboth[[g]][[h]], imporder$maxBothSLG[g, h])
          
          minGap <- min(minGap, min(which(tmpGaps > 0)))
          sumNumGaps <- sumNumGaps + sum(tmpGaps)
          
          sumGaps <- sumGaps + sum(tmpGaps * (1:imporder$maxBothSLG[g, h]))
        }
      }
    }
    MDGapsChart[6, 1] <- minGap
    MDGapsChart[6, 2] <- max(imporder$maxBothSLG)
    MDGapsChart[6, 3] <- sumNumGaps
    MDGapsChart[6, 4] <- sumGaps
  } else {
    MDGapsChart[6, ] <- 0
  }
  
  
  return(MDGapsChart)
}


compute.num.gaps <- function(order, maxGap) {
  numGaps <- rep(0, maxGap)
  numGaps[maxGap] <- nrow(order[[1]])
  if (maxGap > 1) {
    for (i in 1:(maxGap - 1)) {
      numGaps[maxGap - i] <- nrow(order[[i + 1]]) - nrow(order[[i]])
    }
  }
  return(numGaps)
}