#' Function that adds the clustering result to a \code{seqimp} object
#' obtained with the \code{seqimpute} function
#'
#' @param impdata An object of class \code{seqimp} as created by the
#' \code{seqimpute} function
#' @param clustering clustering made on the multiple imputed dataset. Can
#' either be a dataframe or a matrix, where each row correspond to an
#' observation and each column to a multiple imputed dataset
#'
#' @return Returns a \code{seqimp} object containing the cluster to which each
#' sequence in each imputed dataset belongs. Specifically, a column named
#' cluster is added to the imputed datasets.
#'
#' @export
addcluster <- function(impdata, clustering) {
  if (!inherits(impdata, "seqimp")) {
    stop("impdata is not a seqimp object")
  }
  nrows.o <- nrow(impdata$data)
  if (inherits(clustering, "data.frame") | inherits(clustering, "matrix")) {
    if ((nrow(clustering) != nrow(impdata$data) | ncol(clustering) != impdata$m)) {
      stop("clustering should have the same number of rows as the dataset
           that was imputed and a number of columns that corresponds to the
           number of multiple imputations")
    }
  } else {
    if (length(clustering) != (impdata$m * nrow(impdata$data))) {
      stop("if clustering is provided as a vector, it should have the length
           of the number of the number of imputation done times the number
           of rows of the dataset that was imputed")
    } else {
      clustering <- matrix(clustering, nrow = nrow(impdata$data), ncol = impdata$m)
    }
  }

  cluster.data <- rep(NA, nrow(impdata$data))
  for (i in 1:nrow(impdata$data)) {
    if (all(clustering[i, ] == clustering[i, 1])) {
      cluster.data[i] <- clustering[i, 1]
    }
  }
  impdata$data <- cbind(impdata$data, cluster.data)
  colnames(impdata$data)[ncol(impdata$data)] <- "cluster"

  for (i in 1:impdata$m) {
    impdata$imp[[i]] <- cbind(impdata$imp[[i]], clustering[, i])
    colnames(impdata$imp[[i]])[ncol(impdata$imp[[i]])] <- "cluster"
  }


  return(impdata)
}
