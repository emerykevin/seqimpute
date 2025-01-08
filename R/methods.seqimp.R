#' Print a \code{seqimp} object
#' @param x Object of class \code{seqimp}
#' @param ... additional arguments passed to other functions
#' 
#' @author Kevin Emery
#' @export
print.seqimp <- function(x, ...){
  cat("Class: seqimp\n")
  cat("Number of multiple imputations: ",x$m,"\n")
  cat("Method: ",x$method,"\n")
}

#' Summary of a \code{seqimp} object
#' @param object Object of class \code{seqimp}
#' @param ... additional arguments passed to other functions
#' 
#' @author Kevin Emery
#' @export
summary.seqimp <- function(object, ...){
  print(object, ...)
  invisible(object)
}


#' Plot a \code{seqimp} object
#' @description
#' Plot a \code{seqimp} object. The state distribution plot of the first 
#' \code{m} completed datasets is shown, possibly alongside the original 
#' dataset with missing data
#' 
#' @param x Object of class \code{seqimp}
#' @param m Number of completed datasets to show
#' @param include logical that indicates if the original dataset with missing
#' value should be plotted or not
#' @param ... Arguments to be passed to the seqdplot function
#' 
#' @author Kevin Emery
#' @export
plot.seqimp <- function(x, m = 5, include = TRUE, ...){
  tmp <- fromseqimp(x, format="long",include=include)
  tmp <- tmp[tmp$.imp<=m,-c(2)]
  if("cluster"%in%colnames(tmp)){
    tmp <- tmp[,!colnames(tmp)%in%c("cluster")]
  }
  seqtmp <-  suppressMessages(TraMineR::seqdef(tmp[,-c(1)]))
  suppressMessages(TraMineR::seqdplot(seqtmp, group=tmp$.imp,...))
}