#' Transform an object of class \code{seqimp} into a dataframe or a \code{mids}
#' object
#'
#' @description
#' The function converts a \code{seqimp} object into a specified format.
#' 
#'
#' @param data An object of class seqimp as created by the 
#' function \link[seqimpute]{seqimpute}
#' @param format The format in which the seqimp object should be returned. It 
#' could be: \code{"long"}, \code{"stacked"} and \code{"mids"}.
#' See the Details section for the interpretation.
#' @param include logical that indicates if the original dataset with missing
#' value should be included or not. This parameter does not apply 
#' if \code{format="mids"}.
#' @details
#' The argument \code{format} specifies the object that should be returned
#' by the function. It can take the following values
#' \describe{
#'  \item{\code{"long"}}{
#'  
#'  produces a data set in which imputed data sets are stacked vertically. 
#'  The following columns are added: 1) \code{.imp} referring to the 
#'  imputation number, and 2) \code{.id} the row names of the original dataset}
#'  \item{\code{"stacked"}}{
#'  
#'  the same as \code{"long"}, but without the inclusion of 
#'  the two columns \code{.imp} and \code{.id}}
#'  \item{\code{"mids"}}{
#'  
#'  produces an object of class \code{mids}, which is the format
#'  used by the \code{mice} package.}
#' }
#' 
#' @return Transform a \code{seqimp} object into the desired format.
#' 
#' @examples
#' 
#' \dontrun{
#' # Imputation with the MICT algorithm
#' imp <- seqimpute(data = gameadd, var = 1:4)
#' 
#' # The object imp is transformed to a dataframe, where completed datasets are
#' # stacked vertically
#' 
#' imp.stacked <- fromseqimp(data = imp, 
#'     format = "stacked", include = FALSE)
#' }
#'
#' @author Kevin Emery
#'
#' @export
fromseqimp <- function(data, format="long", include=FALSE)
{
  if(!inherits(data,"seqimp")){
    stop("data is not a seqimp object built with the seqimpute() function")
  }
  m <- data$m
  if(format=="long"){
    new_data <- bind_rows(data$imp)
    if(include==TRUE){
      new_data <- rbind(data$data, new_data)
      idx <- 0L:m
    }else{
      idx <- 1L:m
    }
    new_data <- data.frame(.imp = rep(idx, each = nrow(data$data)), 
        .id = rep.int(1L:nrow(data$data), length(idx)), new_data)
    
  }else if(format=="stacked"){
    new_data <- bind_rows(data$imp)
  }else if(format=="mids"){
    new_data <- bind_rows(data$imp)
    
    new_data <- rbind(data$data, new_data)
    idx <- 0L:m
    
    new_data <- data.frame(.imp = rep(idx, each = nrow(data$data)), 
        .id = rep.int(1L:nrow(data$data), length(idx)), new_data)
    new_data <- as.mids(new_data)
  }else{
    stop("format should be equal to either long,
    stacked, or mids'")
  }
  return(new_data)
}

#' Extract all the trajectories without missing value.
#'
#' @param data either a data frame containing sequences of a multinomial 
#' variable with missing data (coded as \code{NA}) or a state sequence 
#' object built with the TraMineR package
#' @param var the list of columns containing the trajectories. 
#' Default is NULL, i.e. all the columns.
#' @return Returns either a data frame or a state sequence object, depending
#' the type of data that was provided to the function
#' 
#' @examples
#' 
#' # Game addiction dataset
#' data(gameadd)
#' # Extract the trajectories without any missing data
#' gameadd.complete <- seqcomplete(gameadd, var = 1:4)
#' 
#' 
#' 
#' @author Kevin Emery
#'
#' @export
seqcomplete <- function(data, var = NULL) {
  data <- dataxtract(data, var)
  if (inherits(data, "stslist")) {
    rowsNA <- rowSums(data == attr(data, "nr"))
  }else{
    rowsNA <- rowSums(is.na(data))
  }
  return(data[rowsNA == 0, ])
}

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
  if(!inherits(impdata, "seqimp")) {
    stop("impdata is not a seqimp object")
  }
  nrows.o <- nrow(impdata$data)
  if(inherits(clustering,"data.frame") | inherits(clustering,"matrix")){
    if((nrow(clustering)!=nrow(impdata$data)|ncol(clustering)!=impdata$m)){
      stop("clustering should have the same number of rows as the dataset
           that was imputed and a number of columns that corresponds to the 
           number of multiple imputations")
    }
  }else{
    if(length(clustering)!=(impdata$m*nrow(impdata$data))){
      stop("if clustering is provided as a vector, it should have the length
           of the number of the number of imputation done times the number
           of rows of the dataset that was imputed")
    }else{      
      clustering <- matrix(clustering,nrow=nrow(impdata$data),ncol=impdata$m)
    }
    
  }
  
  cluster.data <- rep(NA,nrow(impdata$data))
  for(i in 1:nrow(impdata$data)){
    if(all(clustering[i,] == clustering[i,1])){
      cluster.data[i] <- clustering[i,1]
    }
  }
  impdata$data <- cbind(impdata$data,cluster.data)
  colnames(impdata$data)[ncol(impdata$data)] <- "cluster"
  
  for(i in 1:impdata$m){
    impdata$imp[[i]] <- cbind(impdata$imp[[i]],clustering[,i])
    colnames(impdata$imp[[i]])[ncol(impdata$imp[[i]])] <- "cluster"
  }
  
  
  return(impdata)
}

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