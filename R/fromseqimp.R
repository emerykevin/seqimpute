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
