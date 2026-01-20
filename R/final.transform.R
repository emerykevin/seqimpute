final.transform <- function(
    RESULT, data, ODlevels,
    rownamesDataset, nrowsDataset, nr, nc, rowsNA, mi) {
  
  RESULT <- apply(RESULT, 2, as.numeric)
  
  if (length(rowsNA) > 0) {
    RESULT_NEW <- matrix(NA, nrow = nrowsDataset, ncol = ncol(RESULT))
    rownames(RESULT_NEW) <- 1:nrowsDataset
    
    non_NA_rows <- setdiff(1:nrowsDataset, rowsNA)
    RESULT_NEW[non_NA_rows, ] <- RESULT[non_NA_rows, ]
    
    RESULT <- RESULT_NEW
  }
  
  
  RESULT <- as.data.frame(RESULT)
  
  RESULT <- as.data.frame(
    sapply(RESULT, mapvalues,
           from = as.character(as.vector(1:length(ODlevels))),
           to = ODlevels, warn_missing = FALSE
    )
  )
  
  if (inherits(data[1, 1], "factor")) {
    RESULT <- mutate_if(RESULT, is.character, factor, levels = ODlevels)
  }
  
  return(RESULT)
}