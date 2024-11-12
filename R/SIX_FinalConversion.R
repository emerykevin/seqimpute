FinalResultConvert <- function(RESULT, ODClass, ODlevels,
  rownamesDataset, nrowsDataset, nr, nc, rowsNA, mi)
{
  RESULT <- apply(RESULT, 2, as.numeric)

  RESULT <- as.data.frame(RESULT)

    RESULT <- as.data.frame(
      sapply(RESULT, mapvalues,
      from = as.character(as.vector(1:length(ODlevels))),
      to = ODlevels, warn_missing = FALSE)
      )
    
    if (ODClass == "factor") {
      RESULT <- mutate_if(RESULT, is.character, as.factor)
    }



  #### We put again the rows having only NA's discarded at the beginning
  if (length(rowsNA) > 0) {
    if(nrowsDataset%in%rowsNA){
      RESULT <- rbind(RESULT,rep(NA, ncol(RESULT)))
      rownames(RESULT)[nrow(RESULT)] <- nrowsDataset
      if(length(rowsNA)>1){
        for (i in 1:(length(rowsNA)-1)) {
          if (rowsNA[i] == 1) {
            RESULT <- rbind(rep(NA, ncol(RESULT)), RESULT)
          } else {
            RESULT <- rbind(RESULT[1:(rowsNA[i] - 1), ],
                            rep(NA, ncol(RESULT)),
                            RESULT[rowsNA[i]:nrow(RESULT), ])
          }
        }
      }
    }else{
      for (i in 1:length(rowsNA)) {
        if (rowsNA[i] == 1) {
          RESULT <- rbind(rep(NA, ncol(RESULT)), RESULT)
        } else if (rowsNA[i] == nrowsDataset) {
          RESULT <- rbind(RESULT,rep(NA, ncol(RESULT)))
        } else {
          RESULT <- rbind(RESULT[1:(rowsNA[i] - 1), ],
                          rep(NA, ncol(RESULT)),
                          RESULT[rowsNA[i]:nrow(RESULT), ])
        }
      }
    }
  }
    

  return(RESULT)
}
