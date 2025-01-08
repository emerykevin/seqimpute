mict.initial <- function(data, imp, futureDistrib, REFORDI_L, MaxInitGapSize, nr, 
                        nc, ud, nco, ncot, nfi, regr, k, available, noise, ...) 
{
  OD <- data$OD
  covariates <- data$CO
  time.covariates <- data$COt
  COtsample <- data$COtsample

  if (nfi > nc - MaxInitGapSize) {
    nfi <- nc - MaxInitGapSize
    
  }
  
  train <- compute.traindata(data, MaxGap=0, order=0, shift=0, np=0, nc, nr, nfi, k, 
                  pastDistrib=FALSE, futureDistrib, col=0, frame.radius=0, regr,
                  timing=FALSE)
  
  
  reglog <- fitmodel(train, regr, ...)

  for (order in 1:MaxInitGapSize){
    imp <- CreatedModelImputation(order=order, covariates, train, 
            time.covariates, COtsample, OD, imp, pastDistrib=FALSE, futureDistrib, available, 
            REFORDI_L, ncot, nc, np=0, nfi, k, regr, reglog, noise, 
            shift=0, MaxGap=MaxInitGapSize)
    
  }
  return(imp)
}