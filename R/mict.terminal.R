mict.terminal <- function(data, imp, MaxTermGapSize, REFORDT_L, pastDistrib, regr, 
                          npt, nco, ncot, nr, nc, ud, available, k, noise, ...)
{
  OD <- data$OD
  covariates <- data$CO
  time.covariates <- data$COt
  COtsample <- data$COtsample

  if (npt > nc - MaxTermGapSize) {
    npt <- nc - MaxTermGapSize
  }
  
  train <- compute.traindata(data, MaxGap=0, order=0, shift=0, np=npt, nc, nr, nf=0, k, 
                          pastDistrib, futureDistrib=FALSE, col=0, 
                          frame.radius=0, regr, timing=FALSE)
  reglog <- fitmodel(train, regr, ...)
  for (order in 1:MaxTermGapSize){
    imp <- CreatedModelImputation(order=order, covariates, train, 
                                  time.covariates, COtsample, OD, imp, pastDistrib, futureDistrib=FALSE, available, 
                                  REFORDT_L, ncot, nc, np=npt, nf=0, k, regr, reglog, noise, 
                                  shift=0, MaxGap=MaxTermGapSize)
    
  }
  
  return(imp)
}
