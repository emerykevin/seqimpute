# 1. General testing ####
## a. past and future ####
gameadd.traj <- gameadd[,1:4]

test_that("MICT works without error", {
  expect_no_error(seqimpute(gameadd.traj))
})

test_that("MICT-timing works without error", {
  expect_no_error(seqimpute(gameadd.traj, timing = TRUE))
})

test_that("MICT single imputation works without error", {
  expect_no_error(seqimpute(gameadd.traj, m = 1))
})

test_that("MICT-timing multiple imputation works without error", {
  expect_no_error(seqimpute(gameadd.traj, timing = TRUE, m=1))
})

## a.bis more than 2 predictors ####
gameadd.traj.2 <- cbind(gameadd.traj,gameadd.traj)
colnames(gameadd.traj.2) <- c(1:8)

test_that("MICT works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.2, np=2, m=2))
})

test_that("MICT-timing works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.2, timing = TRUE, np=2, m=2))
})

test_that("MICT works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.2, nf=2, m=2))
})

test_that("MICT-timing works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.2, timing = TRUE, nf=2, m=2))
})


test_that("MICT works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.2, np=2, nf=2, m=2))
})

test_that("MICT-timing works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.2, timing = TRUE, np=2, nf=2, m=2))
})
## b. past only ####
test_that("MICT works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, nf=0, m=2))
})

test_that("MICT-timing works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, timing = TRUE, nf=0, m=2))
})

## c. future only ####
test_that("MICT works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, np=0, m=2))
})

test_that("MICT-timing works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, timing = TRUE, np=0, m=2))
})

## d. large nfi + npt####
test_that("MICT works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, nfi=3, m=2))
})

test_that("MICT-timing works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, timing = TRUE, nfi=3, m=2))
})

test_that("MICT works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, npt=3, m=2))
})

test_that("MICT-timing works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, timing = TRUE, npt=3, m=2))
})

# 2. returned object ####
test_that("seqimp object is return by MICT", {
  skip_on_cran()
  expect_s3_class(seqimpute(gameadd.traj, m=2), "seqimp")
})

test_that("seqimp object is return by MICT-timing", {
  skip_on_cran()
  expect_s3_class(seqimpute(gameadd.traj, timing=TRUE, m=2), "seqimp")
})

## d. no missing values ######
imp <- seqimpute(gameadd.traj, m=1)
test_that("MICT", {
  skip_on_cran()
  expect_message(seqimpute(imp$imp[[1]], m=2), 
                 "This dataset has no missing values!")
})

test_that("MICT-timing", {
  skip_on_cran()
  expect_message(seqimpute(imp$imp[[1]], timing=TRUE, m=2), 
                 "This dataset has no missing values!")
})


# 3. type of data ####

## a. numeric - error ####
gameadd.num <- gameadd[,1:4]
for(i in 1:ncol(gameadd.num)){
  levels(gameadd.num[,i]) <- c(1,2)
  gameadd.num[,i] <- as.numeric(gameadd.num[,i])
}

test_that("MICT -- numeric dataset throws error", {
  skip_on_cran()
  expect_error(seqimpute(gameadd.num, m=2))
})

test_that("MICT-t -- numeric dataset throws error", {
  skip_on_cran()
  expect_error(seqimpute(gameadd.num, m=2, timing=TRUE))
})

## b. factor ####
test_that("factor returns factor", {
  skip_on_cran()
  expect_s3_class(seqimpute(gameadd.traj, m=2)$imp[[1]][,1],"factor")
})

test_that("factor returns factor", {
  skip_on_cran()
  expect_s3_class(seqimpute(gameadd.traj, m=2,
      timing=TRUE)$imp[[1]][,1], "factor")
})

## c. character ####
gameadd.traj.ch <- gameadd.traj
for(i in 1:ncol(gameadd.traj.ch)){
  gameadd.traj.ch[,i] <- as.character(gameadd.traj.ch[,i])
}

test_that("character returns character", {
  skip_on_cran()
  expect_equal(inherits(seqimpute(gameadd.traj.ch, 
        m=2)$imp[[1]][,1],"character"),TRUE)
})

test_that("character returns character", {
  skip_on_cran()
  expect_equal(inherits(seqimpute(gameadd.traj.ch, 
        m=2, timing=TRUE)$imp[[1]][,1],"character"),TRUE)
})


# 4. parallel computing ####
## a. General working ####
test_that("MICT works with parallel computing", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, m=2, ParExec = TRUE))
})

test_that("MICT-timing works with parallel computing", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, timing=TRUE, m=2, ParExec = TRUE))
})

## b. ParExec with m=1 ####
test_that("MICT", {
  skip_on_cran()
  expect_message(seqimpute(gameadd.traj, ParExec = TRUE, m=1))
})

test_that("MICT-timing", {
  skip_on_cran()
  expect_message(seqimpute(gameadd.traj, ParExec = TRUE, m=1, timing=TRUE))
})

## c. ncores specified #####
test_that("MICT", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, ParExec = TRUE, ncores=2))
})

test_that("MICT-timing", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, timing=TRUE, ParExec = TRUE, ncores=2))
})

## d. setRNGSeed #####
### i. identical results with identical seed #####
test_that("MICT", {
  skip_on_cran()
  expect_equal(sum(seqimpute(gameadd.traj, ParExec = TRUE, ncores=2, SetRNGSeed=2)$imp[[1]] != seqimpute(gameadd.traj, ParExec = TRUE, ncores=2, SetRNGSeed = 2)$imp[[1]]),0)
})

test_that("MICT-timing", {
  skip_on_cran()
  expect_equal(sum(seqimpute(gameadd.traj, ParExec = TRUE, ncores=2, SetRNGSeed=2, timing=TRUE)$imp[[1]] != seqimpute(gameadd.traj, ParExec = TRUE, ncores=2, SetRNGSeed = 2, timing=TRUE)$imp[[1]]),0)
})

### ii. identical resuls with identiccal seed, single imputation 
test_that("MICT", {
  skip_on_cran()
  expect_equal(sum(seqimpute(gameadd.traj, ParExec = TRUE, ncores=2, SetRNGSeed=2, m=1)$imp[[1]] != seqimpute(gameadd.traj, ParExec = TRUE, ncores=2, SetRNGSeed = 2, m=1)$imp[[1]]),0)
})

test_that("MICT-timing", {
  skip_on_cran()
  expect_equal(sum(seqimpute(gameadd.traj, ParExec = TRUE, ncores=2, SetRNGSeed=2, timing=TRUE, m=1)$imp[[1]] != seqimpute(gameadd.traj, ParExec = TRUE, ncores=2,  SetRNGSeed = 2, timing=TRUE, m=1)$imp[[1]]),0)
})

# 5. var argument works correctly ####
test_that("MICT works", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd, var=1:4, m=2))
})

test_that("MICT-timing works", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd, var=1:4, timing = TRUE, m=2))
})

test_that("MICT single imputation works", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd, var=1:4, m = 1))
})

test_that("MICT-timing multiple imputation works", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd, var=1:4, timing = TRUE, m=1))
})

# 6. values imputed ####
## a. MICT #####
imp.MICT <- seqimpute(gameadd, var=1:4, m=2)

test_that("MICT -- m imputed datasets", {
  skip_on_cran()
  expect_equal(length(imp.MICT$imp), 2)
})


test_that("no missing in imputed 1", {
  skip_on_cran()
  expect_equal(sum(is.na(imp.MICT$imp[[1]])),0)
})

test_that("no missing in imputed 2", {
  skip_on_cran()
  expect_equal(sum(is.na(imp.MICT$imp[[2]])),0)
})

test_that("non missing values from the dataset not changed by the imputation", {
  skip_on_cran()
  expect_equal(imp.MICT$imp[[1]][!is.na(gameadd.traj)],
               gameadd.traj[!is.na(gameadd.traj)])
})

test_that("non missing values from the dataset not changed by the imputation", {
  skip_on_cran()
  expect_equal(imp.MICT$imp[[2]][!is.na(gameadd.traj)],
               gameadd.traj[!is.na(gameadd.traj)])
})

## b. MICTt #####
imp.MICTt <- seqimpute(gameadd, var=1:4, m=2)

test_that("MICTt -- m imputed datasets", {
  skip_on_cran()
  expect_equal(length(imp.MICTt$imp), 2)
})


test_that("no missing in imputed 1", {
  skip_on_cran()
  expect_equal(sum(is.na(imp.MICTt$imp[[1]])),0)
})

test_that("no missing in imputed 2", {
  skip_on_cran()
  expect_equal(sum(is.na(imp.MICTt$imp[[2]])),0)
})

test_that("non missing values from the dataset not changed by the imputation", {
  skip_on_cran()
  expect_equal(imp.MICTt$imp[[1]][!is.na(gameadd.traj)],
               gameadd.traj[!is.na(gameadd.traj)])
})

test_that("non missing values from the dataset not changed by the imputation", {
  skip_on_cran()
  expect_equal(imp.MICTt$imp[[2]][!is.na(gameadd.traj)],
               gameadd.traj[!is.na(gameadd.traj)])
})

# 7. SLG imputation ####
## A. multinom ####
### a. MICT ####
test_that("Left SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, np = 2,m=2))
})

test_that("Left SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, np = 2,m=2,available=FALSE))
})

test_that("Left SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, np = 2, m=2, pastDistrib=TRUE, futureDistrib=TRUE))
})

test_that("Right SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, nf = 2,m=2))
})

test_that("Right SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, nf = 2,m=2,available=FALSE))
})

test_that("Right SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, nf = 2,m=2,pastDistrib=TRUE, futureDistrib=TRUE))
})

test_that("Both SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, nf = 2, np = 2,m=2))
})

### b. MICT-timing ####
test_that("Left SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, np = 2, timing = TRUE,m=2))
})

test_that("Right SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, nf = 2, timing = TRUE,m=2))
})

test_that("Both SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, nf = 2, np = 2, timing = TRUE,m=2))
})

### B. rf ####
## a. MICT ####
test_that("Left SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, regr="rf", np = 2,m=2))
})

test_that("Right SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, regr="rf", nf = 2,m=2))
})

test_that("Both SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, regr="rf", nf = 2, np = 2,m=2))
})

## b. MICT-timing ####
test_that("Left SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, regr="rf", np = 2, timing = TRUE,m=2))
})

test_that("Right SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, regr="rf", nf = 2, timing = TRUE,m=2))
})

test_that("Both SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, regr="rf", nf = 2, np = 2, timing = TRUE,m=2))
})

## bis np=0 or nf=0####
## A. multinom ####
### a. MICT ####
test_that("Left SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, np = 2, nf=0,m=2))
})

test_that("Right SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, nf = 2, np=0, m=2))
})

### b. MICT-timing ####
test_that("Left SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, np = 2, nf=0, timing = TRUE,m=2))
})

test_that("Right SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, nf = 2, np=0, timing = TRUE,m=2))
})


### B. rf ####
## a. MICT ####
test_that("Left SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, regr="rf", np = 2, nf=0, m=2))
})

test_that("Right SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, regr="rf", nf = 2, np= 0, m=2))
})


## b. MICT-timing ####
test_that("Left SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, regr="rf", np = 2, nf=0, timing = TRUE,m=2))
})

test_that("Right SLG", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, regr="rf", nf = 2, np=0, timing = TRUE,m=2))
})
# 8. available ####
## a. past and future ####
gameadd.traj <- gameadd[,1:4]

test_that("MICT works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, available=FALSE, m=2))
})

test_that("MICT-timing works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, available=FALSE, timing = TRUE, m=2))
})

test_that("MICT single imputation works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, available=FALSE, m = 1))
})

test_that("MICT-timing multiple imputation works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, available=FALSE, timing = TRUE, m=1))
})

## b. past only ####
test_that("MICT works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, available=FALSE, nf=0, m=2))
})

test_that("MICT-timing works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, available=FALSE, timing = TRUE, nf=0, m=2))
})

## c. future only ####
test_that("MICT works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, available=FALSE, np=0, m=2))
})

test_that("MICT-timing works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, available=FALSE, timing = TRUE, np=0, m=2))
})


# 9. covariates + time-varying covariates ####
## a. past and future ####
test_that("MICT with covariates and time-varying covariates", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd, var=1:4, m=2, covariates=c("Gender","Age"),
        time.covariates = 8:11))
})

test_that("MICT with one covariate", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd, var=1:4, m=2, covariates=c("Gender")))
})
  
test_that("MICT with covariates and time-varying covariates", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd, var=1:4, m=2, 
                              time.covariates = 8:11))
})

test_that("MICT-timing with covariates and time-varying covariates", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd, var=1:4, m=2, covariates=c("Gender","Age"),
                            time.covariates = 8:11, timing=TRUE))
})

test_that("MICT-timing with one covariate", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd, var=1:4, m=2, covariates=c("Gender"), timing=T))
})
  
test_that("MICT-timing with one time-varying covariate", {
  skip_on_cran()
    expect_no_error(seqimpute(gameadd, var=1:4, m=2, 
                              time.covariates = 8:11, timing=T))
  })

## b. only past ####
test_that("MICT with covariates and time-varying covariates", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd, var=1:4, m=2, covariates=c("Gender","Age"),
                            time.covariates = 8:11, nf=0))
})

test_that("MICT with one covariate", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd, var=1:4, m=2, covariates=c("Gender"),nf=0))
})

test_that("MICT with covariates and time-varying covariates", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd, var=1:4, m=2, 
                            time.covariates = 8:11, nf=0))
})

test_that("MICT-timing with covariates and time-varying covariates", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd, var=1:4, m=2, covariates=c("Gender","Age"),
                            time.covariates = 8:11, timing=TRUE, nf=0))
})

test_that("MICT-timing with one covariate", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd, var=1:4, m=2, covariates=c("Gender"), timing=T, nf=0))
})

test_that("MICT-timing with one time-varying covariate", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd, var=1:4, m=2, 
                            time.covariates = 8:11, timing=T, nf=0))
})

## b. only future ####
test_that("MICT with covariates and time-varying covariates", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd, var=1:4, m=2, covariates=c("Gender","Age"),
                            time.covariates = 8:11, np=0))
})

test_that("MICT with one covariate", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd, var=1:4, m=2, covariates=c("Gender"),np=0))
})

test_that("MICT with covariates and time-varying covariates", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd, var=1:4, m=2, 
                            time.covariates = 8:11, np=0))
})

test_that("MICT-timing with covariates and time-varying covariates", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd, var=1:4, m=2, covariates=c("Gender","Age"),
                            time.covariates = 8:11, timing=TRUE, np=0))
})

test_that("MICT-timing with one covariate", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd, var=1:4, m=2, covariates=c("Gender"), timing=T, np=0))
})

test_that("MICT-timing with one time-varying covariate", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd, var=1:4, m=2, 
                            time.covariates = 8:11, timing=T, np=0))
})


# 10. rows with only NA ####
## A. multiple imputation #####
gameadd.row.NA.1 <- gameadd.traj
gameadd.row.NA.1[1,] <- NA

test_that("MICT with first row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.1,m=2))
})
  
test_that("MICT-t with first row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.1,m=2))
})

gameadd.row.NA.500 <- gameadd.traj
gameadd.row.NA.500[500,] <- NA

test_that("MICT with last row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.500,m=2))
})

test_that("MICT-t with last row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.500,m=2))
})

gameadd.row.NA.100 <- gameadd.traj
gameadd.row.NA.100[100,] <- NA

test_that("MICT with intermediate row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.100,m=2))
})

test_that("MICT-t with intermediate row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.100,m=2))
})

gameadd.row.NA.3r <- gameadd.traj
gameadd.row.NA.3r[c(1,100,500),] <- NA

test_that("MICT with three rows NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.3r,m=2))
})

test_that("MICT-t with three rows NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.3r,m=2))
})


gameadd.row.NA.twolast <- gameadd.traj
gameadd.row.NA.twolast[c(499,500),] <- NA

test_that("MICT with two last row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.twolast,m=2))
})

test_that("MICT-t with two last row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.twolast,m=2))
})


# b. with covariates #####
gameadd.row.NA.1 <- gameadd
gameadd.row.NA.1[1,1:4] <- NA

test_that("MICT with first row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.1,var=1:4,m=2, covariates=c("Gender","Age")))
})

test_that("MICT-t with first row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.1,var=1:4,m=2, covariates=c("Gender","Age")))
})

gameadd.row.NA.500 <- gameadd
gameadd.row.NA.500[500,1:4] <- NA

test_that("MICT with last row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.500,var=1:4,m=2, covariates=c("Gender","Age")))
})

test_that("MICT-t with last row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.500,var=1:4,m=2, covariates=c("Gender","Age")))
})

gameadd.row.NA.100 <- gameadd
gameadd.row.NA.100[100,1:4] <- NA

test_that("MICT with intermediate row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.100,var=1:4,m=2, covariates=c("Gender","Age")))
})

test_that("MICT-t with intermediate row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.100,var=1:4,m=2, covariates=c("Gender","Age")))
})

gameadd.row.NA.3r <- gameadd
gameadd.row.NA.3r[c(1,100,500),1:4] <- NA

test_that("MICT with three rows NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.3r,var=1:4,m=2, covariates=c("Gender","Age")))
})

test_that("MICT-t with three rows NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.3r,var=1:4,m=2, covariates=c("Gender","Age")))
})


gameadd.row.NA.twolast <- gameadd
gameadd.row.NA.twolast[c(499,500),1:4] <- NA

test_that("MICT with two last rows NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.twolast,var=1:4,m=2, covariates=c("Gender","Age")))
})

test_that("MICT-t with two last rows NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.twolast,var=1:4,m=2, covariates=c("Gender","Age")))
})


## B. single imputation #####
gameadd.row.NA.1 <- gameadd.traj
gameadd.row.NA.1[1,] <- NA

test_that("MICT with first row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.1,m=1))
})

test_that("MICT-t with first row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.1,m=1,timing=T))
})

gameadd.row.NA.500 <- gameadd.traj
gameadd.row.NA.500[500,] <- NA

test_that("MICT with last row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.500,m=1))
})

test_that("MICT-t with last row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.500,m=1))
})

gameadd.row.NA.100 <- gameadd.traj
gameadd.row.NA.100[100,] <- NA

test_that("MICT with intermediate row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.100,m=1))
})

test_that("MICT-t with intermediate row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.100,m=1))
})

gameadd.row.NA.3r <- gameadd.traj
gameadd.row.NA.3r[c(1,100,500),] <- NA

test_that("MICT with three rows NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.3r,m=1))
})

test_that("MICT-t with three rows NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.3r,m=1))
})


gameadd.row.NA.twolast <- gameadd.traj
gameadd.row.NA.twolast[c(499,500),] <- NA

test_that("MICT with two last rows NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.3r,m=1))
})

test_that("MICT-t with two last rows NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.3r,m=1))
})

# b. with covariates #####
gameadd.row.NA.1 <- gameadd
gameadd.row.NA.1[1,1:4] <- NA

test_that("MICT with first row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.1,var=1:4,m=1, covariates=c("Gender","Age")))
})

test_that("MICT-t with first row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.1,var=1:4,m=1, covariates=c("Gender","Age")))
})

gameadd.row.NA.500 <- gameadd
gameadd.row.NA.500[500,1:4] <- NA

test_that("MICT with last row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.500,var=1:4,m=1, covariates=c("Gender","Age")))
})

test_that("MICT-t with last row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.500,var=1:4,m=1, covariates=c("Gender","Age")))
})

gameadd.row.NA.100 <- gameadd
gameadd.row.NA.100[100,1:4] <- NA

test_that("MICT with intermediate row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.100,var=1:4,m=1, covariates=c("Gender","Age")))
})

test_that("MICT-t with intermediate row NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.100,var=1:4,m=1, covariates=c("Gender","Age")))
})

gameadd.row.NA.3r <- gameadd
gameadd.row.NA.3r[c(1,100,500),1:4] <- NA

test_that("MICT with three rows NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.3r,var=1:4,m=1, covariates=c("Gender","Age")))
})

test_that("MICT-t with three rows NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.3r,var=1:4,m=1, covariates=c("Gender","Age")))
})


gameadd.row.NA.twolast <- gameadd
gameadd.row.NA.twolast[c(499,500),1:4] <- NA

test_that("MICT with two last rows NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.3r,var=1:4,m=1, covariates=c("Gender","Age")))
})

test_that("MICT-t with two last rows NA", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.row.NA.3r,var=1:4,m=1, covariates=c("Gender","Age")))
})

# 10. pastDistrib #####

## a. past and future ####
gameadd.traj <- gameadd[,1:4]

test_that("MICT works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, pastDistrib=TRUE, m=2))
})

test_that("MICT-timing works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, pastDistrib=TRUE, timing = TRUE, m=2))
})

test_that("MICT single imputation works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, pastDistrib=TRUE, m = 1))
})

test_that("MICT-timing multiple imputation works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, pastDistrib=TRUE, timing = TRUE, m=1))
})

## b. past only ####
test_that("MICT works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, pastDistrib=TRUE, nf=0, m=2))
})

test_that("MICT-timing works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, pastDistrib=TRUE, timing = TRUE, nf=0, m=2))
})

## c. future only ####
test_that("MICT works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, pastDistrib=TRUE, np=0, m=2))
})

test_that("MICT-timing works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, pastDistrib=TRUE, timing = TRUE, np=0, m=2))
})


# 11. futureDistrib #####

## a. past and future ####
gameadd.traj <- gameadd[,1:4]

test_that("MICT works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, futureDistrib=TRUE, m=2))
})

test_that("MICT-timing works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, futureDistrib=TRUE, timing = TRUE, m=2))
})

test_that("MICT single imputation works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, futureDistrib=TRUE, m = 1))
})

test_that("MICT-timing multiple imputation works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, futureDistrib=TRUE, timing = TRUE, m=1))
})

## b. past only ####
test_that("MICT works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, futureDistrib=TRUE, nf=0, m=2))
})

test_that("MICT-timing works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, futureDistrib=TRUE, timing = TRUE, nf=0, m=2))
})

## c. future only ####
test_that("MICT works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, futureDistrib=TRUE, np=0, m=2))
})

test_that("MICT-timing works without error", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, futureDistrib=TRUE, timing = TRUE, np=0, m=2))
})


# 12. deprecated and deleted arguments ####
test_that("deprecated CO", {
  skip_on_cran()
  expect_warning(seqimpute(gameadd.traj, CO=4:8,m=2))
})

test_that("deleted mice.return", {
  skip_on_cran()
  expect_warning(seqimpute(gameadd.traj, mice.return=TRUE,m=2))
})

# 13. absurd choice of arguments ####
## a. wrong number of predictors ####
test_that("nfi<0", {
  skip_on_cran()
  expect_error(seqimpute(gameadd.traj, nfi=-1,m=2))
})

test_that("npt<0", {
  skip_on_cran()
  expect_error(seqimpute(gameadd.traj, npt=-1,m=2))
})

test_that("np<0", {
  skip_on_cran()
  expect_error(seqimpute(gameadd.traj, np=-1,m=2))
})

test_that("nf<0", {
  skip_on_cran()
  expect_error(seqimpute(gameadd.traj, nf=-1,m=2))
})

test_that("np=0 and nf=0", {
  skip_on_cran()
  expect_error(seqimpute(gameadd.traj, np=0, nf=0,m=2))
})

# b. length time.cov wrong ####
test_that("np=0 and nf=0", {
  skip_on_cran()
  expect_error(seqimpute(gameadd, var=1:4, time.covariates = c(8:11,9),m=2))
})

# c. regr lm ####
test_that("np=0 and nf=0", {
  skip_on_cran()
  expect_error(seqimpute(gameadd, var=1:4, regr="lm",m=2))
})

# 14. more than two categories ####
library("TraMineR")
data(mvad)
mvad.miss <- seqaddNA(mvad[,17:86])

# a. regr multinom####
test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(mvad.miss, regr="multinom",m=2))
})

test_that("MICTt - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(mvad.miss, regr="multinom", timing=TRUE,m=2))
})


# b. regr rf ####
test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(mvad.miss, regr="rf",m=2))
})

test_that("MICTt - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(mvad.miss, regr="rf", timing=TRUE,m=2))
})

# 15. num.trees ####
test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, regr="rf",m=1, num.trees=100))
})

# 16. one level ####
gameadd.traj.1lev <- gameadd.traj
levels(gameadd.traj.1lev[,3]) <- c("no","no")

test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.1lev, regr="multinom",m=1))
})

test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.1lev, regr="multinom",m=1, timing=T))
})

test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.1lev, regr="multinom",m=1, np=2))
})

test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.1lev, regr="multinom",m=1, timing=T, np=2))
})

test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.1lev, regr="multinom",m=1, nf=2))
})

test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.1lev, regr="multinom",m=1, timing=T, nf=2))
})

test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.1lev, regr="multinom",m=1, nf=2, np=2))
})

test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.1lev, regr="multinom",m=1, timing=T, nf=2, np=2))
})


gameadd.traj.1lev <- gameadd.traj
levels(gameadd.traj.1lev[,4]) <- c("no","no")

test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.1lev, regr="multinom",m=1))
})

test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.1lev, regr="multinom",m=1, timing=T))
})


gameadd.traj.1lev <- gameadd.traj
levels(gameadd.traj.1lev[,1]) <- c("no","no")

test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.1lev, regr="multinom",m=1))
})

test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.1lev, regr="multinom",m=1, timing=T))
})


gameadd.traj.1lev <- gameadd.traj
levels(gameadd.traj.1lev[,2]) <- c("no","no")

test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.1lev, regr="multinom",m=1, np=2))
})

test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.1lev, regr="multinom",m=1, timing=T, np=2))
})

test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.1lev, regr="multinom",m=1, nf=2))
})

test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.1lev, regr="multinom",m=1, timing=T, nf=2))
})

test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.1lev, regr="multinom",m=1, nf=2, np=2))
})

test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj.1lev, regr="multinom",m=1, timing=T, nf=2, np=2))
})

# 17. end.impute argument ####
test_that("MICT timing - end.impute - one imp", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, regr="multinom", m=1, timing=T, nf=1, np=1, end.impute=TRUE))
})

test_that("MICT timing - end.impute - one imp - sum miss", {
  skip_on_cran()
  expect_equal(sum(is.na(seqimpute(gameadd.traj, regr="multinom", m=1, timing=T, nf=1, np=1, end.impute=FALSE)$imp[[1]])),33)
})

imp <- seqimpute(gameadd.traj, regr="multinom", m=2, timing=T, nf=1, np=1, end.impute=FALSE)

test_that("MICT timing - end.impute - sum missings in imputed 1", {
  skip_on_cran()
  expect_equal(sum(is.na(imp$imp[[1]])),33)
})

test_that("MICT timing - end.impute - sum missings in imputed 2", {
  skip_on_cran()
  expect_equal(sum(is.na(imp$imp[[2]])),33)
})

test_that("MICT timing - end.impute - nf=2,np=2", {
  skip_on_cran()
  expect_equal(sum(is.na(seqimpute(gameadd.traj, regr="multinom", m=1, timing=T, nf=2, np=2, end.impute=FALSE)$imp[[1]])),33)
})



test_that("MICT - end.impute - one imp", {
  skip_on_cran()
  expect_no_error(seqimpute(gameadd.traj, regr="multinom", m=1, timing=F, nf=1, np=1, end.impute=TRUE))
})

test_that("MICT - end.impute - one imp - sum miss", {
  skip_on_cran()
  expect_equal(sum(is.na(seqimpute(gameadd.traj, regr="multinom", m=1, timing=F, nf=1, np=1, end.impute=FALSE)$imp[[1]])),33)
})

imp <- seqimpute(gameadd.traj, regr="multinom", m=2, timing=F, nf=1, np=1, end.impute=FALSE)

test_that("MICT - end.impute - sum missings in imputed 1", {
  skip_on_cran()
  expect_equal(sum(is.na(imp$imp[[1]])),33)
})

test_that("MICT - end.impute - sum missings in imputed 2", {
  skip_on_cran()
  expect_equal(sum(is.na(imp$imp[[2]])),33)
})

test_that("MICT - end.impute - nf=2,np=2", {
  skip_on_cran()
  expect_equal(sum(is.na(seqimpute(gameadd.traj, regr="multinom", m=1, timing=F, nf=2, np=2, end.impute=FALSE)$imp[[1]])),33)
})


# 18. sequence object ####
test_that("MICT - end.impute - one imp - sum miss", {
  skip_on_cran()
  library("TraMineR")
  seqgame <- seqdef(gameadd, var=1:4)
  expect_equal(sum(is.na(seqimpute(seqgame, regr="multinom", m=1, timing=F, nf=1, np=1, end.impute=FALSE)$imp[[1]])),33)
})

test_that("MICT - end.impute - one imp - sum miss", {
  skip_on_cran()
  library("TraMineR")
  seqgame <- seqdef(gameadd, var=1:4)
  expect_equal(sum(is.na(seqimpute(seqgame, regr="multinom", m=1, timing=F, nf=1, np=1)$imp[[1]])),0)
})

test_that("MICT timing - end.impute - one imp - sum miss", {
  skip_on_cran()
  library("TraMineR")
  seqgame <- seqdef(gameadd, var=1:4)
  expect_equal(sum(is.na(seqimpute(seqgame, regr="multinom", m=1, timing=T, nf=1, np=1, end.impute=FALSE)$imp[[1]])),33)
})

test_that("MICT timing - end.impute - one imp - sum miss", {
  skip_on_cran()
  library("TraMineR")
  seqgame <- seqdef(gameadd, var=1:4)
  expect_equal(sum(is.na(seqimpute(seqgame, regr="multinom", m=1, timing=T, nf=1, np=1)$imp[[1]])),0)
})