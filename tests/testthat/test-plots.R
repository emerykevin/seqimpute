# 1. seqmissfplot ####
gameadd.traj <- gameadd[,1:4]

## a. basic functioning ####
test_that("function works", {
  skip_on_cran()
  expect_no_error(seqmissfplot(gameadd.traj))
})

## b. seq object ####
library("TraMineR")
gameadd.seq <- seqdef(gameadd.traj, right=NA)
test_that("function works", {
  skip_on_cran()
  expect_no_error(seqmissfplot(gameadd.seq))
})

## c. without complete trajectories ####
test_that("function works", {
  skip_on_cran()
  expect_no_error(seqmissfplot(gameadd.traj, with.complete = FALSE))
})

## d. with var ####
test_that("function works", {
  skip_on_cran()
  expect_no_error(seqmissfplot(gameadd, var=1:4))
})

## e. seqobject NA as code ####
gameadd.seq <- seqdef(gameadd.traj, right=NA, nr=NA)
test_that("function works", {
  skip_on_cran()
  expect_no_error(seqmissfplot(gameadd.seq))
})

# 2. seqmissIplot ####
## a. basic functioning ####
test_that("function works", {
  skip_on_cran()
  expect_no_error(seqmissIplot(gameadd.traj))
})


## b. seq object ####
library("TraMineR")
gameadd.seq <- seqdef(gameadd.traj, right=NA)
test_that("function works", {
  skip_on_cran()
  expect_no_error(seqmissIplot(gameadd.seq))
})

## c. without complete trajectories ####
test_that("function works", {
  skip_on_cran()
  expect_no_error(seqmissIplot(gameadd.traj, with.complete = FALSE))
})

## d. with var ####
test_that("function works", {
  skip_on_cran()
  expect_no_error(seqmissIplot(gameadd, var=1:4))
})

## e seq with NA code ####
gameadd.seq <- seqdef(gameadd.traj, right=NA, nr=NA)
test_that("function works", {
  skip_on_cran()
  expect_no_error(seqmissIplot(gameadd.seq))
})

# 3. seqmissimplic ####
gameadd.seq <- seqdef(gameadd.traj, right=NA)
test_that("function works", {
  skip_on_cran()
  expect_no_error(seqmissimplic(gameadd, var=1:4))
})

test_that("sequence works", {
  skip_on_cran()
  expect_no_error(seqmissimplic(gameadd.seq, var=1:4))
})