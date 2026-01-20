# 1. General testing ####
## a. past and future ####
gameadd.traj <- gameadd[, 1:4]

test_that("MICT works without error", {
  expect_no_error(seqimpute(gameadd.traj, niter=3))
})

test_that("MICT-timing works without error", {
  expect_no_error(seqimpute(gameadd.traj, timing = TRUE, niter=3))
})

test_that("MICT single imputation works without error", {
  expect_no_error(seqimpute(gameadd.traj, m = 1, niter=3))
})

test_that("MICT-timing multiple imputation works without error", {
  expect_no_error(seqimpute(gameadd.traj, timing = TRUE, m = 1, niter=3))
})


test_that("MICT works without error", {
  skip_on_cran()
  expect_equal(sum(is.na(seqimpute(gameadd, var = 1:4, m = 2, niter=3)$imp[[1]])), 0)
})


test_that("MICT works without error", {
  skip_on_cran()
  expect_equal(sum(is.na(seqimpute(gameadd, var = 1:4, timing = T, niter=3, m = 2)$imp[[1]])), 0)
})



# 14. more than two categories ####
library("TraMineR")
data(mvad)
mvad.miss <- seqaddNA(mvad[, 17:86])

# a. regr multinom####
test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(mvad.miss, regr = "multinom", m = 2, niter=2))
})

test_that("MICTt - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(mvad.miss, regr = "multinom", timing = TRUE, m = 2, niter=2))
})


# b. regr rf ####
test_that("MICT - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(mvad.miss, regr = "rf", m = 2, niter=2))
})

test_that("MICTt - multinom", {
  skip_on_cran()
  expect_no_error(seqimpute(mvad.miss, regr = "rf", timing = TRUE, m = 2, niter=2))
})