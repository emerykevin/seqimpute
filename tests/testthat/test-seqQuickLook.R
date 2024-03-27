gameadd.traj <- gameadd[,1:4]

# 1. basic tests #####
## a. np=1, nf=1 #####
test_that("works without error", {
  skip_on_cran()
  expect_no_error(seqQuickLook(gameadd.traj))
})


test_that("number of NA corresponds to the number in original dataset", {
  skip_on_cran()
  expect_equal(sum(seqQuickLook(gameadd.traj)$sumNAGaps),
      sum(is.na(gameadd.traj)))
})

## b. with Left-SLG ######
test_that("works without error", {
  skip_on_cran()
  expect_no_error(seqQuickLook(gameadd.traj, np=2))
})


## c. with Right-SLG ######
test_that("works without error", {
  skip_on_cran()
  expect_no_error(seqQuickLook(gameadd.traj, nf=2))
})

## d. with Both-SLG ######
test_that("works without error", {
  skip_on_cran()
  expect_no_error(seqQuickLook(gameadd.traj, np=2, nf=2))
})

# 2. use of var #####
test_that("works without error", {
  skip_on_cran()
  expect_no_error(seqQuickLook(gameadd, var=1:4, nf=2))
})

# 3. type of data provided #####
## a. character ####
gameadd.traj.ch <- gameadd.traj
for(i in 1:ncol(gameadd.traj.ch)){
  gameadd.traj.ch[,i] <- as.character(gameadd.traj.ch[,i])
}
test_that("works without error", {
  skip_on_cran()
  expect_no_error(seqQuickLook(gameadd.traj.ch))
})


test_that("number of NA corresponds to the number in original dataset", {
  skip_on_cran()
  expect_equal(sum(seqQuickLook(gameadd.traj.ch)$sumNAGaps),
               sum(is.na(gameadd.traj.ch)))
})

## b. numeric ####
gameadd.num <- gameadd[,1:4]
for(i in 1:ncol(gameadd.num)){
  levels(gameadd.num[,i]) <- c(1,2)
  gameadd.num[,i] <- as.numeric(gameadd.num[,i])
}

test_that("works without error", {
  skip_on_cran()
  expect_no_error(seqQuickLook(gameadd.num))
})


test_that("number of NA corresponds to the number in original dataset", {
  skip_on_cran()
  expect_equal(sum(seqQuickLook(gameadd.num)$sumNAGaps),
               sum(is.na(gameadd.num)))
})


## c. sequence ####
library("TraMineR")
gameadd.seq <- seqdef(gameadd,1:4,right=NA)

test_that("works without error", {
  skip_on_cran()
  expect_no_error(seqQuickLook(gameadd.seq))
})


test_that("number of NA corresponds to the number in original dataset", {
  skip_on_cran()
  expect_equal(sum(seqQuickLook(gameadd.seq)$sumNAGaps),
               sum(is.na(gameadd[,1:4])))
})

# 4. no MD #####
imp <- seqimpute(gameadd.traj, m=1)
test_that("number of NA corresponds to the number in original dataset", {
  skip_on_cran()
  expect_equal(sum(seqQuickLook(imp$imp[[1]])$sumNAGaps), 0)
})

