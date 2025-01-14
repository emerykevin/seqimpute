# 1. fromseqimp #####
gameadd.traj <- gameadd[, 1:4]

# a. basic testing #####
test_that("non seqimp object", {
  skip_on_cran()
  expect_error(fromseqimp(gameadd.traj))
})

test_that("wrong format", {
  skip_on_cran()
  expect_error(fromseqimp(gameadd.traj, format = "hello"))
})

test_that("wrong format", {
  skip_on_cran()
  expect_error(fromseqimp(gameadd.traj, format = 1))
})

## a.1 multiple imputation ####
imp <- seqimpute(data = gameadd, var = 1:4, m = 2)
test_that("stacked - include FALSE", {
  skip_on_cran()
  expect_no_error(fromseqimp(imp, format = "stacked", include = FALSE))
})

test_that("long - include FALSE", {
  skip_on_cran()
  expect_no_error(fromseqimp(imp, format = "long", include = FALSE))
})

test_that("long - include FALSE", {
  skip_on_cran()
  expect_no_error(fromseqimp(imp, format = "mids", include = FALSE))
})


test_that("stacked - include TRUE", {
  skip_on_cran()
  expect_no_error(fromseqimp(imp, format = "stacked", include = TRUE))
})

test_that("long - include TRUE", {
  skip_on_cran()
  expect_no_error(fromseqimp(imp, format = "long", include = TRUE))
})

test_that("long - include TRUE", {
  skip_on_cran()
  expect_no_error(fromseqimp(imp, format = "mids", include = TRUE))
})

# a.2 single imputation #####
imp <- seqimpute(data = gameadd, var = 1:4, m = 1)
test_that("stacked - include FALSE", {
  skip_on_cran()
  expect_no_error(fromseqimp(imp, format = "stacked", include = FALSE))
})

test_that("long - include FALSE", {
  skip_on_cran()
  expect_no_error(fromseqimp(imp, format = "long", include = FALSE))
})

test_that("long - include FALSE", {
  skip_on_cran()
  expect_no_error(fromseqimp(imp, format = "mids", include = FALSE))
})


test_that("stacked - include TRUE", {
  skip_on_cran()
  expect_no_error(fromseqimp(imp, format = "stacked", include = TRUE))
})

test_that("long - include TRUE", {
  skip_on_cran()
  expect_no_error(fromseqimp(imp, format = "long", include = TRUE))
})

test_that("long - include TRUE", {
  skip_on_cran()
  expect_no_error(fromseqimp(imp, format = "mids", include = TRUE))
})

# b. wrong format #####
test_that("broad as format - wrong", {
  skip_on_cran()
  expect_error(fromseqimp(imp, format = "broad", include = TRUE))
})

# 2. seqcomplete #####
## a. data.frame
test_that("test if working", {
  skip_on_cran()
  expect_no_error(seqcomplete(gameadd.traj))
})

test_that("no missing values", {
  skip_on_cran()
  expect_equal(sum(is.na(seqcomplete(gameadd.traj))), 0)
})

test_that("number of trajectories kept correct", {
  skip_on_cran()
  expect_equal(nrow(seqcomplete(gameadd.traj)), 241)
})

test_that("correct object output", {
  skip_on_cran()
  expect_s3_class(seqcomplete(gameadd.traj), "data.frame")
})

## b. seqdef
library("TraMineR")

gameadd.traj.seq <- seqdef(gameadd.traj, right = NA)
test_that("test if working", {
  skip_on_cran()
  expect_no_error(seqcomplete(gameadd.traj.seq))
})

test_that("no missing values", {
  skip_on_cran()
  expect_equal(sum(is.na(seqcomplete(gameadd.traj.seq))), 0)
})

test_that("number of trajectories kept correct", {
  skip_on_cran()
  expect_equal(nrow(seqcomplete(gameadd.traj.seq)), 241)
})

test_that("correct object output", {
  skip_on_cran()
  expect_s3_class(seqcomplete(gameadd.traj.seq), "stslist")
})


# 3. seqwithmiss #####
## a. data.frame
test_that("test if working", {
  skip_on_cran()
  expect_no_error(seqwithmiss(gameadd.traj))
})

test_that("number of trajectories kept correct", {
  skip_on_cran()
  expect_equal(nrow(seqwithmiss(gameadd.traj)), 259)
})

test_that("correct object output", {
  skip_on_cran()
  expect_s3_class(seqwithmiss(gameadd.traj), "data.frame")
})

## b. seqdef
gameadd.traj.seq <- seqdef(gameadd.traj, right = NA)
test_that("test if working", {
  skip_on_cran()
  expect_no_error(seqwithmiss(gameadd.traj.seq))
})

test_that("number of trajectories kept correct", {
  skip_on_cran()
  expect_equal(nrow(seqwithmiss(gameadd.traj.seq)), 259)
})

test_that("correct object output", {
  skip_on_cran()
  expect_s3_class(seqwithmiss(gameadd.traj.seq), "stslist")
})

gameadd.traj.seq.NA <- seqdef(gameadd.traj, right = NA, nr = NA)
test_that("number of trajectories kept correct", {
  skip_on_cran()
  expect_equal(nrow(seqwithmiss(gameadd.traj.seq.NA)), 259)
})

# 4. addcluster ####
cluster <- c(rep(1, 300), rep(2, 300), rep(1, 300), rep(2, 100))

test_that("correct object output", {
  skip_on_cran()
  imp <- seqimpute(data = gameadd, var = 1:4, m = 2)
  expect_s3_class(addcluster(imp, cluster), "seqimp")
})

test_that("error not right length cluster", {
  skip_on_cran()
  imp <- seqimpute(data = gameadd, var = 1:4, m = 1)
  expect_error(addcluster(imp, cluster))
})


test_that("correct object output", {
  skip_on_cran()
  imp <- seqimpute(data = gameadd, var = 1:4, m = 2)
  cluster <- matrix(cluster, nrow = 500, ncol = 2)
  expect_s3_class(addcluster(imp, cluster), "seqimp")
})

test_that("no error plot", {
  skip_on_cran()
  imp <- seqimpute(data = gameadd, var = 1:4, m = 2)
  cluster <- matrix(cluster, nrow = 500, ncol = 2)
  expect_no_error(plot(addcluster(imp, cluster)))
})

test_that("no error print", {
  skip_on_cran()
  imp <- seqimpute(data = gameadd, var = 1:4, m = 2)
  cluster <- matrix(cluster, nrow = 500, ncol = 2)
  expect_no_error(print(addcluster(imp, cluster)))
})

test_that("no error summary", {
  skip_on_cran()
  imp <- seqimpute(data = gameadd, var = 1:4, m = 2)
  cluster <- matrix(cluster, nrow = 500, ncol = 2)
  expect_no_error(summary(addcluster(imp, cluster)))
})

test_that("error not right dimension cluster", {
  skip_on_cran()
  imp <- seqimpute(data = gameadd, var = 1:4, m = 1)
  cluster <- matrix(cluster, nrow = 500, ncol = 2)
  expect_error(addcluster(imp, cluster))
})

test_that("not a seqimp obect", {
  skip_on_cran()
  imp <- seqimpute(data = gameadd, var = 1:4, m = 1)
  cluster <- matrix(cluster, nrow = 500, ncol = 2)
  expect_error(
    addcluster(gameadd.traj, cluster),
    "impdata is not a seqimp object"
  )
})
