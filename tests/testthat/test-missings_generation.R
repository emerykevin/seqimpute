imp <- seqimpute(gameadd, var=1:4, m=1)
game.full <- imp$imp[[1]]


test_that("factor provided", {
  skip_on_cran()
  expect_no_error(seqaddNA(game.full))
})

# b. var ####
test_that("var", {
  skip_on_cran()
  expect_no_error(seqaddNA(game.full, var=1:4))
})

# c. only traj ####
test_that("var", {
  skip_on_cran()
  expect_no_error(seqaddNA(game.full, var=1:3, only.traj = FALSE))
})

test_that("var", {
  skip_on_cran()
  expect_no_error(seqaddNA(game.full, var=1:3, only.traj = TRUE))
})

# d. state high ####
test_that("var", {
  skip_on_cran()
  expect_no_error(seqaddNA(game.full, states.high = "no"))
})