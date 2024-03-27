# 1. type of input ####
## a. factor #####
gameadd.traj <- gameadd[,1:4]
test_that("factor provided", {
  skip_on_cran()
  expect_no_error(seqTrans(gameadd.traj,trans="yes->no"))
})

## b. character #####
gameadd.traj.ch <- gameadd.traj
for(i in 1:ncol(gameadd.traj.ch)){
  gameadd.traj.ch[,i] <- as.character(gameadd.traj.ch[,i])
}

test_that("character provided", {
  skip_on_cran()
  expect_no_error(seqTrans(gameadd.traj.ch,trans="yes->no"))
})

## c. numeric #####
gameadd.num <- gameadd[,1:4]
for(i in 1:ncol(gameadd.num)){
  levels(gameadd.num[,i]) <- c(1,2)
  gameadd.num[,i] <- as.numeric(gameadd.num[,i])
}

test_that("numeric provided", {
  skip_on_cran()
  expect_no_error(seqTrans(gameadd.num,trans="1->2"))
})

# 2. test var ####
test_that("var provided", {
  skip_on_cran()
  expect_no_error(seqTrans(gameadd,var=1:4,trans="yes->no"))
})

# 3. no arrow ####
test_that("var provided", {
  skip_on_cran()
  expect_error(seqTrans(gameadd,var=1:4,trans="tmp-no"))
})

# 4. more than transitions ####
test_that("more than transitions", {
  skip_on_cran()
  expect_no_error(seqTrans(gameadd,var=1:4,trans=c("yes->no","no->yes")))
})

# 5. no impossible transition ####
gameadd.noimp <- gameadd.traj
for(i in 2:4){
  levels(gameadd.noimp[,i]) <- c("no","no")
}

test_that("no transitions", {
  skip_on_cran()
  expect_message(seqTrans(gameadd.noimp,var=1:4,trans=c("no->yes")),"Your dataset has no impossible transitions!")
})