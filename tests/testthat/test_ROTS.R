context("ROTS function")
data(upsSpikeIn)

test_that("Basic tests", {
  rots.out <- ROTS(data=upsSpikeIn, groups=c(0,0,0,1,1,1), B=100, K=500, seed=1)
  expect_that(rots.out, is_a("ROTS"))
  expect_that(length(grep("ups",names(rots.out$d))), equals(36))
  expect_that(all((rots.out$FDR >= 0) & (rots.out$FDR <= 1)), equals(TRUE))
  expect_that(all((rots.out$p >= 0) & (rots.out$p <= 1)), equals(TRUE))
})

test_that("Special cases", {
  case1 <- c(10,10,10,10,10,10)
  case2 <- c(0,0,0,0,0,0)
  case3 <- c(NA,10,10,NA,10,10)
  case4 <- c(NA,0,0,NA,0,0)
  case5 <- c(1,1,1,10,10,10)
  case6 <- c(1,NA,1,10,NA,10)
  rots.out <- ROTS(data=rbind(case1,case2,case3,case4,case5,case6,upsSpikeIn), groups=c(0,0,0,1,1,1), B=100 ,K=500 ,seed=1)
  expect_that(all(rots.out$FDR[1:4] > 0.99), equals(TRUE))
  expect_that(all(rots.out$p[1:4] > 0.99), equals(TRUE))
  expect_that(all(rots.out$logfc[1:4] == 0), equals(TRUE))
  expect_that(all(rots.out$FDR[5:6] < 0.01), equals(TRUE))
  expect_that(all(rots.out$p[5:6] < 0.01), equals(TRUE))
  expect_that(all(rots.out$logfc[5:6] != 0), equals(TRUE))
})

test_that("Comparison to t-test", {
  rots.out <- ROTS(data=upsSpikeIn, groups=c(0,0,0,1,1,1), a1=0, a2=1, paired=FALSE, B=100, K=500, seed=1)$d
  ttest.out <- genefilter::rowttests(as.matrix(upsSpikeIn),fac=factor(c(1,1,1,0,0,0)))$statistic
  expect_that(cor(rots.out,ttest.out), equals(1))
  rots.out <- ROTS(data=upsSpikeIn, groups=c(0,0,0,1,1,1), a1=0, a2=1, paired=TRUE, B=100, K=500, seed=1)$d
  ttest.out <- genefilter::rowttests(as.matrix(upsSpikeIn[,4:6]-upsSpikeIn[,1:3]),fac=factor(c(0,0,0)))$statistic
  expect_that(cor(rots.out,ttest.out), equals(1))
})
