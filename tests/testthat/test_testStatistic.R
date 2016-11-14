context("Test statistic")
data(upsSpikeIn)
upsSpikeIn[1,1:3] <- c(0,0,0)
upsSpikeIn[1,4:6] <- c(0,0,0)
upsSpikeIn[2,1:3] <- c(0,0,0)
upsSpikeIn[2,4:6] <- c(4,4,4)

test_that("Unpaired", {
  stats <- testStatistic(upsSpikeIn[,1:3],upsSpikeIn[,4:6], paired=FALSE)
  expect_that(as.numeric(stats$d[1]), equals(0))
  expect_that(as.numeric(stats$d[2]), equals(4))
  expect_that(as.numeric(stats$s[1]), equals(0))
  expect_that(as.numeric(stats$s[2]), equals(0))
})

test_that("Paired", {
  stats <- testStatistic(upsSpikeIn[,1:3],upsSpikeIn[,4:6], paired=TRUE)
  expect_that(as.numeric(stats$d[1]), equals(0))
  expect_that(as.numeric(stats$d[2]), equals(4))
  expect_that(as.numeric(stats$s[1]), equals(0))
  expect_that(as.numeric(stats$s[2]), equals(0))
})
