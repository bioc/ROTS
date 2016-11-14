context("Bootstrapping")

test_that("Basic tests", {
  expect_that(nrow(bootstrapSamples(data=NULL, B=1000, labels=c(1,1,1,2,2,2), paired=FALSE)), equals(1000))
  expect_that(ncol(bootstrapSamples(data=NULL, B=10, labels=c(1,1,1,1,2,2,2,2), paired=FALSE)), equals(8))
  expect_that(all(apply(bootstrapSamples(data=NULL, B=1000, labels=c(1,1,1,2,2,2), paired=TRUE)[,1:3],1, function(x) all(x<4))), equals(TRUE))
  expect_that(all(apply(bootstrapSamples(data=NULL, B=1000, labels=c(1,1,1,2,2,2), paired=TRUE)[,4:6],1, function(x) all(x>3))), equals(TRUE))
})
