context("FDR")

test_that("Basic tests", {
  expect_that(calculateFDR(0.9, matrix(rep(1,100),10), FALSE), equals(1))
  expect_that(calculateFDR(1.0, matrix(rep(1,100),10), FALSE), equals(1))
  expect_that(calculateFDR(1.1, matrix(rep(1,100),10), FALSE), equals(0))
})
