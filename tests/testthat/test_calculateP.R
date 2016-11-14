context("p-values")

test_that("Single values", {
  expect_that(calculateP(0, matrix(1:100,10)), equals(1))
  expect_that(calculateP(1, matrix(1:100,10)), equals(1))
  expect_that(calculateP(51, matrix(1:100,10)), equals(0.5))
  expect_that(calculateP(100, matrix(1:100,10)), equals(0.01))
  expect_that(calculateP(101, matrix(1:100,10)), equals(0))
})

test_that("Multiple values", {
  expect_that(calculateP(c(0,101,0,51), matrix(1:100,10)), equals(c(1,0,1,0.5)))
})
