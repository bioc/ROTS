context("Bigger N")

test_that("Single values", {
  expect_that(biggerN(0,1:100), equals(100))
  expect_that(biggerN(1,1:100), equals(100))
  expect_that(biggerN(99,1:100), equals(2))
  expect_that(biggerN(100,1:100), equals(1))
  expect_that(biggerN(101,1:100), equals(0))
  expect_that(biggerN(0,100:1), equals(100))
  expect_that(biggerN(1,100:1), equals(100))
  expect_that(biggerN(99,100:1), equals(2))
  expect_that(biggerN(100,100:1), equals(1))
  expect_that(biggerN(101,100:1), equals(0))
})

test_that("Multiple values", {
  expect_that(biggerN(c(1,2,3),1:100), equals(c(98,99,100)))
  expect_that(biggerN(c(3,2,1),1:100), equals(c(98,99,100)))
  expect_that(biggerN(c(1,1,1),1:100), equals(c(100,100,100)))
})

test_that("Special values", {
  expect_that(biggerN(NA,1:100), equals(100))
  expect_that(biggerN(Inf,1:100), equals(0))
  expect_that(biggerN(c(NA,NA,200),1:100), equals(c(0,100,100)))
})
