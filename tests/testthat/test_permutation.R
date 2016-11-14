context("Permutation")

test_that("Basic tests", {
  expect_that(rowSums(permutatedSamples(data=matrix(1:12,2), B=1000)), equals(rep(21,1000)))
  expect_that(ncol(permutatedSamples(data=matrix(1:12,2), B=1000)), equals(6))
})
