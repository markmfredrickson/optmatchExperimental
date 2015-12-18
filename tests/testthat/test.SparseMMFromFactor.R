library(optmatch)
library(SparseM)
library(testthat)
context("SparseMMFromFactor")

test_that("basic input/output", {

  f <- factor(c("a", "a", "b", "b"))

  smm <- SparseMMFromFactor(f)

  expect_true(is(smm, "matrix.csr"))

  expect_error(SparseMMFromFactor(as.numeric(f)))
  expect_error(SparseMMFromFactor(as.character(f)))

  smm2 <- SparseMMFromFactor(as.factor(as.character(f)))
  expect_identical(smm, smm2)

})

test_that("equivalent to model.matrix", {

  f <- factor(c("a", "a", "b", "b"))

  mm <- model.matrix( ~ f + 0, data=f)
  smm <- SparseMMFromFactor(f)

  expect_true(all.equal(mm, as.matrix(smm), check.attributes=FALSE))

  f2 <- factor(c("a", "a"))

  # model.matrix will choke here!
  # mm2 <- model.matrix( ~ f2 + 0, data=f2)
  smm2 <- SparseMMFromFactor(f2)

  expect_identical(as.matrix(smm2), matrix(c(1,1)))

})

test_that("NA handling", {

  f <- factor(c("a", "a", "b", NA))

  # model.matrix's behavior with NA's is different, so not going to do
  # a comparison
  expect_warning(smm <- SparseMMFromFactor(f))

  expect_identical(as.matrix(smm), matrix(c(1,1,0,0,0,0,1,0), ncol=2))

  class(f) <- c("optmatch", "factor")
  expect_silent(smm2 <- SparseMMFromFactor(f))

  expect_error(SparseMMFromFactor(factor(c(NA, NA))))

  expect_identical(smm, smm2)

})
