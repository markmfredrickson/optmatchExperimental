context("InfintySparseMatrix tests")

## test_that("ISM indexing", {

##   data(nuclearplants)
##   m <- match_on(pr ~ cost, data=nuclearplants, caliper=1)

##   expect_equal(dim(m[1:3,1]), c(3,1))
##   expect_equal(dim(m[1,1:3]), c(1,3))
##   expect_equal(dim(m[c("A", "C"), c(4,7,1,2:4)]), c(2, 5))

##   expect_equal(length(m[1:3]), 3)
##   expect_equal(length(m[c("A", "a")]), 2)

##   # These tests fail
##   # expect_equal(dim(m[1:3, ]), c(3, 22))
##   # expect_equal(m, m[])
##   # expect_equal(m, m[,])
## })


test_that("summary for ISM", {
  set.seed(1)
  d <- data.frame(z=rep(0:1, each=5),
                  x=rnorm(10))
  rownames(d) <- letters[1:10]
  m1 <- match_on(z ~ x, data=d)
  sm1 <- summary(m1)

  m2 <- m1 + caliper(m1, width=1)
  sm2 <- summary(m2)
  expect_true(is(sm2, "summary.InfinitySparseMatrix"))
  expect_true(is.list(sm2))
  expect_equal(sm2$total$treatment, 5)
  expect_equal(sm2$total$control, 5)
  expect_equal(sm2$total$matchable, 12)
  expect_equal(sm2$total$unmatchable, 25-12)
  expect_equal(length(sm2$matchable$treatment), 5)
  expect_equal(length(sm2$matchable$control), 4)
  expect_equal(sm2$unmatchable$treatment, character(0))
  expect_equal(sm2$unmatchable$control, "d")
  expect_true(is(sm2$distances, "summaryDefault"))

  m3 <- m2
  m3[1:2] <- Inf
  sm3 <- summary(m3)
  ##  expect_equal(sm3$total$matchable, 10)
  ##  expect_equal(sm3$total$unmatchable, 25-10)
  ##  A bug in num_eligible_matches in optmatch was fixed. Until that
  ##  gets pushed up so that the newest version of optmatch is
  ##  installed, the above two tests will fail on `make test`. They
  ##  should work interactively (with `make load`).
  expect_equal(length(sm3$matchable$treatment), 4)
  expect_equal(length(sm3$matchable$control), 4)
  expect_equal(sm3$unmatchable$treatment, "f")
  expect_equal(sm3$unmatchable$control, "d")
  expect_true(is(sm3$distances, "summaryDefault"))
  expect_true(all(is.finite(sm3$distances)))

  m4 <- m1 + caliper(m1, width=.0001)
  sm4 <- summary(m4)
  expect_equal(sm4$matchable$treatment, character(0))
  expect_equal(sm4$matchable$control, character(0))


})

test_that("summary for BlockedISM", {
  set.seed(1)
  d <- data.frame(z=rep(0:1, each=5),
                  x=rnorm(10),
                  q=rep(c("a", "d"), times=5))
  rownames(d) <- letters[1:10]
  m1 <- match_on(z ~ x + strata(q), data=d, caliper=1)
  sm1 <- summary(m1)

  expect_true(is(sm1, "summary.BlockedInfinitySparseMatrix"))
  expect_true(is.list(sm1))
  expect_equal(length(sm1), 2)
  expect_equal(names(sm1), c("a", "d"))
  expect_true(is(sm1[["a"]], "summary.InfinitySparseMatrix"))
  expect_true(is(sm1[["d"]], "summary.InfinitySparseMatrix"))

})


test_that("summary for DenseMatrix", {
  set.seed(1)
  d <- data.frame(z=rep(0:1, each=5),
                  x=rnorm(10))
  rownames(d) <- letters[1:10]
  m1 <- match_on(z ~ x, data=d)
  sm1 <- summary(m1)
  expect_true(is(sm1, "summary.DenseMatrix"))
  expect_true(is.list(sm1))
  expect_equal(sm1$total$treatment, 5)
  expect_equal(sm1$total$control, 5)
  expect_equal(sm1$total$matchable, 25)
  expect_equal(sm1$total$unmatchable, 0)
  expect_equal(length(sm1$matchable$treatment), 5)
  expect_equal(length(sm1$matchable$control), 5)
  expect_equal(sm1$unmatchable$treatment, character(0))
  expect_equal(sm1$unmatchable$control, character(0))
  expect_true(is(sm1$distances, "summaryDefault"))

  m2 <- m1
  m2[1,] <- Inf
  sm2 <- summary(m2)
  expect_true(is(sm2, "summary.DenseMatrix"))
  expect_true(is.list(sm2))
  expect_equal(sm2$total$treatment, 5)
  expect_equal(sm2$total$control, 5)
  expect_equal(sm2$total$matchable, 20)
  expect_equal(sm2$total$unmatchable, 25-20)
  expect_equal(length(sm2$matchable$treatment), 4)
  expect_equal(length(sm2$matchable$control), 5)
  expect_equal(sm2$unmatchable$treatment, "f")
  expect_equal(sm2$unmatchable$control, character(0))
  expect_true(is(sm2$distances, "summaryDefault"))
})
