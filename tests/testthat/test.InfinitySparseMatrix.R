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
  expect_equal(length(sm1), 4)
  expect_equal(names(sm1), c("a", "d", "overall", "matname"))
  expect_true(is(sm1[["a"]], "summary.InfinitySparseMatrix"))
  expect_true(is(sm1[["d"]], "summary.InfinitySparseMatrix"))

  # The current version of optmatch on CRAN has a bug in
  # num_eligible_matches. This isn't an issue in general, as using
  # load.R allows using the most recent version on the
  # repo. (E.g. works fine in `make interactive` and `make test`.)
  # However, `make check` starts a fresh version of R which I can't
  # figure out how to manipulate, so this will fail. So
  if (compareVersion(as.character(packageVersion("optmatch")), "0.9-5") == 1) {
    expect_true(all.equal(c(5,5,4,21),
                          unlist(sm1$overall$total),
                          check.attributes=FALSE))
  }

  # Alternate ways of calling bloks
  suma1 <- sm1[['a']]
  suma2 <- sm1$`a`
  suma3 <- sm1[[1]]
  expect_identical(suma1, suma2)
  expect_identical(suma1, suma3)


  expect_true(sm1$`a`$blockname == "a")
  expect_true(sm1$`d`$blockname == "d")

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

test_that("distanceSummary suppresses distance", {
  set.seed(1)
  d <- data.frame(z=rep(0:1, each=5),
                  x=rnorm(10),
                  q=rep(c("a", "d"), times=5))
  rownames(d) <- letters[1:10]
  m1 <- match_on(z ~ x, data=d)

  expect_true(!is.null(summary(m1)$distances))
  expect_true(!is.null(summary(m1, distanceSummary=TRUE)$distances))
  expect_true(is.null(summary(m1, distanceSummary=FALSE)$distances))

  m2 <- match_on(z ~ x, data=d, caliper=1)
  expect_true(!is.null(summary(m2)$distances))
  expect_true(!is.null(summary(m2, distanceSummary=TRUE)$distances))
  expect_true(is.null(summary(m2, distanceSummary=FALSE)$distances))

  m3 <- match_on(z ~ x + strata(q), data=d, caliper=1)
  sm3.1 <- summary(m3)
  expect_true(!is.null(sm3.1$overall$distances))
  expect_true(!is.null(sm3.1$a$distances))
  expect_true(!is.null(sm3.1$d$distances))

  sm3.2 <- summary(m3, distanceSummary=TRUE)
  expect_true(!is.null(sm3.2$overall$distances))
  expect_true(!is.null(sm3.2$a$distances))
  expect_true(!is.null(sm3.2$d$distances))

  sm3.3 <- summary(m3, distanceSummary=FALSE)
  expect_true(is.null(sm3.3$overall$distances))
  expect_true(is.null(sm3.3$a$distances))
  expect_true(is.null(sm3.3$d$distances))
})
