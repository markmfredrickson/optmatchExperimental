context("ISM indexing")

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
