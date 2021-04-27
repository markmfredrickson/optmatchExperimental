context("InfintySparseMatrix tests")

test_that("ISM indexing", {

  data(nuclearplants)
  m <- match_on(pr ~ cost, data = nuclearplants, caliper = 1)

  # [X, X]
  expect_equal(dim(m[1:3,2:3]), c(3,2))
  expect_equal(dim(m[3:2,4:2]), c(2,3))
  expect_equal(dim(m[c("A", "C"), c(4,7,1,2:4)]), c(2, 5))

  # [X]
  expect_equal(length(m[1:3]), 3)
  expect_equal(length(m[c("A", "a")]), 2)

  # [X,] or [,X]
  expect_equal(dim(m[1:3, ]), c(3, 22))
  expect_equal(dim(m[, 5:3]), c(10, 3))

  # []
  expect_equal(m, m[])

  # [,]
  m2 <- m[,]
  m@call <- NULL
  m2@call <- NULL
  expect_equal(m, m2)

  # Ignoring drop
  expect_equal(m[1:3, 2:3, drop = TRUE ], m[1:3, 2:3])
  expect_equal(m[1:3, 2:3, drop = FALSE], m[1:3, 2:3])
  expect_equal(m[1:3, , drop = TRUE ], m[1:3, ])
  expect_equal(m[1:3, , drop = FALSE], m[1:3, ])
  expect_equal(m[, 1:3, drop = TRUE ], m[, 1:3])
  expect_equal(m[, 1:3, drop = FALSE], m[, 1:3])
  expect_equal(m[, , drop = TRUE ], m[, ])
  expect_equal(m[, , drop = FALSE], m[, ])
  expect_equal(m[, drop = TRUE ], m[, ])
  expect_equal(m[, drop = FALSE], m[, ])
  expect_equal(m[drop = TRUE ], m[])
  expect_equal(m[drop = FALSE], m[])
})
