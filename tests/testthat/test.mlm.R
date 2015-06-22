## test coercion of optmatch to matrix.csr:
## - gives proper mean diffs w/ matched pairs
## - gives proper mean diffs w/ matches w/ mult controls
## - gives proper mean diffs w/ matches w/ varying num controls
## - gives mean diffs of 0 for matched sets w/ no tx or no ctl
## - associates levels(from)[i] with row i of result

library(optmatch)
library(testthat)
context("Matched Diference Models (mlm)")

test_that("parseMatchingProblem", {

  data(nuclearplants)
  mmm <- fullmatch(pr ~ t1 + t2 + cap, data = nuclearplants)

  n2 <- cbind(nuclearplants, mmm)
  res <- parseMatchingProblem(cost ~ pr*pt + mmm, n2)
  
  expect_is(res$fmla, "formula")
  expect_is(res$mf, "data.frame")
  expect_is(res$match, "optmatch")

  expect_equivalent(res$match, mmm)
  expect_equal(dim(res$mf)[2], 3) # cost, pr, pt, but no mmm

  mmm2 <- fullmatch(pr ~ t1, data = nuclearplants)
  n3 <- cbind(n2, mmm2)
  res2 <- parseMatchingProblem(cost ~ pr + mmm, n3)
  expect_equivalent(res2$match, mmm)
  expect_error(parseMatchingProblem(cost ~ pr + mmm + mmm2, n3), "one")
               
  expect_error(parseMatchingProblem(cost ~ pr, n3), "include")
})

test_that("optmatch -> matrix.csr", {
  data(nuclearplants)

  ### start with the most straightforward case: everyone is matched
  mmm <- fullmatch(pr ~ t1 + t2 + cap, data = nuclearplants)

  csr <- as(mmm, "matrix.csr")

  expect_is(csr, "matrix.csr")
  expect_equal(dim(csr), c(nlevels(mmm), length(mmm)))

  csrm <- as.matrix(csr) # cast it back to dense matrix for testing
  grps <- apply(csrm, 2, function(col) {
    tmp <- which(col != 0)
    if (length(tmp) == 0) {
      return(NA)
    }
    return(tmp)
  })

  tmp <- table(grps, mmm)
  expect_true(all(diag(tmp) != 0))
  diag(tmp) <- 0
  expect_true(all(tmp == 0))

  expect_true(all(rowSums(csrm) == 0))

  expect_equal(as.vector(csr %*% rep(1, length(mmm))), rep(0, nlevels(mmm)))

  ### next case: not everyone is match
  mmm <- pairmatch(pr ~ t1 + t2 + cap, data = nuclearplants)

  csr <- as(mmm, "matrix.csr")

  expect_is(csr, "matrix.csr")
  expect_equal(dim(csr), c(nlevels(mmm), length(mmm)))

  csrm <- as.matrix(csr) # cast it back to dense matrix for testing
  grps <- apply(csrm, 2, function(col) {
    tmp <- which(col != 0)
    if (length(tmp) == 0) {
      return(NA)
    }
    return(tmp)
  })

  tmp <- table(grps, mmm)
  expect_true(all(diag(tmp) != 0))
  diag(tmp) <- 0
  expect_true(all(tmp == 0))

  expect_true(all(rowSums(csrm) == 0))

  expect_equal(as.vector(csr %*% rep(1, length(mmm))), rep(0, nlevels(mmm)))

  ### knock a specific group
  mmm <- fullmatch(pr ~ t1 + t2 + cap, data = nuclearplants)
  mmm <- mmm[mmm != levels(mmm)[1]]

  csr <- as(mmm, "matrix.csr")

  expect_is(csr, "matrix.csr")
  expect_equal(dim(csr), c(nlevels(mmm), length(mmm)))

  csrm <- as.matrix(csr) # cast it back to dense matrix for testing
  grps <- apply(csrm, 2, function(col) {
    tmp <- which(col != 0)
    if (length(tmp) == 0) {
      return(NA)
    }
    return(tmp)
  })

  tmp <- table(grps, mmm)
  expect_true(all(tmp[, 1] == 0))
  expect_true(all(diag(tmp[,-1]) != 0))

  expect_true(all(rowSums(csrm) == 0))

  expect_equal(as.vector(csr %*% rep(1, length(mmm))), rep(0, nlevels(mmm)))
  
  ### remove only a treated member
  mmm <- fullmatch(pr ~ t1 + t2 + cap, data = nuclearplants)
  mmm <- mmm[!(mmm == levels(mmm)[1] & nuclearplants$pr == 1)]

  csr <- as(mmm, "matrix.csr")

  expect_is(csr, "matrix.csr")
  expect_equal(dim(csr), c(nlevels(mmm), length(mmm)))

  csrm <- as.matrix(csr) # cast it back to dense matrix for testing
  grps <- apply(csrm, 2, function(col) {
    tmp <- which(col != 0)
    if (length(tmp) == 0) {
      return(NA)
    }
    return(tmp)
  })

  tmp <- table(grps, mmm)
  expect_true(all(tmp[, 1] == 0))
  expect_true(all(diag(tmp[,-1]) != 0))

  expect_true(all(rowSums(csrm) == 0))

  expect_equal(as.vector(csr %*% rep(1, length(mmm))), rep(0, nlevels(mmm)))
})

test_that("mlm", {
  
  data(nuclearplants)
  mmm <- fullmatch(pr ~ t1 + t2 + cap, data = nuclearplants)

  n2 <- cbind(nuclearplants, mmm)

  # these are non-failure tests. just checking for errors
  mlm(cost ~ mmm, data = n2)
  mlm(cost ~ mmm + 1, data = n2)
  mt1t2 <- mlm(cost ~ t1 + t2 + mmm, data = n2)
  mlm(cost ~ t1 + t2 + mmm, data = n2, fit.type = "robust")
  mlm(cost ~ t1 + t2 + mmm, data = n2, ms.weights = harmonic)

  expect_true(all(names(coef(mt1t2)) %in% c("(Treatment)", "t1", "t2")))
 
  ppp <- pairmatch(pr ~ t1 + t2 + cap, data = nuclearplants)
  n3 <- cbind(nuclearplants, ppp)

  mlm(cost ~ t1 + ppp, data = n3)
  mlm(cost ~ t1 + ppp, data = n3, fit.type = "robust")
  mlm(cost ~ t1 + ppp, data = n3, ms.weights = harmonic)

  # checking to see if factors, interactions work
  nuclearplants$f <- as.factor(sample(c("a","b","c"), dim(nuclearplants)[1], replace = T))
  n4 <- cbind(nuclearplants, mmm)
  m4 <- mlm(cost ~ t1 + t2 + f + mmm, data = n4)

  # using the contrasts.arg
  # use n4 again, which has a factor called if
  m5 <- mlm(cost ~ t1 + t2 + f + mmm, data = n4, contrasts.arg = list(f = contr.treatment(3, contrasts = FALSE)))
  expect_equal(length(coef(m4)) + 1, length(coef(m5))) # expect to have one level for each of level of f: a, b, c; instead of a reference category

  # Treatment should be excluded when formula has -1 in it (ie. no intercept)
  m6 <- mlm(cost ~ t1 + mmm, data = n2)
  expect_equal(length(coef(m6)), 2)

  m7 <- mlm(cost ~ t1 + mmm - 1, data = n2)
  expect_equal(length(coef(m7)), 1)
})

test_that("variances", {
  ppp <- pairmatch(pr ~ t1 + t2 + cap, data = nuclearplants)
  n3 <- cbind(nuclearplants, ppp)

  # simple model
  modelm <- mlm(cost ~ ppp, data = n3)
  modell <- lm(cost ~ pr + ppp, data = n3)

  expect_equal(coef(summary(modelm))["(Treatment)", ],
               coef(summary(modell))["pr", ])

  expect_equal(vcov(modell)["pr", "pr"], vcov(modelm)[1,1])
  
})

test_that("model.matrix.mlm", {
  ppp <- pairmatch(pr ~ t1 + t2 + cap, data = nuclearplants)
  n3 <- cbind(nuclearplants, ppp)

  mm <- mlm(cost ~ pt + ppp, data = n3)
  mx <- model.matrix(mm)

  expect_equal(dim(mx), c(nlevels(ppp), 2))
  
})

test_that("matches without a treated or control unit", {

  data(nuclearplants)

  ### start with the most straightforward case: everyone is matched
  mmm <- pairmatch(pr ~ t1 + t2 + cap, data = nuclearplants)
  expect_true(!is.na(mmm["I"])) # when i ran this at the terminal "I" was matched.

  mmm["I"] <- NA

  # non-crash "tests"
  res <- mlm(cost ~ cap + m, data = cbind(nuclearplants, m = mmm))
  res <- mlm(cost ~ m, data = cbind(nuclearplants, m = mmm))
})

