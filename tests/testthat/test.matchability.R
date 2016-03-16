context("matchability")

test_that("matchability", {
  # Woefully incomplete testing, until we finalize input/output
  set.seed(1)
  d <- data.frame(z=rep(0:1, each=5),
                  x=rnorm(10))
  rownames(d) <- letters[1:10]
  m <- match_on(z ~ x, data=d)
  m1 <- caliper(m, width=1)
  m2 <- caliper(m, width=Inf)
  m3 <- caliper(m, width=.1)
  m4 <- caliper(m, width=.0001)
  m5 <- as.matrix(m1)


  mm1 <- matchability(m1)
  mm2 <- matchability(m2)
  mm3 <- matchability(m3)
  mm4 <- matchability(m4)
  mm5 <- matchability(m5)
  expect_equal(mm1, mm5)
})
