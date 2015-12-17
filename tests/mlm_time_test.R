rm(list=ls())
library(optmatch)
library(optmatchExperimental)
n <- 4000

d <- data.frame(x=rbinom(n, 1, .4),
                y=rbinom(n, 1, .02),
                z=rbinom(n, 1, .5))

d$f <- fullmatch(z ~ x, data=d, max=2, min=1/2)
length(levels(d$f))

system.time(m <- mlm(y ~ x + f, data=d))
