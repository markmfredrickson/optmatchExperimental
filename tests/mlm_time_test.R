rm(list=ls())
library(optmatch)
library(optmatchExperimental)

## Leaving this small for `make check`, make it bigger if actually
## testing speed.
#n <- 4000
n <- 100

d <- data.frame(x=rbinom(n, 1, .4),
                y=rbinom(n, 1, .02),
                z=rbinom(n, 1, .5))

d$f <- fullmatch(z ~ x, data=d, max=2, min=1/2)
length(levels(d$f))

system.time(m <- mlm(y ~ x + f, data=d))
