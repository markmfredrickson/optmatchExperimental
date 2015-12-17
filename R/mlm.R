#' Matched set based treatment effects.
#'
#' Fits a regression model for matched data. Currently only linear
#' regression via `lm` is implemented.
#'
#' The data will be collapsed from individuals units to matched set
#' units. Weighting defined by the size of each matched set (so
#' weights are constant in pair-matching) is calculated as well. The
#' returned object is of class `lm` and should support all relevant
#' follow-up analysis.
#'
#' @title Ordinary least squares for matched differences
#' @param formula The right hand side should contain an `optmatch` object.
#' @param data Data that contains the variables in `formula`, missing data will be imputed using `fill.NAs`
#' @param ms.weights Function of 2 vector args `n.t`, `n.c`, sums of weights from treatment and control group members by matched set, returning vector of matched-set specific weights. `harmonic` returns harmonic means of n.t and n.c; `ett` simply returns n.t.
#' @param fit.type character string indicating type of fit. For now, only "lm", but may expand to include rlm
#' @param fit.control optional list of additional arguments to the fitter
#' @param na.action How NA's are treated. See `model.frame` for details.
#' @param contrasts.arg An optional list of contrast matrices that will be passed to \code{\link{model.matrix}}.
#' @param ... additional arguments passed to `model.frame`
#' @return object of class `lm`
#' @export
#' @author Ben B. Hansen, Mark M. Fredrickson
#' @import SparseM
#' @importFrom MASS rlm
mlm <- function(formula, data, ms.weights = ett, fit.type = "lm", fit.control = list(), na.action = na.pass, contrasts.arg = NULL, ...) {

  cl <- match.call()

  helper <- mlm_helper(formula, data, na.action = na.pass, contrasts.arg, ...)

  fit.weights <- ms.weights(helper$nt, helper$nc)

  if ((fit.type == "robust" || fit.type == "rlm")) {
    return(do.call(MASS::rlm, append(list(helper$X, helper$Y, weights = fit.weights), fit.control)))
  }

  # can't fit with rlm package. fall back to good old lm
  z <- do.call(lm.wfit, append(list(helper$X, helper$Y, w = fit.weights, offset = helper$offset), fit.control))
  class(z) <- c("lm", "mlm")
  z$offset <- helper$offset
  z$terms <- terms(update(helper$parsed$fmla, .~.-1), data = helper$noNAs)
  z$call <- cl
  z$x <- helper$X
  z$y <- helper$Y
  z$weights <- fit.weights
  return(z)
}

# Not exported. Just to avoid duplication between `mlm`, and `model.matrix`.
# See mlm for the arguments
##' Internal function for mlm. See `mlm` for details.
##'
##' @param formula See `mlm`.
##' @param data See `mlm`.
##' @param na.action See `mlm`.
##' @param contrasts.arg See `mlm`.
##' @param ... See `mlm`.
##' @return See `mlm`.
mlm_helper <- function(formula, data, na.action = na.pass, contrasts.arg = NULL, ...) {

  parsed <- parseMatchingProblem(formula, data, na.action, ...)
  outcome <- model.response(parsed$mf, "numeric")
  weights <- model.weights(parsed$mf)
  offset <- model.offset(parsed$mf)

  # this will be a nearly filled in model matrix (ie. all factors expanded), but without an intercept
  noNAs <- fill.NAs(parsed$fmla, parsed$mf, contrasts.arg = contrasts.arg)

  theMatch <- parsed$match

  checkNA <- function(i) {
    if (is.null(i)) {
      return(FALSE)
    }
    return(is.na(i))
  }

  remove <- !with(parsed,
                  checkNA(weights) |
                  is.na(outcome) |
                  is.na(theMatch))

  noNAs <- noNAs[remove, , drop = FALSE]
  if ("(Intercept)" %in% colnames(noNAs)) {
    noNAs <- noNAs[, -which("(Intercept)" %in% colnames(noNAs)), drop = FALSE]
  }

  if(parsed$oname %in% colnames(noNAs)) {
    noNAs <- noNAs[, -which(parsed$oname %in% colnames(noNAs)), drop = FALSE]
  }

  theMatch <- theMatch[remove]
  outcome <- outcome[remove]
  weights <- weights[remove]

  z <- attr(theMatch, "contrast.group")

  # Faster than original sapply; + 0 is to drop the intercept and
  # ensure all indicators are generated.
  nt <- as.vector(z%*%model.matrix(~ theMatch + 0))
  nc <- as.vector((!z)%*%model.matrix(~ theMatch + 0))

  missingTorC <- nt == 0 | nc == 0
  nt <- nt[!missingTorC]
  nc <- nc[!missingTorC]

  matchCsr <- as(theMatch, "matrix.csr")
  matchCsr <- matchCsr[!missingTorC, ]

  if (ncol(noNAs) == 0) {
    X <- matrix(1, ncol = 1, nrow = nrow(matchCsr))
    colnames(X) <- "(Treatment)"
  } else {
    # make the design matrix for the matched sets
    # switching back to dense representation since the we don't expect many zero's in the design matrix
    X <- as.matrix(matchCsr %*% as.matrix(noNAs))
    colnames(X) <- colnames(noNAs)
    # the rows are the matched sets (but we don't need to include those)

    if (has.intercept(parsed$fmla)) {
      X <- cbind("(Treatment)" = 1, X)
    }
  }

  Y <- matchCsr %*% outcome

  return(list(X = X, Y = Y, nt = nt, nc = nc, noNAs = noNAs, weights = weights, offset = offset, parsed = parsed, matchCsr = matchCsr))
}
# This next line should get put in the makeOptmatch.R file. It provides S4 compatability.
setOldClass(c("optmatch", "factor"))

# ROxygen doesn't like this block, so turning it off for now.

##' Sparse matrices with which to assemble treatment minus control differences by matched set
##'
##' @name as
##' @family optmatch
##' @param from An optmatch object
##' @return A matrix.csr object by which to left-multiply vectors
##' and model matrices in order to assemble matched differences.
##' @author Ben B Hansen
##' @import SparseM optmatch
setAs("optmatch", "matrix.csr", function(from) {
  # treatment variable, a logical
  zz <- optmatch:::toZ(attr(from, "contrast.group")) # can remove the explicit namespace when this goes in the optmatch pkg

  # vector of positions of treatment member(s), then
  # control group members
  pos.tc <- order(from, !zz)

  fromNum <- as.numeric(from)
  fromNoNA <- fromNum[!is.na(from)]
  # starting positions for rows of the csr matrix
  rowstarts <- as.integer(c(rep(1, min(fromNoNA)),
                             which(diff(fromNum[pos.tc])==1)+1,
                             length(fromNoNA)+1))

  # each row has 1st tx and then ctl, but we need to know how many of each
  n.t <- as.integer(table(from[zz, drop = FALSE]))
  n.c <- as.integer(table(from))-n.t

  # if either t or c unrepresented in a matched set, null out other group's contrib
  # (this can happen due to `from` having been subsetted, perhaps b/c of NAs elsewhere)
  tscale <- ifelse(n.c&n.t, 1/n.t, 0)
  cscale <- ifelse(n.c&n.t, -1/n.c, 0)

  # multipliers to go in positions pos.tc
  multipliers <- rep(as.vector(rbind(tscale, cscale)),
                     as.vector(rbind(n.t, n.c)) )

  new("matrix.csr",
      ra = multipliers,
      ja = pos.tc[1:length(fromNoNA)],
      ia = rowstarts,
      dimension = c(nlevels(from), length(from)))
})

ett <- function(n.t,n.c) n.t
harmonic <-  function (n.t, n.c) 2*(1/n.t + 1/n.c)^-1

#' Helper to parse a matched analysis from a formula and a data.frame.
#'
#' @param formula The formula.
#' @param data The data.frame containing terms in the formula.
#' @param na.action How NA's are treated. See `model.frame` for details.
#' @param ... Other arguments passed to `model.frame`.
#' @return A list with: `mf` a model frame stripped of the optmatch argument, `match` the matched factor, `oname` the name of the outcome variable, `fmla` an updated formula without the matching vector
parseMatchingProblem <- function(formula, data, na.action = na.pass, ...) {
  mf <- model.frame(formula, data, na.action = na.action, ...)

  isMatch <- sapply(mf, function(i) { inherits(i, "optmatch") })

  if (sum(isMatch) != 1) {
    stop("You must include precisely one matching in the formula.")
  }

  match <- mf[, isMatch, drop = TRUE]
  names(match) <- rownames(data)

  mname <- colnames(mf)[isMatch]

  newf <- update(formula, as.formula(paste(".~. -", mname)), data=mf)

  # now make a new model frame, using the reduce form of the formula
  mf <- model.frame(newf, data, na.action = na.action, ...)

  return(
      list(
          fmla = newf,
          mf = mf,
          match = match,
          oname = as.character(newf)[[2]]))
}

has.intercept <- function(fmla) {
  attr(terms(fmla), "intercept") == 1
}
