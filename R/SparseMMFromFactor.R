##' Turn a factor variable into a sparse matrix of 0's and 1's, such that if observation i
##' has the jth level then there is a 1 at position (i,j) (but nowhere else in row i).
##'
##' NA's give rise to rows with no 1s.
##' As the result is only meaningful in the context of the SparseM package,
##' function requires that SparseM be loaded.
##' @title Sparse matrix dummy coding of a factor variable (omitting the intercept)
##' @param thefactor Factor variable, or object inheriting from class factor
##' @return Sparse csr matrix the columns of which are dummy variables for levels of thefactor
##' @import SparseM
##' @export
##' @author Ben Hansen
##' @examples
##' sparse_mod_matrix <-  SparseMMFromFactor(iris$Species)
##' mod_matrix <- model.matrix(~Species-1, iris)
##' all.equal(SparseM::as.matrix(sparse_mod_matrix),
##'           mod_matrix, check.attributes=FALSE)
SparseMMFromFactor <- function(thefactor) {
  stopifnot(inherits(thefactor, "factor"))
  theNA <- ##if (inherits(thefactor, "optmatch")) !matched(thefactor) else
    is.na(thefactor)

  if (all(theNA)) stop("No non-NA's in thefactor") else {
    if (any(theNA) && !inherits(thefactor, "optmatch")) warning("NA's found in thefactor.")
  }

  nlev <- nlevels(thefactor)
  nobs <- length(thefactor)
  theint <- as.integer(thefactor)
  if (any(theNA)) theint[theNA] <- 1L#nlev + 1L:sum(theNA)
  new("matrix.csr",
      ja=theint,
      ia=1L:(nobs+1L),
      ra=(1L-theNA),
      dimension = c(nobs, nlev) #+sum(theNA)
      )
}
