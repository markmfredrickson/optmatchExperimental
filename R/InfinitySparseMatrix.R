## #' Indexing ISM's
## #'
## #' @param x An InfinitySparseMatrix object.
## #' @param i Row indices.
## #' @param j Col indices. If NULL, then indexing x@.Data occurs.
## #' @param ... Other arguments
## #' @param drop Ignored.
## #' @return The subset.
## #' @export
## setMethod("[", "InfinitySparseMatrix",
##           function(x, i, j=NULL, ..., drop=TRUE) {
##             if (is.null(j)) {
##               return(x@.Data[i, drop=drop, ...])
##             } else {
##               makelogical <- function(index, rowcol) {
##                 switch(class(index),
##                       "numeric" = (1:dim(x)[rowcol]) %in% index,
##                       "integer" = (1:dim(x)[rowcol]) %in% index,
##                       "character" = dimnames(x)[[rowcol]] %in% index,
##                       "logical" = index,
##                       stop("Unrecognized class"))
##               }

##               subi <- makelogical(i, 1)
##               subj <- makelogical(j, 2)

##               subset(x, subset=subi, select=subj)
##             }
##           }
##           )


##' Given a distance matrix which has potentially been calipered,
##' returns which observations are matchable and unmatchable
##'
##' @param object ISM, BISM or DenseMatrix
##' @param ... Ignored.
##' @return List of lists, $matchable$control, $matchable$treatment,
##'   $unmatchable$control andd $unmatchable$treatment
##' @export
##' @name summary.ism
summary.InfinitySparseMatrix <- function(object, ...) {
  if (is(object, "BlockedInfinitySparseMatrix")) {
    out <- lapply(levels(object@groups),
                  function(x) {
                    ism <- subset(object,
                                  subset=object@rownames %in% names(object@groups[object@groups == x]),
                                  select=object@colnames %in% names(object@groups[object@groups == x]))
                    summary(ism)
                  })
    names(out) <- levels(object@groups)
    class(out) <- "summary.BlockedInfinitySparseMatrix"
    return(out)
  }

  finitedata <- is.finite(object@.Data)
  mtreat <- 1:dim(object)[1] %in% sort(unique(object@rows[finitedata]))
  mcontrol   <- 1:dim(object)[2] %in% sort(unique(object@cols[finitedata]))
  distances <- summary(object@.Data[finitedata])

  out <- internal.summary.helper(object, mtreat, mcontrol, distances)

  class(out) <- "summary.InfinitySparseMatrix"
  out
}

##' @export
##' @rdname summary.ism
summary.DenseMatrix <- function(object, ...) {
  mtreat <- apply(object, 1, function(x) any(is.finite(x)))
  mcontrol <- apply(object, 2, function(x) any(is.finite(x)))
  distances <- summary(as.vector(object))

  out <- internal.summary.helper(object, mtreat, mcontrol, distances)

  class(out) <- "summary.DenseMatrix"
  out
}

internal.summary.helper <- function(x,
                                    matchabletxt,
                                    matchablectl,
                                    distances) {
  out <- list()
  d <- dim(x)

  # Size of treatment and control groups
  out$total$treatment <- d[1]
  out$total$control <- d[2]

  # Count of eligble and ineligible pairs.
  num_elig <- num_eligible_matches(x)[[1]]
  num_inelig <- prod(d) - num_elig
  out$total$matchable <- num_elig
  out$total$unmatchable <- num_inelig

  out$matchable$treatment <- rownames(x)[matchabletxt]
  out$matchable$control <- colnames(x)[matchablectl]
  out$unmatchable$treatment <- rownames(x)[!matchabletxt]
  out$unmatchable$control <- colnames(x)[!matchablectl]

  out$distances <- distances
  return(out)
}

##' @export
print.summary.InfinitySparseMatrix <- function(x, ...) {
  ### NOT UPDATED
  if (x$total$unmatchable == 0) {
    cat(paste("All", sum(unlist(x$total)), "matches are eligible.\n"))
  } else {
    cat(paste("Out of", sum(unlist(x$total)),
              "total potential eligible matches,", x$total$matchable,
              "are eligible for matching and", x$total$unmatchable,
              "matchings are prohibited.\n\n"))
  }

  if (x$total$unmatchable > 0) {
    if (length(x$unmatchable$treatment) > 0) {
      cat(paste("The following treatment group members are ineligible for any matches:\n"))
      cat(paste(x$unmatchable$treatment, collapse=", "))
      cat("\n\n")
    }
    if (length(x$unmatchable$control) > 0) {
      cat(paste("The following control group members are ineligible for any matches:\n"))
      cat(paste(x$unmatchable$control, collapse=", "))
      cat("\n\n")
    }
  }
}

##' @export
print.summary.DenseMatrix <- function(x, ...) {
  ### NOT UPDATED
  cat(paste0("All ", x$total, " potential matches (", x$dim[1],
             " treatment members and ", x$dim[2],
             " control members) are eligible.\n"))
}
