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

  finitedata <- is.finite(object@.Data)
  mtreat <- 1:dim(object)[1] %in% sort(unique(object@rows[finitedata]))
  mcontrol  <- 1:dim(object)[2] %in% sort(unique(object@cols[finitedata]))
  distances <- summary(object@.Data[finitedata])

  out <- internal.summary.helper(object, mtreat, mcontrol, distances)
  out$matname <- deparse(substitute(object))

  class(out) <- "summary.InfinitySparseMatrix"
  out
}

##' @export
##' @rdname summary.ism
summary.BlockedInfinitySparseMatrix <- function(object, ...) {
  out <- lapply(levels(object@groups),
                function(x) {
                  ism <- subset(object,
                                subset=object@rownames %in% names(object@groups[object@groups == x]),
                                select=object@colnames %in% names(object@groups[object@groups == x]))
                  s <- summary(ism)
                  s$matname <- deparse(substitute(object))
                  s$blockname <- x
                  return(s)
                })
  names(out) <- levels(object@groups)

  out$overall <- summary.InfinitySparseMatrix(object)

  out$matname <- list(matname = deparse(substitute(object)),
                       blocknames = levels(object@groups))
  out$overall$matname <- out$matname$matname

  class(out) <- "summary.BlockedInfinitySparseMatrix"
  return(out)
}

##' @export
##' @rdname summary.ism
summary.DenseMatrix <- function(object, ...) {
  mtreat <- apply(object, 1, function(x) any(is.finite(x)))
  mcontrol <- apply(object, 2, function(x) any(is.finite(x)))
  distances <- summary(as.vector(object))

  out <- internal.summary.helper(object, mtreat, mcontrol, distances)

  out$matname <- deparse(substitute(object))

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
  out$total$matchable <- Reduce("+", num_eligible_matches(x))
  out$total$unmatchable <- prod(d) - out$total$matchable

  out$matchable$treatment <- rownames(x)[matchabletxt]
  out$matchable$control <- colnames(x)[matchablectl]
  out$unmatchable$treatment <- rownames(x)[!matchabletxt]
  out$unmatchable$control <- colnames(x)[!matchablectl]

  out$distances <- distances
  return(out)
}

##' Print method for InfinitySparseMatrix, BlockedInfinitySparseMatrix
##' and DenseMatrix
##'
##' @param x Output from a summary call of an InfinitySparseMatrix,
##'   BlockedInfinitySparseMatrix or DenseMatrix (of class
##'   summary.infinitysparsematrix,
##'   summary.blockedinfinitysparsematrix or summary.densematrix
##'   respectively).
##' @param printAllBlocks If `x` is a BlockedInfinitySparseMatrix,
##'   should summaries of all blocks be printed alongside the overall
##'   summary?
##' @param ... Additional arguments passed to `print` calls.
##' @return `x`
##' @export
##' @name print.summary.ism
print.summary.InfinitySparseMatrix <- function(x, ...) {
  cat(paste("Membership:", x$total$treatment, "treatment,",
            x$total$control, "control\n"))
  cat(paste("Total eligible potential matches:", x$total$matchable,
            "\n"))
  cat(paste("Total ineligible potential matches:", x$total$unmatchable,
            "\n"))
  cat("\n")

  numunmatch <- sapply(x$unmatchable, length)
  for (i in 1:2) {
    if (numunmatch[i] > 0) {
      cat(paste0(numunmatch[i], " unmatchable ", names(numunmatch)[i],
                 " member", if(numunmatch[i] > 1) { "s" } , ":\n"))
      cat("\t")
      cat(paste(x$unmatchable[[i]][1:min(5, numunmatch[i])],
                collapse=", "))
      if (numunmatch[i] > 5) {
        cat(", ...\n")

        cat(paste0("See summary(", x$matname, ")",
                   if (!is.null(x$blockname)) {
                     paste0("$`", x$blockname, "`")
                   }, "$unmatchable$",
                   names(numunmatch)[i], " for a complete list."))
        }
      cat("\n\n")
    }
  }

  if (any(!is.na(x$distances))) {
    cat("Summary of distances:\n")
    print(x$distances, ...)
    cat("\n")
  }
  return(x)
}

##' @export
##' @rdname print.summary.ism
print.summary.BlockedInfinitySparseMatrix <- function(x, ..., printAllBlocks=FALSE) {

  cat("Summary across all blocks:\n")
  print(x$overall, ...)
  blockentries <- names(x) %in% x$matname$blocknames

  if (!printAllBlocks) {
    cat("Block structure:\n")
    blocksummary <- matrix(unlist(lapply(x[blockentries], "[", "total")),
                           byrow=TRUE, ncol=4)
    rownames(blocksummary) <- x$matname$blocknames
    colnames(blocksummary) <- c("#Treatment", "#Control", "Matchable",
                                "Unmatchable")
    print(blocksummary)

    cat("\n")

    cat(paste0("To see summaries for individual blocks,",
               " call for example summary(",
               x$matname$matname, ")$`",
               x$matname$blocknames[1], "`.\n"))
  } else {
    cat("Indiviual blocks:\n\n")
    print(x[blockentries], ...)
  }

  cat("\n")
  return(x)
}

##' @export
##' @rdname print.summary.ism
print.summary.DenseMatrix <- function(x, ...) {
  print.summary.InfinitySparseMatrix(x, ...)
  return(x)
}
