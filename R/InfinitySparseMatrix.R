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
##' @param distanceSummary Default TRUE. Should a summary of minimum
##'   distance per treatment member be calculated? May be slow on
##'   larger data sets.
##' @param printAllBlocks If `object` is a
##'   BlockedInfinitySparseMatrix, should summaries of all blocks be
##'   printed alongside the overall summary? Default FALSE.
##' @param blockStructure If `object` is a BlockedInfinitySparseMatrix
##'   and `printAllBlocks` is false, print a quick summary of each
##'   individual block. Default TRUE. If the number of blocks is high,
##'   consider suppressing this.
##' @return List of lists, $matchable$control, $matchable$treatment,
##'   $unmatchable$control andd $unmatchable$treatment
##' @export
##' @name summary.ism
summary.InfinitySparseMatrix <- function(object, ..., distanceSummary=TRUE) {

  finitedata <- is.finite(object@.Data)
  mtreat <- 1:dim(object)[1] %in% sort(unique(object@rows[finitedata]))
  mcontrol  <- 1:dim(object)[2] %in% sort(unique(object@cols[finitedata]))

  if (distanceSummary & length(object@.Data[finitedata])) {
    distances <- summary(tapply(object@.Data[finitedata],
                                object@rows[finitedata],
                                min))
  } else {
    distances <- NULL
  }

  out <- internal.summary.helper(object, mtreat, mcontrol, distances)
  attr(out, "ismname") <- deparse(substitute(object))

  class(out) <- "summary.InfinitySparseMatrix"
  out
}

##' @export
##' @rdname summary.ism
summary.BlockedInfinitySparseMatrix <- function(object, ...,
                                                distanceSummary=TRUE,
                                                printAllBlocks=FALSE,
                                                blockStructure=TRUE) {

  ismname <- deparse(substitute(object))

  out <- lapply(levels(object@groups),
                function(x) {
                  ism <- subset(object,
                                subset=object@rownames %in% names(object@groups[object@groups == x]),
                                select=object@colnames %in% names(object@groups[object@groups == x]))
                  s <- summary(ism, ..., distanceSummary=distanceSummary)
                  attr(s, "ismname") <- ismname
                  attr(s, "blockname") <- x
                  return(s)
                })
  names(out) <- levels(object@groups)

  out$overall <- summary.InfinitySparseMatrix(object, ...,
                                              distanceSummary=distanceSummary)

  attr(out, "ismname") <- ismname
  attr(out, "blocknames") <- levels(object@groups)

  attr(out$overall, "ismname") <- attr(out, "ismname")

  attr(out, "printAllBlocks") <- printAllBlocks
  attr(out, "blockStructure") <- blockStructure

  class(out) <- "summary.BlockedInfinitySparseMatrix"
  return(out)
}

##' @export
##' @rdname summary.ism
summary.DenseMatrix <- function(object, ..., distanceSummary=TRUE) {
  mtreat <- apply(object, 1, function(x) any(is.finite(x)))
  mcontrol <- apply(object, 2, function(x) any(is.finite(x)))
  if (distanceSummary & length(object@.Data[is.finite(object@.Data)])) {
    distances <- summary(apply(object, 1, min))
  } else {
    distances <- NULL
  }

  out <- internal.summary.helper(object, mtreat, mcontrol, distances)

  attr(out, "ismname") <- deparse(substitute(object))

  class(out) <- "summary.DenseMatrix"
  out
}

internal.summary.helper <- function(x,
                                    matchabletxt,
                                    matchablectl,
                                    distances=NULL) {
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

        cat(paste0("See summary(", attr(x, "ismname"), ")",
                   if (!is.null(attr(x, "blockname"))) {
                     paste0("$`", attr(x, "blockname"), "`")
                   }, "$unmatchable$",
                   names(numunmatch)[i], " for a complete list."))
        }
      cat("\n\n")
    }
  }

  if (!is.null(x$distances) && any(!is.na(x$distances))) {
    cat("Summary of minimum matchable distance per treatment member:\n")
    print(x$distances, ...)
    cat("\n")
  }
  return(invisible(x))
}

##' @export
##' @rdname print.summary.ism
print.summary.BlockedInfinitySparseMatrix <- function(x, ...) {

  cat("Summary across all blocks:\n")
  print(x$overall, ...)
  blockentries <- names(x) %in% attr(x, "blocknames")

  if (!attr(x, "printAllBlocks")) {
    if (attr(x, "blockStructure")) {
      cat("Block structure:\n")
      blocksummary <- matrix(unlist(lapply(x[blockentries], "[", "total")),
                             byrow=TRUE, ncol=4)
      blocksummary <- cbind(sapply(sapply(sapply(x[blockentries], "[", "matchable"), "[", "treatment"), length),
                            sapply(sapply(sapply(x[blockentries], "[", "matchable"), "[", "control"), length),
                            sapply(sapply(sapply(x[blockentries], "[", "unmatchable"), "[", "treatment"), length),
                            sapply(sapply(sapply(x[blockentries], "[", "unmatchable"), "[", "control"), length))
      rownames(blocksummary) <- paste0("`",attr(x, "blocknames"),"`")
      colnames(blocksummary) <- c("Matchable Txt",
                                  "Matchable Ctl",
                                  "Unmatchable Txt",
                                  "Unmatchable Ctl")
      print(blocksummary)

      cat("\n")
    }

    cat(paste0("To see summaries for individual blocks,",
               " call for example summary(",
               attr(x, "ismname"), ")$`",
               attr(x, "blocknames")[1], "`.\n"))
  } else {
    cat("Indiviual blocks:\n\n")
    for (i in attr(x, "blocknames")) {
      cat(paste0("`",i,"`\n"))
      print(x[[i]])
    }
  }

  cat("\n")
  return(invisible(x))
}

##' @export
##' @rdname print.summary.ism
print.summary.DenseMatrix <- function(x, ...) {
  print.summary.InfinitySparseMatrix(x, ...)
  return(invisible(x))
}
