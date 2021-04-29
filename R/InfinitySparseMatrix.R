#' Indexing ISM's
#'
#' @param x An InfinitySparseMatrix object.
#' @param i Row indices.
#' @param j Col indices.
#' @param ... Other arguments
#' @param drop Ignored.
#' @return The subset.
#' @export
setMethod("[", "InfinitySparseMatrix",
          function(x, i, j =NULL, ..., drop = TRUE) {
            if ( "drop" %in% names(match.call())) {
              warning("'drop' argument ignored for InfinitySparseMatrix object subsetting ")
            }
            # Handles [X] cases
            if (nargs() < 3) {
              if (missing(i)) {
                return(x)
              }
              return(x@.Data[i, ...])
            } else {
              # At this point we have two arguments, but one could be
              # null (e.g. [X,] or [,X], as opposed to [X])

              # when missing, replace with NULL to be handled below
              if (missing(i)) i <- NULL
              if (missing(j)) j <- NULL

              makelogical <- function(index, rowcol) {
                switch(class(index),
                       "numeric" = (1:dim(x)[rowcol]) %in% index,
                       "integer" = (1:dim(x)[rowcol]) %in% index,
                       "character" = dimnames(x)[[rowcol]] %in% index,
                       "logical" = index,
                       "NULL" = rep(TRUE, dim(x)[rowcol]),
                       stop("Unrecognized class"))
              }

              subi <- makelogical(i, 1)
              subj <- makelogical(j, 2)

              subset(x, subset=subi, select=subj)
            }
          })
