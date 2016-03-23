#' Indexing ISM's
#'
#' @param x An InfinitySparseMatrix object.
#' @param i Row indices.
#' @param j Col indices. If NULL, then indexing x@.Data occurs.
#' @param ... Other arguments
#' @param drop Ignored.
#' @return The subset.
#' @export
setMethod("[", "InfinitySparseMatrix",
          function(x, i, j=NULL, ..., drop=TRUE) {
            if (is.null(j)) {
              return(x@.Data[i, drop=drop, ...])
            } else {
              makelogical <- function(index, rowcol) {
                switch(class(index),
                      "numeric" = (1:dim(x)[rowcol]) %in% index,
                      "integer" = (1:dim(x)[rowcol]) %in% index,
                      "character" = dimnames(x)[[rowcol]] %in% index,
                      "logical" = index,
                      stop("Unrecognized class"))
              }

              subi <- makelogical(i, 1)
              subj <- makelogical(j, 2)

              subset(x, subset=subi, select=subj)
            }
          }
          )
