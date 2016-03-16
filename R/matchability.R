##' Given a distance matrix which has potentially been calipered,
##' returns which observations are matchable and unmatchable
##'
##' @param distance ISM, BISM or Matrix
##' @return List of lists, $matchable$control, $matchable$treatment,
##'   $unmatchable$control andd $unmatchable$treatment
##' @export
matchability <- function(distance) {
  UseMethod("matchability")
}

matchability.InfinitySparseMatrix <- function(distance) {
  out <- list()
  d <- dim(distance)

  matchabletreatment <- 1:d[1] %in% sort(unique(distance@rows))
  matchablecontrol   <- 1:d[2] %in% sort(unique(distance@cols))
  out$matchable$treatment <- distance@rownames[matchabletreatment]
  out$matchable$control <- distance@colnames[matchablecontrol]
  out$unmatchable$treatment <- distance@rownames[!matchabletreatment]
  out$unmatchable$control <- distance@colnames[!matchablecontrol]
  out
}

matchability.BlockedInfinitySparseMatrix <- function(distance) {
}

matchability.matrix <- function(distance) {
  out <- list()

  matchabletreatment <- apply(distance, 1, function(x) any(x != Inf))
  matchablecontrol   <- apply(distance, 2, function(x) any(x != Inf))
  out$matchable$treatment <- rownames(distance)[matchabletreatment]
  out$matchable$control <- colnames(distance)[matchablecontrol]
  out$unmatchable$treatment <- rownames(distance)[!matchabletreatment]
  out$unmatchable$control <- colnames(distance)[!matchablecontrol]
  out
}
