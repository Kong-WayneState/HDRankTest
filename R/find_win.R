#' @name find_win
#' @title Find the Lag Window Width
#' @description Select the lag window width via cross-validation.
#' Candidates lag window width are defined as \code{0, step, 2*step, ..., p},
#' where \code{p} is the number of variables. Cross-validation is performed by
#' splitting the data into \code{cv.fold} parts. For each fold, the covariance
#' matrix is computed on the training data, banded according to the candidate
#' lag window width, and compared to the covariance matrix of the test data
#' using the specified norm.
#'
#' @param sam A numeric \eqn{(n \times p)} matrix of scaled observations.
#' @param step A positive integer giving the spacing of candidate bandwidths.
#'   Default is 5.
#' @param cv.fold An integer giving the number of cross-validation folds.
#'   Default is 10.
#' @param norm.type A character string specifying the matrix norm used in the
#'   loss calculation (e.g., "F" for Frobenius norm).
#'
#' @return An integer giving the selected lag window width.
#' @export
#'
#' @importFrom stats cov
#'
#'
find_win <- function(sam, step = 1, cv.fold = 10, norm.type = "E") {
  n <- nrow(sam)
  p <- ncol(sam)
  bandwidth <- seq(0, p, by = step)
  n.bandwidth <- length(bandwidth)
  if (n.bandwidth == 1)
    return(bandwidth)

  diff.norm <- matrix(0, cv.fold, n.bandwidth)
  fold.size <- floor(n / cv.fold)
  sam.idx <- sample.int(n)

  for (i in 1:cv.fold) {
    if (i < cv.fold) {
      temp.idx <- sam.idx[((i - 1) * fold.size + 1):(i * fold.size)]
    } else {
      temp.idx <- sam.idx[((i - 1) * fold.size + 1):n]
    }

    sam.train <- sam[-temp.idx, ]
    sam.test  <- sam[temp.idx, ]

    sam.train.cov <- cov(sam.train)
    sam.test.cov  <- cov(sam.test)

    for (j in seq_len(n.bandwidth)) {
      sam.train.cov.band <- sam.train.cov
      sam.train.cov.band[
        abs(row(sam.train.cov.band)-col(sam.train.cov.band)) > bandwidth[j]
      ] <- 0

      diff.norm[i, j] <-
        norm(sam.train.cov.band - sam.test.cov,
             type = norm.type)
    }
  }

  return(bandwidth[which.min(colMeans(diff.norm))])
}
