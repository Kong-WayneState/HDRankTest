#' @name ovlRankTest
#' @title Overall Rank-Based Test
#' @description High dimensional Rank-based inference
#'              by Kong and Harrar (2020)
#'
#' @param X1 A numeric \eqn{(n_1 \times p)} matrix of data from sample 1.
#' @param X2 A numeric \eqn{(n_2 \times p)} matrix of data from sample 2.
#'
#' @return values of test statistic and p-value.
#'
#' @export
#'
#' @importFrom stats cov pnorm
#'
#' @examples
#' # generating 2-sample data
#' set.seed(123)
#' p = 300
#' n1 = 60
#' n2 = 80
#' X1 = matrix(rnorm(n1*p), nrow = n1, ncol = p)
#' X2 = matrix(rnorm(n2*p), nrow = n2, ncol = p)
#'
#' ovlRankTest(X1, X2)
#'
#'
ovlRankTest <- function(X1, X2){
  p  <- ncol(X1)
  n1 <- nrow(X1)
  n2 <- nrow(X2)

  # Rank transformation
  Y <- matrix((rank(t(rbind(X1, X2))) - 1/2) / (p * (n1 + n2)), nrow = p)
  Y1 <- t(Y[ , 1 : n1])
  Y2 <- t(Y[ , (n1 + 1) : (n1 + n2)])
  S1 <- cov(Y1)
  S2 <- cov(Y2)
  TS <- sum((colMeans(Y1) - colMeans(Y2))^2) -
                                       (sum(diag(S1)) / n1 + sum(diag(S2)) / n2)

  c12 <- sum(S1 * S2)
  A1 <- Y1 %*% t(Y1)
  diag(A1) <- 0
  A2 <- Y2 %*% t(Y2)
  diag(A2) <- 0
  c1 <- (sum(A1^2) * (n1 - 1) * (n1 - 2) - 2 * sum((rowSums(A1))^2) * (n1 - 1) +
           (sum(A1))^2) / (n1 * (n1 - 1) * (n1 - 2) * (n1 - 3))
  c2 <- (sum(A2^2) * (n2 - 1) * (n2 - 2) - 2 * sum((rowSums(A2))^2) * (n2 - 1) +
           (sum(A2))^2) / (n2 * (n2 - 1) * (n2 - 2) * (n2 - 3))

  var.est <- 2 * c1 / n1 / (n1 - 1) + 2 * c2 / n2 / (n2 - 1) + 4 * c12 / n1 / n2

  stat <- TS / sqrt(var.est)
  pval <- pnorm(stat, lower.tail = FALSE)

  return(list(Statistic = stat, p.value = pval))
}
