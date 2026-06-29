# Two-Sample Inference for High-Dimensional Mean Vectors with Unequal Covariance Matrices

#' @name cmpTest
#' @title Composite Test
#' @description Gregory, Carroll, Baladandayuthapani and Lahiri (2015, JASA):
#' A two-sample test for equality of means in high dimension;
#' and Zhang and Wang (2021, CSDA): A more powerful test of equality of
#' high-dimensional two-sample means.
#'
#' @param X1 A numeric \eqn{(n_1 \times p)} matrix of data from sample 1.
#' @param X2 A numeric \eqn{(n_2 \times p)} matrix of data from sample 2.
#' @param W The lag window width. If \code{NA}, defaults to \code{ceiling(p^(3/8))}
#' @param order The order of the centering correction.
#'    Available options are 0, 1, and 2.
#' @param method The scale estimation method.
#'    Available options are \code{"mpt"} and \code{"gct"}.
#' @param weight.fun window weight function, \code{"parz"} for Parzen weight and
#'     \code{"trapez"} for trapezoid weight.
#'
#' @return values of test statistic, p-value and window width.
#'
#' @export
#'
#' @importFrom stats var cov pnorm acf
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
#' cmpTest(X1, X2)
#' cmpTest(X1, X2, method = "gct")
#'
#'

# Note that both GCT.test() in the highD2pop package and
# zwl_test() in the highDmean package use a two-sided p-value,
# whereas the present implementation uses a one-sided p-value.
cmpTest = function(X1, X2, W = NA, order = 0, method = "mpt",
                   weight.fun = "parz"){
  p <- ncol(X1)
  n1 <- nrow(X1)
  n2 <- nrow(X2)

  X1.mean <- colMeans(X1)
  X2.mean <- colMeans(X2)
  X1.var <- apply(X1, 2, var)
  X2.var <- apply(X2, 2, var)

  D.sq <- (X1.mean - X2.mean)^2/(X1.var / n1 + X2.var / n2)
  TS <- mean(D.sq)

  # Estimate mean
  center.est <- mean(mapply(est_mean_1d, as.data.frame(X1), as.data.frame(X2),
                            MoreArgs = list(ntoorderminus = order)))
  if (is.na(W)){
    W <- ceiling(p^(3/8))
  }

  if (method == "mpt"){
    X1.cov <- cov(X1)
    X2.cov <- cov(X2)
    lambda1 <- n1 / (n1 + n2)
    lambda2 <- n2 / (n1 + n2)
    numerator <- 2 * (lambda2 * X1.cov + lambda1 *  X2.cov)^2
    a <-  lambda2 * diag(X1.cov)
    b <-  lambda1 * diag(X2.cov)
    denom <-  (a + b) %*% t(a + b)
    idx = outer(1:p, 1:p, function(x, y) {abs(x - y) <= W})
    var.est = sum(numerator[idx] / denom[idx]) / p^2
  } else{
    gamma <- as.vector(
      stats::acf(D.sq, lag.max = W, plot = FALSE, type = "covariance")$acf
    )
    var.est <- sum(c(1, 2*get(weight.fun)(W + 1)[-(W+1)]) * gamma/(p : (p - W)))
  }

  stat = (TS - center.est)/sqrt(var.est)
  pval <- pnorm(stat, lower.tail = FALSE) # even power: right-tail

  return(list(Statistic = stat, p.value = pval, window = W))
}


#' Center Estimation (marginal)
#' @keywords internal
#' @noRd
est_mean_1d <- function(x, y, ntoorderminus = 2){
  n <- length(x)
  m <- length(y)
  if (ntoorderminus == 0) {
    return(1)
  } else if (ntoorderminus == 1) {
    sig.sq.x.hat <- mean(x^2) - mean(x)^2
    sig.sq.y.hat <- mean(y^2) - mean(y)^2
    tau.sq.hat <- sig.sq.x.hat + (n/m) * sig.sq.y.hat
    mu.3.hat <- mean((x - mean(x))^3)
    eta.3.hat <- mean((y - mean(y))^3)
    a1 <- tau.sq.hat^(-1) * (sig.sq.x.hat + (n/m)^2 * sig.sq.y.hat)
    a2 <- tau.sq.hat^(-3) * 2 * (mu.3.hat + (n/m)^2 * eta.3.hat)^2
    c <- 1 + n^(-1) * (a1 + a2)
    return(c)
  } else if (ntoorderminus == 2) {
    sig.sq.x.hat <- mean(x^2) - mean(x)^2
    sig.sq.y.hat <- mean(y^2) - mean(y)^2
    tau.sq.hat <- sig.sq.x.hat + (n/m) * sig.sq.y.hat
    mu.3.hat <- mean((x - mean(x))^3)
    eta.3.hat <- mean((y - mean(y))^3)
    mu.4.hat <- mean((x - mean(x))^4)
    eta.4.hat <- mean((y - mean(y))^4)
    mu.5.hat <- mean((x - mean(x))^5)
    eta.5.hat <- mean((y - mean(y))^5)
    a1 <- tau.sq.hat^(-1) * (sig.sq.x.hat + (n/m)^2 * sig.sq.y.hat)
    a2 <- tau.sq.hat^(-3) * 2 * (mu.3.hat + (n/m)^2 * eta.3.hat)^2
    b1 <- tau.sq.hat^(-2) * ((sig.sq.x.hat + (n/m)^2 * sig.sq.y.hat) -
                               ((mu.4.hat - 3 * sig.sq.x.hat^2) + (n/m)^4 *
                                  (eta.4.hat - 3 * sig.sq.y.hat^2)))
    b2 <- tau.sq.hat^(-3) * ((sig.sq.x.hat + (n/m)^2 * sig.sq.y.hat) *
                               ((mu.4.hat - sig.sq.x.hat^2) + (n/m)^3 *
                                  (eta.4.hat - sig.sq.y.hat^2)) -
                               4 * (mu.3.hat + (n/m)^2 * eta.3.hat) *
                               (mu.3.hat + (n/m)^3 * eta.3.hat) -
                               2 * (mu.3.hat^2 + (n/m)^5 * eta.3.hat^2))
    b3 <- tau.sq.hat^(-4) * (6 * (sig.sq.x.hat + (n/m)^2 * sig.sq.y.hat) *
                               (mu.3.hat + (n/m)^2 * eta.3.hat)^2 -
                               6 * (mu.3.hat + (n/m)^2 * eta.3.hat) *
                               (mu.5.hat - 2 * mu.3.hat * sig.sq.x.hat +
                                  (n/m)^4 * (eta.5.hat -
                                               2 * eta.3.hat * sig.sq.y.hat)) -
                               3 * ((mu.4.hat - sig.sq.x.hat^2) +
                                      (n/m)^3 * (eta.4.hat - sig.sq.y.hat^2))^2)
    b4 <- tau.sq.hat^(-5) * (3 * (sig.sq.x.hat + (n/m) *  sig.sq.y.hat) *
                               ((mu.4.hat - sig.sq.x.hat^2) +
                                  (n/m)^3 * (eta.4.hat - sig.sq.y.hat^2))^2 +
                               12 * (mu.3.hat + (n/m)^2 * eta.3.hat)^2 *
                               ((mu.4.hat - sig.sq.x.hat^2) +
                                  (n/m)^3 * (eta.4.hat - sig.sq.y.hat^2)))
    c <- 1 + n^(-1) * (a1 + a2) + n^(-2) * (b1 + b2 + b3 + b4)
    return(c)
  }
}



#' Input a positive integer r and return a list of Parzen weights
#' @keywords internal
#' @noRd
parz <- function(r) {
  sapply(1:r, function(i){ifelse(i < (r/2),
                                 1 - 6 * (i / r)^2 + 6 * (i / r)^3,
                                 2 * (1 - (i / r))^3)})
}



#' Input a positive integer r and return a list of trapezoid weights
#' @keywords internal
#' @noRd
trapez <- function(r){
  sapply(1:r, function(i){ifelse(i <= ceiling(r/2), 1,
                                 1 - (i-ceiling(r/2))/(r-ceiling(r/2)))})
  # in the paper, they used floor()
  # in the code, they used ceiling()
}


