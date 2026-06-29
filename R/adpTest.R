# Two-Sample Inference for High-Dimensional Mean Vectors with Unequal Covariance Matrices

#' @name adpTest
#' @title Adaptive Test
#' @description XU, Lin, Wei and Pan (2016, Biometrika): An adaptive two-sample
#' test for high-dimensional means.
#'
#' @param X1 A numeric \eqn{(n_1 \times p)} matrix of data from sample 1.
#' @param X2 A numeric \eqn{(n_2 \times p)} matrix of data from sample 2.
#' @param W The lag window width. If \code{NA}, the elbow method is used to
#'   determine the value. If \code{"CV"}, the value is selected via
#'   cross-validation using the tuning parameters \code{step} and \code{fold}.
#' @param pow A vector of powers used to construct the adaptive test.
#' @param ... Arguments passed to the \code{find_win} function.
#' Default options \code{step = 1}, \code{cv.fold = 10},
#'   and \code{norm.type = "E"}
#'
#' @return values of test statistic, p-value and window width.
#'
#' @export
#'
#' @importFrom stats var pnorm cov2cor acf
#' @importFrom pathviewr find_curve_elbow
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
#' adpTest(X1, X2, W = "CV")
#'
#'
adpTest <- function(X1, X2, W = NA, pow = c(1:6, Inf), ...){
  p  <- ncol(X1)
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  pow <- sort(unique(c(pow, 2, Inf)))  # Include powers 2, and Inf

  X1.mean <- colMeans(X1)
  X2.mean <- colMeans(X2)
  X1.var  <- apply(X1, 2, var)
  X2.var  <- apply(X2, 2, var)

  D <- X2.mean - X1.mean                # Difference under H0
  D.sq <- (D^2)/(X1.var/n1 + X2.var/n2)

  # ------------- Supremum-Type Test -------------
  L.inf.stat <- max(D.sq)
  pval.inf <- 1 - exp(-exp(-(L.inf.stat - 2*log(p) + log(log(p)))/2)/sqrt(pi))

  # Find the lag window width if needed
  if (is.na(W)) {
    tacf <- as.vector(
      stats::acf(D.sq, lag.max = p - 1, plot = FALSE, type = "covariance")$acf
    )
    W <- pathviewr::find_curve_elbow(data.frame(index = 1:p, acf = tacf))  - 1
  }else if (W == "CV") {
    X_scale <- rbind(scale(X1, center = TRUE, scale = FALSE),
                     scale(X2, center = TRUE, scale = FALSE))
    W <- find_win(X_scale, ...)
  }

  # ------------- Sum-of-Power-Type Test -------------
  pow.fin <- pow[is.finite(pow)]                       # Finite powers only
  L.pow <- sapply(pow.fin, function(gam){sum(D^gam)})
  L.pow.mean <- est_sp_mean(X1, X2, pow = pow.fin)        # Estimate mean
  L.pow.cov <- est_sp_cov(X1, X2, win = W, pow = pow.fin) # Estimate var-cov matrix
  L.pow.var <- diag(L.pow.cov)                         # Estimate var
  L.pow.stat <- (L.pow - L.pow.mean) / sqrt(L.pow.var)
  pval.pow <- ifelse(pow.fin %% 2 == 1,
                     2 * (1 - pnorm(abs(L.pow.stat))), # odd power: two-tails
                     1 - pnorm(L.pow.stat)             # even power: right-tail
  )

  # ------------- Adaptive Test -------------
  # Joint of odd powers
  odd.pow.id <- which(pow.fin %% 2 == 1)
  n.odd.pow <- length(odd.pow.id)
  if(n.odd.pow == 0){
    pval.odd <- NA
  }else if(n.odd.pow == 1){
    pval.odd <- pval.pow[odd.pow.id]
  }else{
  L.odd.max <- max(abs(L.pow.stat[odd.pow.id]))
  Cov.odd <- cov2cor(L.pow.cov[odd.pow.id, odd.pow.id])
  pval.odd <- ifelse(is.na(L.odd.max), NA,
                     1 - mvtnorm::pmvnorm(lower = rep(-L.odd.max, n.odd.pow),
                                          upper = rep(L.odd.max, n.odd.pow),
                                          mean = rep(0, n.odd.pow),
                                          sigma = Cov.odd))
  }
  # Joint of even powers
  even.pow.id <- which(pow.fin %% 2 == 0)
  n.even.pow <- length(even.pow.id)
  if(n.even.pow == 0){
    pval.even <- NA
  }else if(n.even.pow == 1){
    pval.even <- pval.pow[even.pow.id]
  }else{
  L.even.max <- max(L.pow.stat[even.pow.id])
  Cov.even <- cov2cor(L.pow.cov[even.pow.id, even.pow.id])
  pval.even <- ifelse(is.na(L.even.max), NA,
                      1 - mvtnorm::pmvnorm(lower = rep(-Inf, n.even.pow),
                                           upper = rep(L.even.max, n.even.pow),
                                           mean = rep(0, n.even.pow),
                                           sigma = Cov.even))
  }

  pval.min <- min(c(pval.odd, pval.even, pval.inf), na.rm = TRUE)

  nadp <- 3 - sum(is.na(c(pval.odd, pval.even, pval.inf)))

  # Adaptive P-value via Tippett’s method
  pval.adp <- 1 - (1 -  pval.min)^nadp

  # ------------- Output -------------
  stat = c(L.pow.stat, L.inf.stat)
  pval = c(pval.pow, pval.inf, pval.adp)

  names(stat) = c(paste("SPU_", pow, sep = ""))
  names(pval)= c(paste("SPU_", pow, sep = ""), "aSPU")

  return(list(Statistic = stat, p.value = pval, window = W))

}

#' Center Estimation
#' @keywords internal
#' @noRd
est_sp_mean <- function(X1, X2, pow){
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  X1.var <- apply(X1, 2, var)
  X2.var <- apply(X2, 2, var)

  sk <- X1.var / n1 + X2.var / n2

  n.pow <- length(pow)
  L.mean <- numeric(n.pow)
  for(i in 1:n.pow){
    if(pow[i] %% 2 == 1){
      L.mean[i] = 0
    }else{
      a <- pow[i] / 2
      L.mean[i] <- factorial(pow[i]) / (factorial(a) * (2^a)) * sum(sk^a)
    }
  }
  return(L.mean)
}


#' Variance-Covariance Matrix Estimation
#' @keywords internal
#' @noRd
est_sp_cov <- function(X1, X2, win = 0, pow = 1:6){
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  S1 = cov(X1) / n1
  S2 = cov(X2) / n2

  idx = abs(row(S1) - col(S1)) > win # truncate
  S1[idx] <- 0
  S2[idx] <- 0

  sk <- diag(S1) + diag(S2)

  # ---------- unimean ----------
  unimean <- function(gam) {
    ifelse (gam %% 2 == 1, 0,
    factorial(gam)/factorial(gam/2) * (sum(sk^(gam/2)) / (2^(gam/2))))
  }

  mu_0 <- sapply(pow, unimean)

  # ---------- mixmean ----------
  mixmean <- function(gam, eta){
    if ((gam + eta) %% 2 == 1) return(0)
    d <- gam %% 2
    temp <- 0
    while (d <= min(gam, eta)) {
      mat <- outer(sk, sk, function(a, b)
        a^((gam - d) / 2) * b^((eta - d) / 2)) * (S1 + S2)^d
      diag(mat) <- 0
      dem <- factorial((gam - d) / 2) *
        factorial((eta - d) / 2) * factorial(d) * 2^((gam + eta) / 2 - d)
      temp <- temp + sum(mat) / dem
      d <- d + 2
    }
    return(factorial(gam) * factorial(eta) * temp)
  }

  n.pow <- length(pow)
  L.cov <- matrix(0, n.pow, n.pow)

  for(i in 1 : n.pow){
    for(j in 1 : i){
      P1 <- unimean(pow[i] + pow[j])
      P2 <- mu_0[i] * mu_0[j]
      P3 <- mixmean(pow[i], pow[j])
      val <- P1 - P2 + P3
      L.cov[i, j] <- val
      L.cov[j, i] <- val
    }}
  return(L.cov)

  }



