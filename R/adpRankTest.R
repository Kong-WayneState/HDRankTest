#' @name adpRankTest
#' @title Adaptive Rank-Based Test
#' @description Kong, Ye, and Weng (2026): Adaptive rank-based inference for
#' high-dimensional two-sample problems.
#'
#' @param X1 A numeric \eqn{(n_1 \times p)} matrix of data from sample 1.
#' @param X2 A numeric \eqn{(n_2 \times p)} matrix of data from sample 2.
#' @param omega The true value of the relative effect.
#' @param W The lag window width. If \code{NA}, the elbow method is used to
#'   determine the value. If \code{"CV"}, the value is selected via
#'   cross-validation using the tuning parameters \code{step} and \code{fold}.
#' @param pow A vector of powers used to construct the adaptive test.
#' @param ... Arguments passed to the \code{find_win} function.
#' Default options \code{step = 1}, \code{cv.fold = 10},
#'   and \code{norm.type = "E"}
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
#' adpRankTest(X1, X2, W = "CV")
#'
#'
adpRankTest <- function(X1, X2, omega = 0.5, W = NA, pow = c(1:6, Inf), ...){
  p  <- ncol(X1)
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  pow <- sort(unique(c(pow, 2, Inf)))  # Include powers 2, and Inf

  # Rank transformation
  Y <- mapply(function(x1, x2){rank(c(x1, x2)) - c(rank(x1), rank(x2))},
    as.data.frame(X1), as.data.frame(X2))

  Y1 <- Y[1:n1, ] / n2
  Y2 <- Y[-(1:n1), ] / n1

  # Column means and variances
  # Y1.mean <- colMeans(Y1)
  Y2.mean <- colMeans(Y2)                 # Y1.mean + Y2.mean = 1
  Y1.var  <- apply(Y1, 2, var)
  Y2.var  <- apply(Y2, 2, var)

  D <- Y2.mean - omega                    # Difference under H0
  D.sq <- (D^2)/(Y1.var/n1 + Y2.var/n2)

  # ------------- Rank-based Supremum-Type Test -------------
  R.inf.stat <- max(D.sq)
  pval.inf <- 1 - exp(-exp(-(R.inf.stat - 2*log(p) + log(log(p)))/2)/sqrt(pi))

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

  # ------------- Rank-based Sum-of-Power-Type Test -------------
  pow.fin <- pow[is.finite(pow)]                       # Finite powers only
  R.pow <- sapply(pow.fin, function(gam){sum(D^gam)})
  R.pow.mean <- est_rsp_mean(Y1, Y2, pow = pow.fin)        # Estimate mean
  R.pow.cov <- est_rsp_cov(Y1, Y2, win = W, pow = pow.fin) # Estimate var-cov matrix
  R.pow.var <- diag(R.pow.cov)                         # Estimate var
  R.pow.stat <- (R.pow - R.pow.mean) / sqrt(R.pow.var)
  pval.pow <- ifelse(pow.fin %% 2 == 1,
                     2 * (1 - pnorm(abs(R.pow.stat))), # odd power: two-tails
                     1 - pnorm(R.pow.stat)             # even power: right-tail
                     )

  # ------------- Rank-based Adaptive Test -------------
  # Joint of odd powers
  odd.pow.id <- which(pow.fin %% 2 == 1)
  n.odd.pow <- length(odd.pow.id)
  if(n.odd.pow == 0){
    pval.odd <- NA
  }else if(n.odd.pow == 1){
    pval.odd <- pval.pow[odd.pow.id]
  }else{
    R.odd.max <- max(abs(R.pow.stat[odd.pow.id]))
    Cov.odd <- cov2cor(R.pow.cov[odd.pow.id, odd.pow.id])
    pval.odd <- ifelse(is.na(R.odd.max), NA,
                       1 - mvtnorm::pmvnorm(lower = rep(-R.odd.max, n.odd.pow),
                                            upper = rep(R.odd.max, n.odd.pow),
                                            mean = rep(0, n.odd.pow),
                                            sigma = Cov.odd))
  }

  #Joint of even powers
  even.pow.id <- which(pow.fin %% 2 == 0)
  n.even.pow <- length(even.pow.id)
  if(n.even.pow == 0){
    pval.even <- NA
  }else if(n.even.pow == 1){
    pval.even <- pval.pow[even.pow.id]
  }else{
    R.even.max <- max(R.pow.stat[even.pow.id])
    Cov.even <- cov2cor(R.pow.cov[even.pow.id, even.pow.id])
    pval.even <- ifelse(is.na(R.even.max), NA,
                      1 - mvtnorm::pmvnorm(lower = rep(-Inf, n.even.pow),
                                           upper = rep(R.even.max, n.even.pow),
                                           mean = rep(0, n.even.pow),
                                           sigma = Cov.even))
  }

  pval.min <- min(c(pval.odd, pval.even, pval.inf), na.rm = TRUE)
  nadp <- 3 - sum(is.na(c(pval.odd, pval.even, pval.inf)))
  # Adaptive P-value via Tippett’s method
  pval.adp <- 1 - (1 -  pval.min)^nadp

  # ------------- Output -------------
  stat = c(R.pow.stat, R.inf.stat)
  pval = c(pval.pow, pval.inf, pval.adp)

  names(stat) = c(paste("RSPU_", pow, sep = ""))
  names(pval)= c(paste("RSPU_", pow, sep = ""), "aRSPU")

  return(list(Statistic = stat, p.value = pval, window = W))
}


#' Center Estimation
#' @keywords internal
#' @noRd
est_rsp_mean <- function(Y1, Y2, pow){
  n1 <- nrow(Y1)
  n2 <- nrow(Y2)
  Y1.var <- apply(Y1, 2, var)
  Y2.var <- apply(Y2, 2, var)

  sk <- Y1.var / n1 + Y2.var / n2

  if (any((pow %% 2 == 1) & (pow > 1))){
    m1k <- apply(Y1, 2, function(y){mean((y - mean(y))^3)}) # m1
    m2k <- apply(Y2, 2, function(y){mean((y - mean(y))^3)}) # m2
    m12k <- mapply(function(y1, y2){
      mean(outer(y1, y2, function(r1, r2){
        ((r2 > r1) + (r2 == r1) / 2 - mean(y2)) * (r2 - mean(y2))*(r1 -mean(y1))
      }))},
      as.data.frame(Y1), as.data.frame(Y2)) # m12
    mk <- m2k / (6 * n2^2) - m1k / (6*n1^2) - m12k / (n1 * n2)
  }else {mk <- NULL}

  n.pow <- length(pow)
  R.mean <- numeric(n.pow)
  for(i in 1:n.pow){
    if(pow[i] == 1){
      R.mean[i] = 0
    }else if(pow[i] %% 2 == 1){
      a <- (pow[i] - 3) / 2
      R.mean[i] <- factorial(pow[i]) / (factorial(a) * (2^a)) * sum(mk * sk^a)
    }else if(pow[i] %% 2 == 0){
      a <- pow[i] / 2
      R.mean[i] <- factorial(pow[i]) / (factorial(a) * (2^a)) * sum(sk^a)
    }
  }
  return(R.mean)
}

#' Variance-Covariance Matrix Estimation
#' @keywords internal
#' @noRd
est_rsp_cov <- function(Y1, Y2, win = 0, pow = 1:6) {
  n1 <- nrow(Y1)
  n2 <- nrow(Y2)
  S1 = cov(Y1) / n1
  S2 = cov(Y2) / n2

  idx = abs(row(S1) - col(S1)) > win # truncate
  S1[idx] <- 0
  S2[idx] <- 0

  sk <- diag(S1) + diag(S2)

  # ---------- unimean ----------
  unimean <- function(gam) {
    if (gam %% 2 == 1) return(0)
    factorial(gam)/factorial(gam/2) * (sum(sk^(gam/2)) / (2^(gam/2)))
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
  R.cov <- matrix(0, n.pow, n.pow)

  for(i in 1 : n.pow){
    for(j in 1 : i){
      P1 <- unimean(pow[i] + pow[j])
      P2 <- mu_0[i] * mu_0[j]
      P3 <- mixmean(pow[i], pow[j])
      val <- P1 - P2 + P3
      R.cov[i, j] <- val
      R.cov[j, i] <- val
    }}
  return(R.cov)
}


