#' @name rData
#' @title Generate Random Data
#' @description Generates simulated data from an ARMA model, a multivariate
#' normal distribution, or a multivariate t distribution.
#'
#' @param p An integer specifying the number of variables (dimension).
#' @param n An integer specifying the number of observations (sample size).
#' @param mu A numeric vector of the means.
#' @param type A character string indicating the data generation type
#' (e.g., \code{"arma.norm"}, \code{"arma.t"}, \code{"arma.cauchy"}, \code{"arma.gamma"},
#' \code{"ld.norm"}, \code{"ld.t"}, \code{"ld.cauchy"}, \code{"ld.gamma"},
#' \code{"cov.norm"}, \code{"cov.t"}, \code{"cov.cauchy"}).\describe{
#' \item{"normal"}{Normal(mean = 0, sd = 1);}
#' \item{"t"}{scaled t(df = 3);}
#' \item{"cauchy"}{Cauchy(location = 0, scale = 0.1);}
#' \item{"gamma"}{centered Gamma(shape = 4, scale = 0.5).}
#' }
#' @param arma A list containing ARMA model parameters, such as
#' IND: \code{list(order=c(0,0,0))};
#' WD: \code{list(ar = c(.4, -.1), ma = c(0.2, 0.3))};
#' SD: \code{list(ar = 0.9)}.
#' @param H A numeric value representing the Hurst exponent for long-range dependence.
#' Default for LD: 0.7.
#' @param Sigma A \code{p} by \code{p} covariance matrix. Default is the identity matrix.
#'
#' @return A matrix of dimensions \code{n} by \code{p} containing the generated random data.
#' @export
#'
#' @importFrom stats arima.sim rnorm rt rgamma rcauchy
#' @importFrom mvtnorm rmvnorm rmvt
#'
#'
#' @examples
#' set.seed(123)
#' X1 <- rData(p = 300, n = 60)
#' X2 <- rData(p = 300, n = 80)
#'

rData <- function(p = 300, n = 60, mu = rep(0, p), type = "arma.norm",
                     arma=list(order=c(0, 0, 0)), H=0.7, Sigma = diag(p)){
  if (type == "arma.norm"){
    X <- rarma_normal(n=n, mu=mu, arma=arma)  #ARMA(2,2) N(0,1)
  }else if (type == "arma.t"){
    X <- rarma_t(n=n, mu=mu, arma=arma)   #ARMA(2,2) scaled t_3
  }else if (type == "arma.gamma"){
    X <- rarma_gamma(n=n, mu=mu, arma=arma)  #ARMA(2,2) centered Gamma(4,0.2)
  }else if (type == "arma.cauchy"){
    X <- rarma_cauchy(n=n, mu=mu, arma=arma)  #ARMA(2,2) Cauchy(0,0.1)
  }else  if (type == "ld.norm"){
    X <- rld_normal(n=n, mu=mu, H=H)  # long range dep N(0,1)
  }else if (type == "ld.t"){
    X <- rld_t(n=n, mu=mu, H=H)  # long range dep t_3
  }else if (type == "ld.gamma"){
    X <- rld_gamma(n=n, mu=mu, H=H)  # long range dep centered Gamma(4,2)
  }else if (type == "ld.cauchy"){
    X <- rld_cauchy(n=n, mu=mu, H=H)  # long range dep Cauchy(0,0.1)
  } else if (type == "cov.norm"){
    X <- mvtnorm::rmvnorm(n = n, mean = mu, sigma = Sigma)
  }else if (type == "cov.t"){
    X <- mvtnorm::rmvt(n = n, sigma = Sigma, df = 3, delta = mu)
  }else if (type == "cov.cauchy"){
    X <- mvtnorm::rmvt(n = n, sigma = Sigma, df = 1, delta = mu)
  }else(return("no such type"))

  return(X)
}



#############################################################

#' @keywords internal
#' @noRd
# ARMA(2,2) with N(0, 1) Innovations
rarma_normal <- function(n, mu, arma, center=0, sig=1){
  arma_normal<- function(nsim, p, center, sig, arma){
    e = stats::arima.sim(n = p, model = arma,
                         rand.gen = function(n, ...) stats::rnorm(n, mean=center, sd=sig))
    return(e)
  }
  e = sapply(rep(1, n), arma_normal, p=length(mu), center=center, sig=sig, arma=arma)
  X = t(mu + e)
  return(X)
}


#' @keywords internal
#' @noRd
# ARMA(2,2) with scaled t(3) Innovations
rarma_t <- function(n, mu, arma, q = 3){
  arma_t<- function(nsim, p, q, arma){
    e = stats::arima.sim(n = p, model = arma,
                         rand.gen = function(n, ...) stats::rt(n, df=q)/sqrt(q/(q-2)))
    return(e)
  }
  e = sapply(rep(1, n), arma_t, p=length(mu), q=q, arma=arma)
  X = t(mu + e)
  return(X)
}


#' @keywords internal
#' @noRd
# ARMA(2,2) with Cauchy(0, 0.1) Innovations
rarma_cauchy <- function(n, mu, arma, center=0, sig=0.1){
  arma_cauchy <- function(nsim, p, center, sig, arma){
    e = stats::arima.sim(n = p, model = arma,
                         rand.gen = function(n, ...) stats::rcauchy(n, location=center, scale=sig))
    return(e)
  }

  e = sapply(rep(1,n), arma_cauchy, p=length(mu), center=center, sig=sig, arma=arma)
  X = t(mu + e)
  return(X)
}


#' @keywords internal
#' @noRd
# ARMA(2,2) with centered Gamma(4, 0.5) Innovations
rarma_gamma <- function(n, mu, arma, a=4, s=0.5){
  arma_gamma <- function(nsim, p, a, s, arma){
    e = stats::arima.sim(n = p, model = arma,
                         rand.gen = function(n, ...) (stats::rgamma(n, shape=a, scale = s) - a*s)/sqrt(a*s^2))
    return(e)
  }

  e = sapply(rep(1,n), arma_gamma, p=length(mu), a=a, s=s, arma=arma)
  X = t(mu + e)
  return(X)
}


#' @keywords internal
#' @noRd
# Long-range dependent N(0, 1)
rld_normal <- function(n, mu, H = .7, center=0, sig=1){
  p = length(mu)
  R <- matrix(0, p, p)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      k <- abs(i - j)
      R[i, j] = R[j, i] = .5 * ((k + 1) ^ (2 * H) + (k - 1) ^ (2 * H) - 2 * k ^ (2 * H))
    }
  }
  diag(R) <- 1
  U <- chol(R)
  rand.gen = stats::rnorm(n * p, mean = center, sd = sig)
  e <- matrix(rand.gen, nrow = n) %*% U
  X = t(mu + t(e))
  return(X)
}



#' @keywords internal
#' @noRd
# Long-range dependent scaled t(3)
rld_t <- function(n, mu, H = .7, q = 3){
  p = length(mu)
  R <- matrix(0, p, p)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      k <- abs(i - j)
      R[i, j] = R[j, i] = .5 * ((k + 1) ^ (2 * H) + (k - 1) ^ (2 * H) - 2 * k ^ (2 * H))
    }
  }
  diag(R) <- 1
  U <- chol(R)
  rand.gen = stats::rt(n * p, df = q)/sqrt(q/(q-2))
  e <- matrix(rand.gen, nrow = n) %*% U
  X = t(mu + t(e))
  return(X)
}



#' @keywords internal
#' @noRd
# Long-range dependent Cauchy(0, 0.1)
rld_cauchy <- function(n, mu, H = .7, center=0, sig=0.1){
  p = length(mu)
  R <- matrix(0, p, p)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      k <- abs(i - j)
      R[i, j] = R[j, i] = .5 * ((k + 1) ^ (2 * H) + (k - 1) ^ (2 * H) - 2 * k ^ (2 * H))
    }
  }
  diag(R) <- 1
  U <- chol(R)
  rand.gen = stats::rcauchy(n * p, location=center, scale=sig)
  e <- matrix(rand.gen, nrow = n) %*% U
  X = t(mu + t(e))
  return(X)
}


#' @keywords internal
#' @noRd
# Long-range dependent centered Gamma(4, 0.5)
rld_gamma <- function(n, mu, H = .7, a=4, s=0.5){
  p = length(mu)
  R <- matrix(0, p, p)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      k <- abs(i - j)
      R[i, j] = R[j, i] = .5 * ((k + 1) ^ (2 * H) + (k - 1) ^ (2 * H) - 2 * k ^ (2 * H))
    }
  }
  diag(R) <- 1
  U <- chol(R)
  rand.gen = (stats::rgamma(n * p, shape=a, scale = s) - a*s)/sqrt(a*s^2)
  e <- matrix(rand.gen, nrow = n) %*% U
  X = t(mu + t(e))
  return(X)
}
