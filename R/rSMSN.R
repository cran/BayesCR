#' @export
#'
#' @title Generate SMSN Random Variables
#'
#' @description \code{rSMSN} Generate random variables with one of the following distributions: 
#' Normal, Student-t, Contaminated Normal, Slash, Skew-Normal, Skew-t and Skew-Slash.
#'
#' @param n Number of observations to be generated. 
#' @param mu Location parameter. 
#' @param sigma2 Scale parameter. 
#' @param lambda Shape parameter (control skewness). Only must be provided for Skew-Normal, 
#' Skew-t and Skew-Slash distributions. 
#' @param nu Degree of freedom. Must not be provided for Normal and Skew-Normal distribution. 
#' Must be a vector of length 2 for Contaminated-Normal distribution. 
#' @param dist Distribution to be used: "Normal" for Normal model, "T" for Student-t model, 
#' "Slash" for slash model, "NormalC" for contaminated Normal model, "SN" for Skew-Normal model, 
#' "ST" for Skew-t model and "SSL" for Skew-Slash model. 
#' 
#' @details If Y follows a Contaminated Normal model, than a observation y comes from a normal 
#' distribution with mean "mu"" and variance "sigma2/rho" with probabilty "nu" 
#' and comes from a normal distribution with mean "mu" and variace "sigma2" 
#' with probability "1-nu". 
#'
#' @seealso \code{\link{Bayes.CR}}, \code{\link{motorettes}}
#'
#' @examples # Generate a sample with 100 observations of a symmetric Student-t distribution
#' 
#' sample <- rSMSN(n=100,mu=5,sigma2=2,lambda=0,nu=3,dist="T")

rSMSN <- function(n, mu, sigma2, lambda, nu, dist) {
  if (dist == "T" | dist == "Slash" | dist == "NormalC" | dist == "ST" | 
      dist == "SSL") {
    if (length(nu) == 0) 
      stop("nu must be provided.")
  }
  if (dist == "NormalC") {
    if (length(nu) != 2) 
      stop("nu must be a vector of size 2.")
  }
  if (dist == "SN" | dist == "SSL" | dist == "ST") {
    if (length(lambda) == 0) 
      stop("lambda must be provided.")
  }
  
  q <- length(beta)
  if (dist == "Normal") {
    y <- mu + rnorm(n, 0, sd = sqrt(sigma2))
  }
  if (dist == "T") {
    y <- mu + sqrt(sigma2) * rt(n, df = nu)
  }
  if (dist == "NormalC") {
    p <- runif(n)
    u <- rep(1, n)
    u[p < nu[1]] <- nu[2]
    y <- mu + rnorm(n, mean = 0, sd = sqrt(sigma2))/sqrt(u)
  }
  if (dist == "Slash") {
    u <- rbeta(n = n, shape1 = nu, shape2 = 1)
    y <- mu + rnorm(n, mean = 0, sd = sqrt(sigma2))/sqrt(u)
  }
  if (dist == "SN" | dist == "ST" | dist == "SSL") {
    y <- rep(0, n)
    if (dist == "SN") {
      u <- rep(1, n)
    }
    if (dist == "ST") {
      u <- rgamma(n = n, shape = nu/2, rate = nu/2)
    }
    if (dist == "SSL") {
      u <- rbeta(n = n, shape1 = nu, shape2 = 1)
    }
    deltinha <- lambda/sqrt(1 + lambda^2)
    Delta <- sqrt(sigma2) * deltinha
    tau <- sigma2 * (1 - deltinha^2)
    
    T1 <- rnorm(n)
    T0 <- rnorm(n)
    T <- abs(T0) * u^(-1/2)
    y <- mu + Delta * T + u^(-1/2) * sqrt(tau) * T1
  }
  return(y)
}




