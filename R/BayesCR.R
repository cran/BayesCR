#' @title BayesCR: Bayesian Analysis of Censored Regression Models Under Scale Mixture of Skew Normal Distributions
#'
#' @description Propose a parametric fit for censored linear regression models based on SMSN distributions, from a Bayesian perspective. Also, generates SMSN random variables.
#'
#' @references
#' Monique B. Massuia, Aldo M. Garay, Celso R. Cabral and Victor H. Lachos. "Bayesian analysis of censored linear regression models with scale mixtures of skew-normal distributions". Statistics and Its Interface, 2017, vol. 10, pages 425-439
#'
#' @author
#' Aldo M. Garay \email{medina_garay@yahoo.com}, Monique B. Massuia \email{moniquemassuia@gmail.com}, 
#' Victor H. Lachos \email{hlachos@ime.unicamp.br} and Eraldo B. Anjos Filho \email{ebdaf1@de.ufpe.br}
#'
#' Maintainer: Aldo M. Garay \email{medina_garay@yahoo.com}
#'
#' @seealso \code{\link{Bayes.CR}}, \code{\link{rSMSN}}, \code{\link{motorettes}}
#' 
#' @importFrom stats coefficients cov deriv dnorm dt integrate lm pnorm pt qnorm rbeta rbinom rlnorm rnorm rt runif sd var
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats pgamma qgamma rgamma
#' @import rootSolve 
#' @import truncdist 
#' @import mvtnorm 
#' @import mnormt 
#'
#' @docType package
#' @name BayesCR
#' 
"_PACKAGE"