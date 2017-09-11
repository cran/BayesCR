#' Accelerated Life Tests On Electrical Insulation
#'
#' Accelerated life tests on electrical insulation in motorettes with censoring times. 
#'
#' @format A data frame with 40 observed times of life tests on electrical insulation in motorettes at 
#' four different temperatures (150C, 170c, 190c and 200c). y corresponds to log10 of the failure time 
#' (or end of study time, in case of right censored observations), x corresponds to (100/(temperature + 273.2)) 
#' and cc is a indicator of censoring (1 if censored, 0 if not).
#' 
#' @source Tan, M., Tian, G. L. and Ng, K. W. (2009). Bayesian Missing Data Problems: EM, Data Augmentation and Noniterative Computation.
#'
"motorettes"