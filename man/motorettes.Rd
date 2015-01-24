\name{motorettes}
\docType{data}
\alias{motorettes}
\title{ Accelerated Life Tests On Electrical Insulation }
\description{ 
Accelerated life tests on electrical insulation in motorettes with censoring times. 
}
\usage{data(motorettes)}
\format{ 
A data frame with 40 observed times of life tests on electrical insulation in motorettes at four different temperatures (150C, 170c, 190c and 200c). y corresponds to log10 of the failure time (or end of study time, in case of right censored observations), x corresponds to (100/(temperature + 273.2)) and cc is a indicator of censoring (1 if censored, 0 if not).  
}

\encoding{latin1}

\source{ 
Tan, M., Tian, G. L. and Ng, K. W. (2009). Bayesian Missing Data Problems: EM, Data Augmentation and Noniterative Computation. 
}

\examples{
##Load the data
data(motorettes)
attach(motorettes)

##Set design matrix
x <- cbind(1,x)

##Fits a right censored normal model
Normal <- Bayes.CR(cc,x,y,cens="right",dist="Normal",n.thin=10,burnin=200,n.iter=800,
n.chains=1,chain="TRUE")
}

\keyword{datasets}
