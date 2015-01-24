\name{rSMSN}
\alias{rSMSN}
\title{ Generate SMSN Random Variables }
\description{
Generate random variables with one of the following distributions: Normal, Student-t, Contaminated Normal, Slash, Skew-Normal, Skew-t and Skew-Slash.
}
\usage{ rSMSN(n,mu,sigma2,lambda,nu,dist) }
\arguments{
\item{n}{ Number of observations to be generated. }
\item{mu}{ Location parameter. }
\item{sigma2}{ Scale parameter. }
\item{lambda}{ Shape parameter (control skewness). Only must be provided for Skew-Normal, Skew-t and Skew-Slash distributions. }
\item{nu}{Degree of freedom. Must not be provided for Normal and Skew-Normal distribution. Must be a vector of length 2 for Contaminated-Normal distribution. }
\item{dist}{ Distribution to be used: "Normal" for Normal model, "T" for Student-t model, "Slash" for slash model, "NormalC" for contaminated Normal model, "SN" for Skew-Normal model, "ST" for Skew-t model and "SSL" for Skew-Slash model. }

}
\details{ If Y follows a Contaminated Normal model, than a observation y comes from a normal distribution with mean "mu"" and variance "sigma2/rho" with probabilty "nu" and comes from a normal distribution with mean "mu" and variace "sigma2" with probability "1-nu". }

\references{
 Aldo M. Garay, Victor H. Lachos, Heleno Bolfarine, Celso R. Cabral. 2013. "Bayesian analysis of censored regression models with scale mixture of normal distributions".
 Vicente G. Cancho, Dipak K. Dey, Victor H. Lachos and Marinho G. Andrade. 2010. "Bayesian nonlinear regression models with scale mixtures of skew-normal distributions: estimation and case influence diagnostics".
}
\author{ Aldo M. Garay \email{medina_garay@yahoo.com}, Monique Bettio Massuia \email{moniquemassuia@gmail.com}, Victor Hugo Lachos \email{hlachos@ime.unicamp.br} }
\examples{
# Generate a sample with 100 observations of a symmetric Student-t distribution
sample <- rSMSN(n=100,mu=5,sigma2=2,lambda=0,nu=3,dist="T")
}