\name{Bayes.CR}
\alias{Bayes.CR}
\title{ Bayesian Analysis of Censored Linear Regression Models with Scale Mixtures of Normal (SMN)
 Distributions }
\description{
Propose a robust parametric  model of censored linear regression models based on SMN distributions from a Bayesian perspective.
}
\usage{Bayes.CR(cc, x,y,cens="1",type="Normal",influence="FALSE",prior=NULL,
hyper=NULL,n.thin=10,burnin=100,n.iter=2000, n.chains=2,chain="TRUE")}
\arguments{
\item{cc}{ Vector of censoring indicators. For each observation: 0 if non-censored, 1 if censored }
\item{x}{ Matrix or vector of covariates }
\item{y}{ Vector of responses in case of right/left censoring }
\item{cens}{ "1" for left censoring, "2" for right censoring }
\item{type}{ Distribution to be used: "Normal" for Normal model, "T" for Student-t model, "Slash" for slash model and "NormalC" for contaminated Normal model }
\item{influence}{ "TRUE" or "FALSE". Indicates if the divergence measures should be plotted }
\item{prior}{ Prior distribution to be used for the degrees of freedom under the T model: "Exp" for exponential distribution, "Jeffreys" for Jeffreys prior, "Unif" for Uniforme distribution and "Hierar" for Hierarchical prior. }
\item{hyper}{ Value of hyperparameter for the exponential prior. Must not be provided in case of others prior distributions }
\item{n.thin}{ Lag for posterior sample }
\item{burnin}{ Burn-in for posterior sample }
\item{n.iter}{ The number of iterations to be considered. }
\item{n.chains}{ The number of chains to be considered. It must be less than 5. }
\item{chain}{ If "TRUE", then all the chain is stored for posterior analysis. }
}
\details{ The following prior distributions were considered: Normal for the regression coefficients, inverse gamma for the scale parameter. }
\value{ 
\item{Mean}{ Posterior mean for the parameters }
\item{Sd}{ Standard deviations for the parameters }
\item{HPD}{ HPD(95\%) interval for the parameters }
\item{LPML}{ Log-marginal pseudo likelihood for model selection. }
\item{DIC}{ DIC criterion for model selection. }
\item{EAIC}{ EAIC criterion for model selection. }
\item{EBIC}{ EBIC criterion for model selection. }
}
\references{
 Aldo M. Garay, Victor H. Lachos, Heleno Bolfarine, Celso R. Cabral. "Bayesian Analysis of Censored Regression Models with Scale Mixture of Normal Distributions".
 }
\author{ Aldo M. Garay \email{aldowmg@ime.usp.br}, Victor Hugo Lachos \email{hlachos@ime.unicamp.br} }
\examples{ ## see examples in \code{\link{wage.rate}} }
\keyword{ censored linear regression }


