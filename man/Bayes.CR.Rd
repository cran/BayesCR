\name{Bayes.CR}
\alias{Bayes.CR}
\title{ Bayesian Analysis of Censored Regression Models Under Scale Mixture of Skew Normal Distributions }
\description{
Propose a parametric fit for censored linear regression models based on SMSN distributions, from a Bayesian perspective.
}
\usage{Bayes.CR(cc,x,y,cens="left",dist="Normal",criteria="FALSE",influence="FALSE",
spacing=NULL,prior=NULL,hyper=NULL,n.thin=10,burnin=100,n.iter=2000,
n.chains=2,chain="TRUE")}
\arguments{
\item{cc}{ Vector of censoring indicators. For each observation: 0 if non-censored, 1 if censored. }
\item{x}{ Matrix or vector of covariates. }
\item{y}{ Vector of responses in case of right/left censoring. }
\item{cens}{ "left" for left censoring, "right" for right censoring. }
\item{dist}{ Distribution to be used: "Normal" for Normal model, "T" for Student-t model, "Slash" for slash model, "NormalC" for contaminated Normal model, "SN" for Skew-Normal model, "ST" for Skew-t model and "SSL" for Skew-Slash model. }
\item{criteria}{ "TRUE" or "FALSE". Indicates if model selection criteria (LPML, DIC, EAIC, EBIC and WAIC) should be computed. }
\item{influence}{ "TRUE" or "FALSE". Indicates if the divergence measures (KL divergence, J, L and Chi Distance) should be computed. }
\item{spacing}{ Should only be specified if at least one of "influence" or "criteria" is TRUE. This is the lag between observations of the final chain (after burn-in and thinning) used to compute these measures. If spacing=1, all the chain is used. }
\item{prior}{ Prior distribution to be used for the degrees of freedom under Student-t model: "Exp" for exponential distribution, "Jeffreys" for Jeffreys prior, "Unif" for Uniforme distribution and "Hierar" for Hierarchical prior (exponential with a parameter that follows a uniform distribution). Must be "NULL" for other models. }
\item{hyper}{ Value of hyperparameter for the exponential prior. Must not be provided in case of others prior distributions. }
\item{n.thin}{ Lag for posterior sample. }
\item{burnin}{ Burn-in for posterior sample. }
\item{n.iter}{ The number of iterations to be considered (before burnin and thinning). }
\item{n.chains}{ The number of chains to be considered. It must be less than 5. }
\item{chain}{ If "TRUE", all the posterior chains are stored for posterior analysis. }
}
\details{Specification of the priors distributions is given in reference papers (Garay et. al 2013 and Cancho et. al 2010). See Gelman et. al for the difference between the two versions of WAIC criterion. Calculations under the Skew-slash model may take a while, as it involves numerical integrations - you may want to specify big values to "spacing" under this model. For the Contaminated Normal model, a observation y comes from a normal distribution with mean "x beta" and variance "sigma2/rho" with probabilty "nu" and comes from a normal distribution with mean "x beta" and variance "sigma2" with probability 1-"nu". }
\value{ 
\item{Mean}{ Posterior mean for the parameters. }
\item{Sd}{ Standard deviations for the parameters. }
\item{HPD}{ HPD(95\%) interval for the parameters. }
\item{LPML}{ Log-marginal pseudo likelihood for model selection. }
\item{DIC}{ DIC criterion for model selection. }
\item{EAIC}{ EAIC criterion for model selection. }
\item{EBIC}{ EBIC criterion for model selection. }
\item{WAIC1}{First version of Watanabe-Akaike information criterion. }
\item{WAIC2}{Second version of Watanabe-Akaike information criterion. }
}
\references{
 Aldo M. Garay, Victor H. Lachos, Heleno Bolfarine, Celso R. Cabral. 2013. "Bayesian analysis of censored regression models with scale mixture of normal distributions".
 Vicente G. Cancho, Dipak K. Dey, Victor H. Lachos and Marinho G. Andrade. 2010. "Bayesian nonlinear regression models with scale mixtures of skew-normal distributions: estimation and case influence diagnostics".
 Andrew Gelman, Jessica Hwang, Aki Vehtari. 2013. "Understanding predictive information criteria for Bayesian models". 
}
\author{ Aldo M. Garay \email{medina_garay@yahoo.com}, Monique Bettio Massuia \email{moniquemassuia@gmail.com}, Victor Hugo Lachos \email{hlachos@ime.unicamp.br} }
\examples{ ##see examples in \code{\link{motorettes}} }
\keyword{ censored linear regression }


