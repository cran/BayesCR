\name{wage.rate}
\docType{data}
\alias{wage.rate}
\title{ Wage rate of 753 women}
\description{ Wage rate of 753 married white women with left censoring 
}
\usage{data(wage.rate)}
\format{ 
A data frame with 753 observed wage rate of married white women in 1975. Some wage rate are set equal to zero, this means that those wives did not work in 1975, therefore, these observations are considered left censored at zero. 
}

\encoding{latin1}

\source{ 
Mroz, T.A. 1987. "The sensitivity of an empirical model of married women´s hours of work to economic and statistical assumptions". 
}

\examples{
\dontrun{
## Load the data
data(wage.rate)

## Set the response y and covariate x

   attach(wage.rate)
   cc <- (inlf==0)+0
    x <- cbind(1,age,educ,kidslt6,kidsge6)
    y <- wage

## Fits a left censored Normal model 

Normal<- Bayes.CR(cc, x,y,cens="1",type="Normal",influence="FALSE",prior=NULL,
hyper=NULL,n.thin=10,burnin=100,n.iter=2000,n.chains=2,chain="TRUE")


## Fits a left censored T model, using Hierarchical prior for degree freedon

T<- Bayes.CR(cc, x,y,cens="1",type="T",influence="FALSE",prior="Hierar",
hyper=NULL,n.thin=10,burnin=100,n.iter=2000, n.chains=2,chain="TRUE")


## Fits a left censored Slash model

Slash<- Bayes.CR(cc, x,y,cens="1",type="Slash",influence="FALSE",prior=NULL,
hyper=NULL,n.thin=15,burnin=100,n.iter=2000, n.chains=2,chain="TRUE")

}
}
\keyword{datasets}



