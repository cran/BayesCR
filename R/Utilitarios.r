###### Distribuições Acumuladas das Normal Independ (NI)
################################################################################
   
PearsonVII <- function(y,mu,sigma2,nu,delta)
     {
     Acum <- z <- vector(mode = "numeric", length = length(y))
  sigma2a <- sigma2*(delta/nu)
     for (i in 1:length(y))
     { 
      z[i] <- (y[i]-mu)/sqrt(sigma2a)
     Acum[i] <- pt(z[i],df=nu)
     }
    return(Acum)
     }


AcumSlash <- function(y,mu,sigma2,nu)
  {
  Acum <- z <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y))
   { 
  z[i] <- (y[i]-mu)/sqrt(sigma2)
    f1 <- function(u) nu*u^(nu-1)*pnorm(z[i]*sqrt(u))
Acum[i]<- integrate(f1,0,1)$value	 	
   }
return(Acum)
     }

AcumNormalC <- function(y,mu,sigma2,nu)
   {
  Acum <- vector(mode = "numeric", length = length(y))
 for (i in 1:length(y))
   { 
  eta  <- nu[1]    
  gama <- nu[2]
Acum[i] <- eta*pnorm(y[i],mu,sqrt(sigma2/gama)) + (1-eta)*pnorm(y[i],mu,sqrt(sigma2))
   }
return(Acum)
   }


###### Densidades das Normais Independentes
################################################################################

dPearsonVII<- function(y,mu,sigma2,nu,delta)
    {
     f <-  z <- vector(mode = "numeric", length = length(y))
  sigma2a <- sigma2*(delta/nu)
     for (i in 1:length(y))
     { 
     z[i] <- (y[i]-mu)/sqrt(sigma2a)
     f[i] <- dt(z[i],df=nu)/sqrt(sigma2a)
     }
     return(f)
     }

dSlash <- function(y,mu,sigma2,nu)
     {
     resp <- z <- vector(mode = "numeric", length = length(y))
     for (i in 1:length(y)) 
     {
     z[i] <- (y[i]-mu)/sqrt(sigma2)
     f1 <- function(u) nu*u^(nu-0.5)*dnorm(z[i]*sqrt(u))/sqrt(sigma2)
     resp[i] <- integrate(f1,0,1)$value	 	
     }
     return(resp)
     }

dNormalC <- function(y,mu,sigma2,nu)
    {
     Acum <- vector(mode = "numeric", length = length(y))
     for (i in 1:length(y))
     { 
    eta  <- nu[1]    
    gama <- nu[2]
   Acum[i]  <- eta*dnorm(y[i],mu,sqrt(sigma2/gama)) + (1-eta)*dnorm(y[i],mu,sqrt(sigma2))
     }
   return(Acum)
    }

################################################################################
# Calculo da Log Verossimilhança para a Normal Contaminada
################################################################################

verosCN <- function(auxf,cc,cens,sigma2,nu)
{
  auxf1 <- sqrt(auxf)
  if (cens=="1")
  {
  ver1 <-(sum(log(dNormalC(auxf1[cc==0],0,1,nu)/sqrt(sigma2)))+ sum(log(AcumNormalC(auxf1[cc==1],0,1,nu))))
  }
  if (cens=="2")
  {
  ver1 <-(sum(log(dNormalC(auxf1[cc==0],0,1,nu)/sqrt(sigma2)))+ sum(log(AcumNormalC(auxf1[cc==1],0,1,nu))))
  }
  return(ver1)
}