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


Rhat1 <- function(param,n.iter,burnin,n.chains,n.thin) 
{
 param <- as.matrix(param)
 efect <- (n.iter-burnin)/n.thin
     p <- ncol(param)
   mat <- matrix(param,nrow=efect,ncol=n.chains*p)    
  rhat <- matrix(0,nrow=p,ncol=1)                                          
   for(i in 1:p)
   {
    l1 <- 2*(i-1)+1
    c1 <- 2*i
    rhat[i,1] <- Rhat(mat[,l1:c1]) 
   }
return(rhat=rhat)
}


Rhat <- function(mat) 
{
  m <- ncol(mat)
  n <- nrow(mat)
  b <- apply(mat,2,mean)
  B <- sum((b-mean(mat))^2)*n/(m-1)
  w <- apply(mat,2,var)
  W <- mean(w)
  s2hat <- (n-1)/n*W + B/n
  Vhat <- s2hat + B/m/n 
  covWB <- n /m * (cov(w,b^2)-2*mean(b)*cov(w,b))
  varV <- (n-1)^2 / n^2 * var(w)/m +
          (m+1)^2 / m^2 / n^2 * 2*B^2/(m-1) +
          2 * (m-1)*(n-1)/m/n^2 * covWB
  df <- 2 * Vhat^2 / varV
  R <- sqrt((df+3) * Vhat / (df+1) / W)
  return(R)
 }
 
 
hpd <- function(x, alpha)
{
   n <- length(x)
   m <- max(1, ceiling(alpha * n))

   y <- sort(x)
   a <- y[1:m]
   b <- y[(n - m + 1):n]

   i <- order(b - a)[1]
   structure(c(a[i], b[i]), names = c("Lower Bound", "Upper Bound"))
}
