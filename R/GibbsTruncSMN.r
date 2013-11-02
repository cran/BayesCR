#################################################################################
########## ALGORITMO GIBBS SAMPLER PARA MODELOS CENSURADOS SMN     ###############
# algoritmo Gibbs para o caso normal, T, Slash, Normal Contaminado para os dois 
# tipos de censuras. 
# Para o caso t-Student usamos diferentes dist. a prioris para \nu (desconhecido)       
################################################################################
 
### Apriori Jeffreys
##########################

MHnu<-function(last,U,lambda,prior="Jeffreys",hyper)
{
     n <- length(U) 
     if(prior=="Jeffreys")         # Apriori Jeffreys
     {
            gJeffreys <- function(nu,U)
            {
            ###USADO NO MH
               n<-length(U)
              ff<- log(sqrt(nu/(nu+3))*sqrt(trigamma(nu/2)-trigamma((nu+1)/2)-2*(nu+3)/((nu)*(nu+1)^2)))+0.5*n*nu*log(nu/2)+(0.5*nu)*sum(log(U)-U)-n*log(gamma(nu/2))
              return(ff)
            }
 Fonseca1 <-deriv(~log(sqrt(nu/(nu+3))*sqrt(trigamma(nu/2)-trigamma((nu+1)/2)-2*(nu+3)/((nu)*(nu+1)^2))) +0.5*n*nu*log(nu/2)-n*log(gamma(nu/2)),c("nu"),function(nu){},hessian=TRUE)
 Fonseca2 <- deriv(~(0.5*nu)*(log(U)-U),c("nu"),function(U,nu){},hessian=TRUE)
 aux1 <- Fonseca1(last)
 aux2 <- Fonseca2(U,last)                                                   
   q1 <- attr(aux1,"gradient")[1]+sum(attr(aux2,"gradient"))
   q2 <- attr(aux1,"hessian")[1]+sum(attr(aux2,"hessian"))
   aw <- last-q1/q2
   bw <- max(0.001,-1/q2)
 cand <- rtrunc(1, spec="norm",a=2.1, b=100, mean = aw, sd = sqrt(bw)) 
 alfa <- (exp(gJeffreys(cand,U))/exp(gJeffreys(last,U)))*(dtrunc(last,spec="norm",a=2.1, b=100, mean = aw, sd = sqrt(bw))/dtrunc(cand,spec="norm", a=2.1, b=100, mean = aw, sd = sqrt(bw)))
      }

      if(prior=="Exp")           # Exponencial (hiperpar)
      {
            gExp <- function(nu,U,hyper)
            {
            ###USADO NO MH
               n<-length(U)
              ff<- (-hyper*nu) + 0.5*n*nu*log(nu/2)+(0.5*nu)*sum(log(U)-U)-n*log(gamma(nu/2))
              return(ff)
            }
            Fonseca1<- deriv(~ (-hyper*nu) + 0.5*n*nu*log(nu/2)-n*log(gamma(nu/2)),c("nu"),function(nu){},hessian=TRUE)
            Fonseca2<- deriv(~(0.5*nu)*(log(U)-U),c("nu"),function(U,nu){},hessian=TRUE)
                aux1<-Fonseca1(last)
                 aux2<-Fonseca2(U,last)                                                   
                    q1<-attr(aux1,"gradient")[1]+sum(attr(aux2,"gradient"))
                    q2<-attr(aux1,"hessian")[1]+sum(attr(aux2,"hessian"))
                    aw<- last-q1/q2
                    bw<- max(0.001,-1/q2)
                  cand <- rtrunc(1, spec="norm",a=2.1, b=100, mean = aw, sd = sqrt(bw))                  
                  alfa<-(exp(gExp(cand,U,hyper))/exp(gExp(last,U,hyper)))*(dtrunc(last,spec="norm",a=2.1, b=100, mean = aw, sd = sqrt(bw))/dtrunc(cand,spec="norm", a=2.1, b=100, mean = aw, sd = sqrt(bw)))
      }

      if(prior=="Unif")           # Uniforme(2,100)
      {
            gUnif <- function(nu,U)
            {
            ###USADO NO MH
               n<-length(U)
              ff<- 0.5*n*nu*log(nu/2)+(0.5*nu)*sum(log(U)-U)-n*log(gamma(nu/2))
              return(ff)
            }
            Fonseca1<- deriv(~ 0.5*n*nu*log(nu/2)-n*log(gamma(nu/2)),c("nu"),function(nu){},hessian=TRUE)
            Fonseca2<- deriv(~(0.5*nu)*(log(U)-U),c("nu"),function(U,nu){},hessian=TRUE)
                aux1<-Fonseca1(last)
                aux2<-Fonseca2(U,last)                                                   
                  q1<-attr(aux1,"gradient")[1]+sum(attr(aux2,"gradient"))
                  q2<-attr(aux1,"hessian")[1]+sum(attr(aux2,"hessian"))
                  aw<- last-q1/q2
                  bw<- max(0.001,-1/q2)
                cand <- rtrunc(1, spec="norm",a=2.1, b=100, mean = aw, sd = sqrt(bw))                                  
                alfa<-(exp(gUnif(cand,U))/exp(gUnif(last,U)))*(dtrunc(last,spec="norm",a=2.1, b=100, mean = aw, sd = sqrt(bw))/dtrunc(cand,spec="norm", a=2.1, b=100, mean = aw, sd = sqrt(bw)))

      }
        
      if(prior=="Hierar")           # Hierárquica
      {
            gHierar <- function(nu,U)
            {
            ###USADO NO MH
               n<-length(U)
              ff<- (-lambda*nu) + 0.5*n*nu*log(nu/2)+(0.5*nu)*sum(log(U)-U)-n*log(gamma(nu/2))
              return(ff)
            }
            Fonseca1<- deriv(~ (-lambda*nu) + 0.5*n*nu*log(nu/2)-n*log(gamma(nu/2)),c("nu"),function(nu){},hessian=TRUE)
            Fonseca2<- deriv(~(0.5*nu)*(log(U)-U),c("nu"),function(U,nu){},hessian=TRUE)
                aux1<- Fonseca1(last)
                 aux2<-Fonseca2(U,last)                                                   
                   q1<-attr(aux1,"gradient")[1]+sum(attr(aux2,"gradient"))
                   q2<-attr(aux1,"hessian")[1]+sum(attr(aux2,"hessian"))
                   aw<- last-q1/q2
                   bw<- max(0.001,-1/q2)
                 cand <- rtrunc(1, spec="norm",a=2.1, b=100, mean = aw, sd = sqrt(bw))                                   
                 alfa<-(exp( gHierar(cand,U))/exp( gHierar(last,U)))*(dtrunc(last,spec="norm",a=2.1, b=100, mean = aw, sd = sqrt(bw))/dtrunc(cand,spec="norm", a=2.1, b=100, mean = aw, sd = sqrt(bw)))
      }
      
      ifelse(runif(1) < min(alfa, 1), last <- cand, last<-last)
      return(last)
}

MHrhoCN <-function(last,U,cc,sigma2,nu,s0,s1,cens)
{
last1 <- last/(1-last)
 desv <- .001
 cand <- rlnorm(1,meanlog=log(last1), sdlog=sqrt(desv))

  ver <- verosCN(U,cc,cens=cens,sigma2,nu)
    g <- function(r,ver)
    {
   r1 <- r/(1+r)
    ff<- (r1)^(s0-1)*(1-r1)^(s1-1)*ver*(1/(1+r)^2)
    return(ff)
    }
alfa <- g(cand,ver)*cand/(g(last1,ver)*last1)
ifelse(runif(1) < min(alfa, 1), last <- cand, last<-last)
return(last)
}

################################################################################
# ALGORITMO DE GIBBS
################################################################################

GibbsTruncSMN<-function(cc,y,x,n.iter,n.thin,burnin,type="normal", cens="1",prior="Jeffreys",hyper=hyper, n.chains=n.chains)
{
        cad <- 0
          x <- as.matrix(x)
          p <- ncol(x)                                            
          n <- nrow(x)                                                
      sigma2<- lambda<- Aux <- V <- rho <- matrix(0,n.iter*n.chains,1)
       media<- matrix(0,n.iter*n.chains,n)
        beta<- matrix(0,n.iter*n.chains,p)
          y1<-y
   
          ###### Hiper parâmetros: comum ######
          mu0<-matrix(0,p,1)
       Sigma0<-1000*diag(p)
        alfa1<-0.1
        alfa2<-0.01
           a <-  r1 <- 0.02       ### Hiper Parametros de Lambda ~ Unif(0.02,0.5)
           b <-  s1 <- 0.5        ### Hiper Parametros de Lambda do Truncamento

      cat('% of iterations \n')
    if(type=="T")
    {
        U  <-matrix(0,n.iter*n.chains,n)
        nu <-matrix(0,n.iter*n.chains,1)
       set1<- set <- c()
   while(cad< n.chains)
   {
         cad <- cad +1
  #### Atualizador 
   Cont <- ceiling(n.iter/10)
 n.iter <- Cont*10
    Lin <- n.iter*(cad-1)
      W <- seq(Cont+Lin,n.iter*cad,Cont)
      W 
          reg <- lm(y ~ x[,2:p]) # Consideramos depois da seugna coluna porque "lm"  inlui a colunas de 1
            z <- 1 
    lambda[1+ Lin] <- 1
      beta[1+ Lin,]<- as.vector(coefficients(reg),mode="numeric")
    sigma2[1+ Lin] <- sum((y-x%*%(beta[1+Lin,]))^2)/(n-p)
     media[1+ Lin,]<-x%*%(beta[1+Lin,])

      U[1+Lin,]<-rinvgamma(1,1)       
      nu[1+Lin]<- 6
      cat("\n")
       Fil <- 2 + Lin
       Col <- n.iter*cad
      for (j in Fil:Col)
      {
    if(cens=="1")
    {
    y[cc==1] <- rtruncnorm(1, a=-Inf, b=y1[cc==1], mean = media[j-1,cc==1], sd = sqrt((U[j-1,cc==1])^(-1)*sigma2[j-1]))## Densid de Y/Resto - Cens Right
    }
   if(cens=="2")
    {          
    y[cc==1] <- rtruncnorm(1, a=y1[cc==1], b=Inf, mean = media[j-1,cc==1], sd = sqrt((U[j-1,cc==1])^(-1)*sigma2[j-1]))  ### Densid de Y/Resto - Cens Left
    }
    sigma2[j] <-rinvgamma(1, alfa1+n/2, scale = (sum(U[j-1,]*(y-media[j-1,])^2)/2+alfa2))                       ### Densid de sigma2/Resto   
        	xast <-sqrt(U[j-1,])*x
        	yast <-sqrt(U[j-1,])*y
       	SigmaA <-solve(solve(Sigma0)+t(xast)%*%(xast)/sigma2[j])
           muA <- SigmaA%*%(solve(Sigma0)%*%mu0+t(xast)%*%yast/sigma2[j])
       beta[j,]<-rmvnorm(1,muA,SigmaA)                                                                           ### Densid de beta/Resto   
     	media[j,]<-x%*%(beta[j,])
     	   meany2<-(y-media[j,])^2/sigma2[j-1]
        par2gam<- (meany2+nu[j-1])/2
       	  U[j,]<-rgamma(n,nu[j-1]/2+0.5,par2gam)                                                                  ### Densid de U/Resto    

        Aux[j] <- runif(1,0,1)
          V[j] <- pgamma(a,nu[j-1]) + (pgamma(b,2,nu[j-1])- pgamma(a,2,nu[j-1]))*Aux[j]
     lambda[j] <- qgamma(V[j],2,nu[j-1])                                                                                                                
         nu[j] <- MHnu(nu[j-1],U[j,],lambda[j],prior,hyper)                                                     ### Densid de Nu/Resto    

    ###### Cont Atualizador
    if(j==W[z])
    {
     z <- z +1
     BayCens(Fil,Col,j,cad)
    }
    ######
    }
    cat('\r')
   initial <- burnin 
       set1 <- seq(initial+n.thin +n.iter*(cad-1) ,n.iter*cad,n.thin)
        set <- c(set,set1)
   }
    return(list(beta=beta[set,],sigma2=sigma2[set],U=U[set,],nu=nu[set],mu=media[set,]))
    }
if(type=="Normal")
{

   set1<- set <- c()
   while(cad< n.chains)
   {
   
    cad <- cad +1
  #### Atualizador 
  
   Cont <- ceiling(n.iter/10)
 n.iter <- Cont*10
    Lin <- n.iter*(cad-1)
      W <- seq(Cont+Lin,n.iter*cad,Cont)
      W 
    reg <- lm(y ~ x[,2:p]) # Consideramos depois da seugna coluna porque "lm"  inlui a colunas de 1
      z <- 1 
#    lambda[1+ Lin] <- 1
      beta[1+ Lin,]<- as.vector(coefficients(reg),mode="numeric")
    sigma2[1+ Lin] <- sum((y-x%*%(beta[1+Lin,]))^2)/(n-p)
     media[1+ Lin,]<-x%*%(beta[1+Lin,])
    cat("\n")
    Fil <- 2 + Lin
    Col <- n.iter*cad
   for (j in Fil:Col)
   {
     if(cens=="1")  
     {
     y[cc==1] <- rtruncnorm(1, a=-Inf, b=y1[cc==1], mean = media[j-1,cc==1], sd = sqrt(sigma2[j-1]))                 ### Densid de Y/Resto - Cens Left
     }
     if(cens=="2")
     {
     y[cc==1] <- rtruncnorm(1, a=y1[cc==1], b=Inf, mean = media[j-1,cc==1], sd = sqrt(sigma2[j-1]))                  ### Densid de Y/Resto - Cens Right
     }
 sigma2[j]<- rinvgamma(1, alfa1+n/2, scale = (sum((y-media[j-1,])^2)/2+alfa2))                                        ### Densid de sigma2/Resto    
     xast <- x
   	 yast <- y
   SigmaA <- solve(solve(Sigma0)+t(xast)%*%(xast)/sigma2[j])
     	muA <- SigmaA%*%(solve(Sigma0)%*%mu0+t(xast)%*%yast/sigma2[j])
 beta[j,] <- rmvnorm(1,muA,SigmaA)                                                                                    ### Densidade de beta/Resto
media[j,] <- x%*%(beta[j,])

    ###### Cont Atualizador
    if(j==W[z])
    {
     z <- z +1
     BayCens(Fil,Col,j,cad)
    }
    ######
    }
    cat('\r')
   initial <- burnin 
       set1 <- seq(initial+n.thin +n.iter*(cad-1) ,n.iter*cad,n.thin)
        set <- c(set,set1)    
}
return(list(beta=beta[set,],sigma2=sigma2[set],mu=media[set,]))
}
if(type=="Slash")                     
{
# r1 <- s1 <- 100000
  U <-matrix(0,n.iter*n.chains,n)
 nu <-matrix(0,n.iter*n.chains,1)
set1<- set <- c()
   while(cad< n.chains)
   {
    r1 <- 0.01                             
    s1 <- 1
   cad <- cad +1
  #### Atualizador 
   Cont <- ceiling(n.iter/10)
 n.iter <- Cont*10
    Lin <- n.iter*(cad-1)
      W <- seq(Cont+Lin,n.iter*cad,Cont)
      W 
    reg <- lm(y ~ x[,2:p]) # Consideramos depois da seugna coluna porque "lm"  inlui a colunas de 1
      z <- 1 
      beta[1+ Lin,]<- as.vector(coefficients(reg),mode="numeric")
    sigma2[1+ Lin] <- sum((y-x%*%(beta[1+Lin,]))^2)/(n-p)
     media[1+ Lin,]<- x%*%(beta[1+Lin,])
#        U[1+Lin,] <- runif(n,0,1)       
         U[1+Lin,] <- rbeta(n,1,1)       
         nu[1+Lin] <- 1.5
      cat("\n")
       Fil <- 2 + Lin
       Col <- n.iter*cad
 for (j in Fil:Col)
 {
 if(cens=="1")
 {
  y[cc==1] <- rtruncnorm(1, a=-Inf, b=y1[cc==1], mean = media[j-1,cc==1], sd = sqrt((U[j-1,cc==1])^(-1)*sigma2[j-1])) ### Densid de Y/Resto - Cens Lefth
 }
 if(cens=="2")                                
 {
  y[cc==1] <- rtruncnorm(1, a=y1[cc==1], b=Inf, mean = media[j-1,cc==1], sd = sqrt((U[j-1,cc==1])^(-1)*sigma2[j-1]))  ### Densid de Y/Resto - Cens Right
 }
sigma2[j] <- rinvgamma(1, alfa1+n/2, scale = (sum(U[j-1,]*(y-media[j-1,])^2)/2+alfa2))                               ### Densid de sigma2/Resto   
    xast <- sqrt(U[j-1,])*x                     
    yast <- sqrt(U[j-1,])*y                                                                    
 	SigmaA <- solve(solve(Sigma0)+t(xast)%*%(xast)/sigma2[j])
     muA <- SigmaA%*%(solve(Sigma0)%*%mu0+t(xast)%*%yast/sigma2[j])
 beta[j,] <- rmvnorm(1,muA,SigmaA)                                                                                   ### Densid de beta/Resto   
media[j,] <- x%*%(beta[j,])      
   meany2 <- (y-media[j,])^2/sigma2[j-1]
  par2gam <- (meany2)/2                                                                                            
  	U[j,] <- rtrunc(n, "gamma", a =0.001 , b =1 , shape=nu[j-1]+0.5, scale=1/par2gam)                                    ### Densid de U/Resto    
    
#   print(U[j,])         
    rho[j] <- rtrunc(1, "gamma", a = r1, b =s1 , shape=2, scale=1/nu[j-1])
     nu[j] <- rgamma(1,shape= n+1,rate = rho[j]-sum(log(U[j,])))                                                     ### Densidade de Nu/Resto    

    ###### Cont Atualizador
    if(j==W[z])
    {
     z <- z +1
     BayCens(Fil,Col,j,cad)
    }
    ######
}
   cat('\r')
   initial <- burnin 
       set1 <- seq(initial+n.thin +n.iter*(cad-1) ,n.iter*cad,n.thin)
        set <- c(set,set1)   
}
return(list(beta=beta[set,],sigma2=sigma2[set],U=U[set,],nu=nu[set],mu=media[set,]))
}
if(type=="NormalC")
{
  U <-matrix(0,n.iter*n.chains,n)
 nu <-rho <- rho1 <- matrix(0,n.iter*n.chains,1)
 p1 <- p2 <- ptot <- A <- matrix(0,n.iter*n.chains,n)
set1<- set <- c()

   while(cad< n.chains)
   {
 AuxCN <- matrix(0,n,1)
 v0 <- s0 <- 2 
 v1 <- s1 <- 2
cad <- cad +1
  #### Atualizador 
   Cont <- ceiling(n.iter/10)
 n.iter <- Cont*10
    Lin <- n.iter*(cad-1)
      W <- seq(Cont+Lin,n.iter*cad,Cont)
      W 
    reg <- lm(y ~ x[,2:p]) # Consideramos depois da seugna coluna porque "lm"  inlui a colunas de 1
      z <- 1 
      beta[1+ Lin,]<- as.vector(coefficients(reg),mode="numeric")
    sigma2[1+ Lin] <- sum((y-x%*%(beta[1+Lin,]))^2)/(n-p)
     media[1+ Lin,]<- x%*%(beta[1+Lin,])
         U[1+Lin,] <- 1
       rho[1+Lin] <- nu[1+Lin] <- 0.4
      cat("\n")
       Fil <- 2 + Lin
       Col <- n.iter*cad
     cont  <- c()       
 for (j in Fil:Col)
 {
 if(cens=="1")
 {
  y[cc==1] <- rtruncnorm(1, a=-Inf, b=y1[cc==1], mean = media[j-1,cc==1], sd = sqrt((U[j-1,cc==1])^(-1)*sigma2[j-1])) ### Densid de Y/Resto - Cens Right
 }
 if(cens=="2")
 {
 y[cc==1] <- rtruncnorm(1, a=y1[cc==1], b=Inf, mean = media[j-1,cc==1], sd = sqrt((U[j-1,cc==1])^(-1)*sigma2[j-1]))  ### Densid de Y/Resto - Cens Left
 }
sigma2[j] <- rinvgamma(1, alfa1+n/2, scale = (sum(U[j-1,]*(y-media[j-1,])^2)/2+alfa2))                               ### Densid de sigma2/Resto   
    xast <- sqrt(U[j-1,])*x
    yast <- sqrt(U[j-1,])*y                                                                    
 	SigmaA <- solve(solve(Sigma0)+t(xast)%*%(xast)/sigma2[j])
     muA <- SigmaA%*%(solve(Sigma0)%*%mu0+t(xast)%*%yast/sigma2[j])
 beta[j,] <- rmvnorm(1,muA,SigmaA)                                                                                   ### Densid de beta/Resto   
media[j,] <- x%*%(beta[j,])
   meany2 <- (y-media[j,])^2/sigma2[j-1]
  par2gam <- (meany2)/2

#### Para este casso 
   p1[j,] <- nu[j-1]*sqrt(rho[j-1])*exp(-0.5*rho[j-1]*meany2)
   p2[j,] <- (1-nu[j-1])*exp(-0.5*meany2)
 ptot[j,] <- p1[j,]+p2[j,]
    A[j,] <- p1[j,]/ptot[j,]
   AuxCN <- rbern(n,A[j,])           
  U[j,] <- rho[j-1]*AuxCN + (1-AuxCN)                                                                                ### Densid de U/Resto    
             
   cont <- count(U[j,]==rho[j-1])
 nu[j] <- rbeta(1,v0+cont,v1+n-cont)                                                                                 ### Densid de Nu/Resto

#### Para a priori de de nu~Beta(v0,v1)

 rho1[j] <- MHrhoCN(rho[j-1],meany2,cc,sigma2[j],c(nu[j],rho[j-1]),s0,s1,cens)
 rho[j] <- rho1[j]/(1 + rho1[j])  
 
###### Cont Atualizador
if(j==W[z])
{
z <- z +1
BayCens(Fil,Col,j,cad)
}
######
}
   cat('\r')
   initial <- burnin 
       set1 <- seq(initial+n.thin +n.iter*(cad-1) ,n.iter*cad,n.thin)
        set <- c(set,set1)   
}
return(list(beta=beta[set,],sigma2=sigma2[set],U=U[set,],nu=nu[set],rho=rho[set],mu=media[set,]))
}
}

####################
## FIM do Gibbs
###################

## Updating Iterations
##########################

BayCens <- function(start.iter,end.iter,j,cad, ...)
{
    pb <- txtProgressBar(start.iter,end.iter,
                     initial = start.iter, style=3, width=10,
                     char=ifelse((cad ==1||cad==3),"+","*"))
   Sys.sleep(0.5); setTxtProgressBar(pb, j)
    cat("\r")
    cat("\r")
}


###########################

################################################################################
###                    CRITERIOS: CPO, DIC, EAIC, KL
################################################################################
      
LogVerosCens<-function(cc,y,mu,sigmae,nu,type="Normal", cens="2")
{
#Cens="1" : Left
#Cens= "2" : Rigth 

  m <- length(y)
  ver<- matrix(0,m,1)

  auxy<- matrix(0,m,1)

if(type=="Normal") 
{
     for (j in 1:m )
     {
         auxy[j]<- (y[j]-mu[j])/sqrt(sigmae)
         if (cens=="2")
         {
             if (cc[j]==0) {ver[j]<-dnorm(auxy[j])/sqrt(sigmae)}
             if (cc[j]==1) {ver[j]<-pnorm(-auxy[j])}
         }
         if (cens=="1")
         {
             if (cc[j]==0) {ver[j]<-dnorm(auxy[j])/sqrt(sigmae)}
             if (cc[j]==1) {ver[j]<-pnorm(auxy[j])}
         }
     }
 }    
 if(type=="T") 
 {                                                    
       for (j in 1:m )
       {
           auxy[j]<- (y[j]-mu[j])/sqrt(sigmae)
           if (cens=="2")
           {
               if (cc[j]==0) {ver[j]<-dt(auxy[j],df=nu)/sqrt(sigmae)}
               if (cc[j]==1) {ver[j]<-pt(-auxy[j],df=nu)}
           }
          if (cens=="1")
           {
               if (cc[j]==0) {ver[j]<-dt(auxy[j],df=nu)/sqrt(sigmae)}
               if (cc[j]==1) {ver[j]<-pt(auxy[j],df=nu)}
           }
       }
 }              
 if(type=="Slash") 
 {                                                    
       for (j in 1:m )
       {
           auxy[j]<- (y[j]-mu[j])/sqrt(sigmae)
           if (cens=="2")
           {
               if (cc[j]==0) {ver[j]<-dSlash(auxy[j],0,1,nu)/sqrt(sigmae)}
               if (cc[j]==1) {ver[j]<-AcumSlash(-auxy[j],0,1,nu)}
           }
          if (cens=="1")
           {
               if (cc[j]==0) {ver[j]<-dSlash(auxy[j],0,1,nu)/sqrt(sigmae)}
               if (cc[j]==1) {ver[j]<-AcumSlash(auxy[j],0,1,nu)}
           }
       }
 }  
 if(type=="NormalC") 
 {                                                    
       for (j in 1:m )
       {
           auxy[j]<- (y[j]-mu[j])/sqrt(sigmae)
           if (cens=="2")
           {
               if (cc[j]==0) {ver[j]<-dNormalC(auxy[j],0,1,nu)/sqrt(sigmae)}
               if (cc[j]==1) {ver[j]<-AcumNormalC(-auxy[j],0,1,nu)}
           }
          if (cens=="1")
           {
               if (cc[j]==0) {ver[j]<-dNormalC(auxy[j],0,1,nu)/sqrt(sigmae)}
               if (cc[j]==1) {ver[j]<-AcumNormalC(auxy[j],0,1,nu)}
           }
       }
 }              
return(ver)
}

criterios<-function(cc,y,espac=20,cadeia,type="T", cens="2", p=p){
m<-length(y)
if (type == "Normal")
{
ver<-LogVerosCens(cc,y,apply(cadeia$mu,2,mean),mean(cadeia$sigma2),1,type=type, cens=cens)
n.iter<-length(cadeia$sigma2)/espac
Loglikaux<-matrix(0,m,n.iter)
CPOaux<-matrix(0,m,n.iter)
for(k in 1:n.iter){
i<-espac*k
fss<-cadeia$mu[i,]
sigma2es<-cadeia$sigma2[i]
Loglikaux[,k]<- LogVerosCens(cc,y,fss,sigma2es,1,type=type, cens=cens)
CPOaux[,k]<-1/Loglikaux[,k]
}
Np<-length(c(p,sigma2es))
}

if (type == "T")
{
ver<-LogVerosCens(cc,y,apply(cadeia$mu,2,mean),mean(cadeia$sigma2),mean(cadeia$nu),type=type, cens=cens)
n.iter<-length(cadeia$sigma2)/espac
Loglikaux<-matrix(0,m,n.iter)
CPOaux<-matrix(0,m,n.iter)
for(k in 1:n.iter){
i<-espac*k
fss<-cadeia$mu[i,]
sigma2es<-cadeia$sigma2[i]
nuss<-cadeia$nu[i]
Loglikaux[,k]<- LogVerosCens(cc,y,fss,sigma2es,nuss,type=type, cens=cens)
CPOaux[,k]<-1/Loglikaux[,k]
}
Np<-length(c(p,sigma2es,1))
}

if (type == "Slash")
{
ver<-LogVerosCens(cc,y,apply(cadeia$mu,2,mean),mean(cadeia$sigma2),mean(cadeia$nu),type=type, cens=cens)
n.iter<-length(cadeia$sigma2)/espac
Loglikaux<-matrix(0,m,n.iter)
CPOaux<-matrix(0,m,n.iter)
for(k in 1:n.iter){
i<-espac*k
fss<-cadeia$mu[i,]
sigma2es<-cadeia$sigma2[i]
nuss<-cadeia$nu[i]
Loglikaux[,k]<- LogVerosCens(cc,y,fss,sigma2es,nuss,type=type, cens=cens)
CPOaux[,k]<-1/Loglikaux[,k]
}
Np<-length(c(p,sigma2es,1))
}

if (type == "NormalC")
{
ver<-LogVerosCens(cc,y,apply(cadeia$mu,2,mean),mean(cadeia$sigma2),c(mean(cadeia$nu),mean(cadeia$rho)),type=type, cens=cens)
n.iter<-length(cadeia$sigma2)/espac
Loglikaux<-matrix(0,m,n.iter)
CPOaux<-matrix(0,m,n.iter)
for(k in 1:n.iter){
i<-espac*k
fss<-cadeia$mu[i,]
sigma2es<-cadeia$sigma2[i]
nuss<-c(cadeia$nu[i],cadeia$rho[i])
Loglikaux[,k]<- LogVerosCens(cc,y,fss,sigma2es,nuss,type=type, cens=cens)
CPOaux[,k]<-1/Loglikaux[,k]
}
Np<-length(c(p,sigma2es,2))
}

CPO<-sum(log(1/(apply(CPOaux,1,mean))))
DIC<-2*mean(apply(-2*log(Loglikaux),2,sum))-sum(2*log(ver))
EAIC<- -2*sum(log(ver))+2*Np
EBIC<- -2*sum(log(ver))+log(m)*Np
return(list(CPO=CPO,DIC=DIC,EAIC=EAIC,EBIC=EBIC))
}

diverg<-function(cc,y,espac=50,cadeia,type, cens, p=p)
{
      m<-length(y)
      if (type == "Normal")
          {
               ver <- LogVerosCens(cc,y,apply(cadeia$mu,2,mean),mean(cadeia$sigma2),1,type=type, cens=cens)
              n.iter <- length(cadeia$sigma2)/espac
         Loglikaux <- matrix(0,m,n.iter)
            CPOaux <- matrix(0,m,n.iter)
            for(k in 1:n.iter)
            {
              i <- espac*k
            fss <- cadeia$mu[i,]
       sigma2es <- cadeia$sigma2[i]
  Loglikaux[,k] <- LogVerosCens(cc,y,fss,sigma2es,1,type=type, cens=cens)
     CPOaux[,k] <- 1/Loglikaux[,k]
            }
             Np <- p + 1
          }

if (type == "T")
{
ver<-LogVerosCens(cc,y,apply(cadeia$mu,2,mean),mean(cadeia$sigma2),mean(cadeia$nu),type=type, cens=cens)
n.iter<-length(cadeia$sigma2)/espac
Loglikaux<-matrix(0,m,n.iter)
CPOaux<-matrix(0,m,n.iter)
for(k in 1:n.iter){
i<-espac*k
fss<-cadeia$mu[i,]
sigma2es<-cadeia$sigma2[i]
nuss<-cadeia$nu[i]
Loglikaux[,k]<- LogVerosCens(cc,y,fss,sigma2es,nuss,type=type, cens=cens)
CPOaux[,k]<-1/Loglikaux[,k]
}
Np<- p + 2
}

if (type == "Slash")
{
ver<-LogVerosCens(cc,y,apply(cadeia$mu,2,mean),mean(cadeia$sigma2),mean(cadeia$nu),type=type, cens=cens)
n.iter<-length(cadeia$sigma2)/espac
Loglikaux<-matrix(0,m,n.iter)
CPOaux<-matrix(0,m,n.iter)
for(k in 1:n.iter){
i<-espac*k
fss<-cadeia$mu[i,]
sigma2es<-cadeia$sigma2[i]
nuss<-cadeia$nu[i]
Loglikaux[,k]<- LogVerosCens(cc,y,fss,sigma2es,nuss,type=type, cens=cens)
CPOaux[,k]<-1/Loglikaux[,k]
}
Np<- p + 2
}

if (type == "NormalC")
{
ver<-LogVerosCens(cc,y,apply(cadeia$mu,2,mean),mean(cadeia$sigma2),c(mean(cadeia$nu),mean(cadeia$rho)),type=type, cens=cens)
n.iter<-length(cadeia$sigma2)/espac
Loglikaux<-matrix(0,m,n.iter)
CPOaux<-matrix(0,m,n.iter)
for(k in 1:n.iter){
i<-espac*k
fss<-cadeia$mu[i,]
sigma2es<-cadeia$sigma2[i]
nuss<-c(cadeia$nu[i],cadeia$rho[i])
Loglikaux[,k]<- LogVerosCens(cc,y,fss,sigma2es,nuss,type=type, cens=cens)
CPOaux[,k]<-1/Loglikaux[,k]
}
Np<- p + 3
}

aa<- exp(log(apply(CPOaux,1,mean))+ log(Loglikaux))
IL<- apply(log(aa),1,mean)
JL<- apply((aa-1)*log(aa),1,mean)
LL<- apply(abs(aa-1),1,mean)
CHL<- apply((aa-1)^2,1,mean)

 casIL <- min(sum(IL>0.143),5)
 casJL <- min(sum(JL>0.274),5)
 casLL <- min(sum(LL>0.500),5)
casCHL <- min(sum(CHL>0.250),5)

### Diagnostico
#par(mfrow = c(2,2 ))
par(mfrow = c(1,3 ))

plot(IL, ylab="K-L divergence", xlab="Index",type="h",main="")
abline(h=0.143)
if(casIL>0)
{
identify(IL,n=casIL)
}
plot(JL, ylab="J-distance", xlab="Index",type="h",main="")
abline(h=0.27)
if(casJL>0)
{
identify(JL,n=casJL)
}
plot(LL, ylab="L_1-distance", xlab="Index",type="h",main="")
abline(h=0.5)
if(casLL>0)
{
identify(LL,n=casLL)
}
#plot(CHL, ylab="Chi divergence", xlab="Index",type="h",main="")
#abline(h=0.250)
#if(casCHL>0)
#{
#identify(CHL,n=casCHL)
#}
}

