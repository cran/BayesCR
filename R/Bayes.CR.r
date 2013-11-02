################################################################################
####     Bayesian Analysis of Normal/independent Censored Regression Model #####
### Encontrando as estimativas dos parâmetros nos Modelos Censurados para os NI
################################################################################
### Models: 
### 1- Normal               
### 2- T            ===>  Para o "nu" consideramos diferentes Aprioris
### 3- Slash        ===>  Para o "nu" consideramos diferentes Aprioris
### 5- NormalC      ===>  Para "nu1" ~ Beta(a1,b1) e "nu2" ~ Beta(a2,b2)
##############################
# algoritmo Gibbs para o caso normal, T, Slash, Normal Contaminado para os dois 
# tipos de censuras. 
# Para o caso t-Student usamos diferentes dist. a prioris para \nu (desconhecido)       #
# influence - burn    prior -  hyper  - chain
################################################################################                     

Bayes.CR <- function(cc, x,y,cens="1",type="Normal",influence="FALSE",prior=NULL,hyper=NULL,n.thin=10,burnin=100,n.iter=2000,n.chains=2,chain="TRUE")
{
	## Verify error at parameters specification

	namesx <- ('x1     ')
  if(ncol(as.matrix(x))>1){
  	for(i in 2:ncol(as.matrix(x))){namesx <- cbind(namesx, paste("x",i,"     ",sep=""))}
  }
	if(n.chains > 4) stop("The number of chains must be less than 5")
  if(ncol(as.matrix(y)) > 1) stop("Only univariate linear regression supported!")
  if(ncol(as.matrix(cc)) > 1) stop("Only univariate linear regression supported!")
	if( length(y) != nrow(as.matrix(x)) ) stop("X variable does not have the same number of lines than y")
	if( (length(x) == 0) | (length(y) == 0) ) stop("All parameters must be provided.")
	if( (type != "T") && (type != "Normal") && (type != "Slash") && (type != "NormalC")) stop("Distribution family not supported. Check documentation!")
  if( (cens != "1") && (cens != "2")) stop("Censored type not supported. 1 for left censoring and 2 for right censoring.")
  if ((burnin >= n.iter)|| burnin < 1) 
  {
  stop("Invalid burnin number")
  }
  if (!is.numeric(n.iter) || n.iter < 1) 
  {
  stop("Invalid number of iterations")
  }
  if (n.thin >= n.iter) 
  {
  stop("Invalid number of lag (n.thin) for posterior sample")
  }
  
  if(type == "T")
  {
   if(length(prior) == 0) stop("The nu prior must be provided.")
   if( (prior != "Exp") && (prior != "Jeffreys") && (prior != "Hierar") && (prior != "Unif")) stop("Prior distribution not supported. Check documentation.")
   if(prior == "Exp")
  { 
  if(length(hyper) == 0) stop("The hyper parameter must be provided.")
  if(!is.numeric(hyper)) stop("Invalid value of hyper parameter")
  if(hyper<=0) stop("The hyper parameter must be positive.")  
  }
  }
  
      out <- GibbsTruncSMN(cc,y,x,n.iter,n.thin,burnin,type=type, cens=cens, prior, hyper, n.chains=n.chains)
        p <- ncol(as.matrix(x))
    betas <- as.matrix(apply(out$beta,2,mean))
   sigma2 <- mean(out$sigma2)
    param <- rbind(betas,sigma2)
  sebetas <- as.matrix(apply(out$beta,2,sd))
 sesigma2 <- sd(out$sigma2)
       se <- rbind(sebetas,sesigma2)
      HPD <- matrix(0,nrow=p,ncol=2)
   for (i in 1:p)
   {
    HPD[i,]<- hpd(out$beta[,i],alpha=0.05) 
     
   }
  HPDS2  <- hpd(out$sigma2,alpha=0.05)    
  HPDTot <- rbind(HPD,HPDS2)
  paramT <- round(cbind(param, se, HPDTot),digits=5)

###
namespar <- colnames(x) 
    colx <- ncol(as.matrix(x))
if(length(namespar)==0)namespar <- namesx[1:colx]
#namespar <- colnames(x) 
if(n.chains>1)
{
   RBeta <- Rhat1(out$beta,n.iter,burnin,n.chains,n.thin) 
  Rsigma <- Rhat1(out$sigma2,n.iter,burnin,n.chains,n.thin)
 RhatFin <- rbind(RBeta,Rsigma)
  paramT <- round(cbind(param, se, HPDTot,RhatFin),digits=5)
dimnames(paramT) <- list(c(namespar,expression(sigma^2)),c("Mean", "Sd", " HPD(95%)","","Rhat"))
}
else{
dimnames(paramT) <- list(c(namespar,expression(sigma^2)),c("Mean", "Sd", " HPD(95%)",""))
}
 
  if( (type=="T") || (type=="Slash"))
  {
    nuFin <- mean(out$nu)
     senu <- sd(out$nu)    
   HPDnu  <- hpd(out$nu,alpha=0.05)  
    param <- rbind(betas,sigma2,nuFin)
       se <- rbind(sebetas,sesigma2,senu)
   HPDTot <- rbind(HPD,HPDS2,HPDnu)
   paramT <- round(cbind(param, se, HPDTot),digits=5)
 namespar <- colnames(x)
     colx <- ncol(as.matrix(x))
if(length(namespar)==0)namespar <- namesx[1:colx]
 
if(n.chains>1)
{
   RBeta <- Rhat1(out$beta,n.iter,burnin,n.chains,n.thin) 
  Rsigma <- Rhat1(out$sigma2,n.iter,burnin,n.chains,n.thin)
     Rnu <- Rhat1(out$nu,n.iter,burnin,n.chains,n.thin)  
 RhatFin <- rbind(RBeta,Rsigma, Rnu)
  paramT <- round(cbind(param, se, HPDTot,RhatFin),digits=5)
dimnames(paramT) <- list(c(namespar,expression(sigma^2),expression(nu)),c("Mean", "Sd", " HPD(95%)","", "Rhat"))
}
else{
dimnames(paramT) <- list(c(namespar,expression(sigma^2),expression(nu)),c("Mean", "Sd", " HPD(95%)",""))
}
  }
  if( type=="NormalC")
  {
      nu1 <- mean(out$nu)
     rho1 <- mean(out$rho)
     senu <- sd(out$nu)
    serho <- sd(out$rho)
   HPDnu  <- hpd(out$nu,alpha=0.05)
   HPDrho <- hpd(out$rho,alpha=0.05)
       se <- rbind(sebetas,sesigma2,senu,serho)
    param <- rbind(betas,sigma2,nu1,rho1)
   HPDTot <- rbind(HPD,HPDS2,HPDnu,HPDrho)
   paramT <- round(cbind(param, se, HPDTot),digits=5)
 namespar <- colnames(x)
     colx <- ncol(as.matrix(x))
if(length(namespar)==0)namespar <- namesx[1:colx]
if(n.chains>1)
{
   RBeta <- Rhat1(out$beta,n.iter,burnin,n.chains,n.thin) 
  Rsigma <- Rhat1(out$sigma2,n.iter,burnin,n.chains,n.thin)
     Rnu <- Rhat1(out$nu,n.iter,burnin,n.chains,n.thin)  
    Rrho <- Rhat1(out$rho,n.iter,burnin,n.chains,n.thin)  
 RhatFin <- rbind(RBeta,Rsigma, Rnu,Rrho)
  paramT <- round(cbind(param, se, HPDTot,RhatFin),digits=5)
dimnames(paramT) <- list(c(namespar,expression(sigma^2),expression(nu),expression(rho)),c("Mean","Sd", " HPD(95%)","","Rhat"))
}
else{
dimnames(paramT) <- list(c(namespar,expression(sigma^2),expression(nu),expression(rho)),c("Mean","Sd", " HPD(95%)",""))
}
  }
 cat('\n') 
 cat('-------------------------------------------------------------\n')
 cat('Posterior mean(Mean), standard deviation(Sd) and HPD interval\n')
 cat('-------------------------------------------------------------\n')
 print(paramT)
 cat('-------------------------------------------------------------\n')
 cat('\r \n')
      x <- as.matrix(x)
  param <- ncol(x) 
   crit <- criterios(cc,y,espac=50,cadeia=out,type=type, cens=cens, p=param) 
critFin <- c(crit$CPO, crit$DIC, crit$EAIC, crit$EBIC)
critFin <- round(t(as.matrix(critFin)),digits=3)
dimnames(critFin) <- list(c("Value"),c("LPML", "DIC", "EAIC","EBIC"))
 cat('\n') 
 cat('Model selection criteria\n')
 cat('-------------------------------------------------------------\n')
 print(critFin)
 cat('-------------------------------------------------------------\n')
 cat('\r \n')

  if(influence=="TRUE")
  {
      x <- as.matrix(x)
  param <- ncol(x) 
    Div <- diverg(cc,y,espac=50,cadeia=out,type=type, cens=cens, p=param)
  }
  if(chain=="TRUE")
  {
  return(list(beta=out$beta[,],sigma2=out$sigma2,nu=out$nu,rho=out$rho))
  }       
}




