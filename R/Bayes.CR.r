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
# criteria: TRUE ou FALSE - should compute EDC, EAIC, CPO and P_value criterias for model selection
################################################################################                     influence - burn    prior -  hyper  - chain

Bayes.CR <- function(cc, x,y,cens="1",type="Normal",criteria="FALSE",influence="FALSE",prior=NULL,hyper=NULL,n.thin=10,burnin=0.2,n.iter=25000,chain="TRUE")
{
	## Verify error at parameters specification
	if(ncol(as.matrix(y)) > 1) stop("Only univariate linear regression supported!")
  if(ncol(as.matrix(cc)) > 1) stop("Only univariate linear regression supported!")
	if( length(y) != nrow(as.matrix(x)) ) stop("X variable does not have the same number of lines than y")
	if( (length(x) == 0) | (length(y) == 0) ) stop("All parameters must be provided.")
	if( (type != "T") && (type != "Normal") && (type != "Slash") && (type != "NormalC")) stop("Distribution family not supported. Check documentation!")
  if( (cens != "1") && (cens != "2")) stop("Censored type not supported. 1 for left censoring and 2 for right censoring.")
  if (burnin < 0 || burnin > 1) 
  {
  stop("Invalid burnin proportion")
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
  
      out <- GibbsTruncSMN(cc,y,x,n.iter,n.thin,burnin,type=type, cens=cens, prior, hyper)
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
    HPD[i,]<- boa.hpd(out$beta[,i],alpha=0.05)  
   }
  HPDS2  <- boa.hpd(out$sigma2,alpha=0.05)  
  HPDTot <- rbind(HPD,HPDS2)
  paramT <- round(cbind(param, se, HPDTot),digits=5)
namespar <- colnames(x) 
 dimnames(paramT) <- list(c(namespar,expression(sigma^2)),c("Mean", "Sd", " HPD(95%)",""))
 
  if( (type=="T") || (type=="Slash"))
  {
    nuFin <- mean(out$nu)
     senu <- sd(out$nu)    
   HPDnu  <- boa.hpd(out$nu,alpha=0.05)
    param <- rbind(betas,sigma2,nuFin)
       se <- rbind(sebetas,sesigma2,senu)
   HPDTot <- rbind(HPD,HPDS2,HPDnu)
   paramT <- round(cbind(param, se, HPDTot),digits=5)
 namespar <- colnames(x)
 dimnames(paramT) <- list(c(namespar,expression(sigma^2),expression(nu)),c("Mean", "Sd", " HPD(95%)",""))
  }
  if( type=="NormalC")
  {
      nu1 <- mean(out$nu)
     rho1 <- mean(out$rho)
     senu <- sd(out$nu)
    serho <- sd(out$rho)
   HPDnu  <- boa.hpd(out$nu,alpha=0.05)
   HPDrho <- boa.hpd(out$rho,alpha=0.05)
       se <- rbind(sebetas,sesigma2,senu,serho)
    param <- rbind(betas,sigma2,nu1,rho1)
   HPDTot <- rbind(HPD,HPDS2,HPDnu,HPDrho)
   paramT <- round(cbind(param, se, HPDTot),digits=5)
 namespar <- colnames(x)
 dimnames(paramT) <- list(c(namespar,expression(sigma^2),expression(nu),expression(rho)),c("Mean","Sd", " HPD(95%)",""))
  }
 cat('\n') 
 cat('Posterior mean(Mean), standard deviation(Sd) and HPD interval\n')
 cat('--------------------------------------------------\n')
 print(paramT)
 cat('--------------------------------------------------\n')
 cat('\r')
  if(criteria=="TRUE")
  {
      x <- as.matrix(x)
  param <- ncol(x) 
   crit <- criterios(cc,y,espac=50,cadeia=out,type=type, cens=cens, p=param) 
critFin <- c(crit$CPO, crit$DIC, crit$EAIC, crit$EBIC)
critFin <- round(t(as.matrix(critFin)),digits=3)
dimnames(critFin) <- list(c("Value"),c("CPO", "DIC", "EAIC","EBIC"))
 cat('\n') 
 cat('Model selection criteria\n')
 cat('--------------------------------------------------\n')
 print(critFin)
 cat('--------------------------------------------------\n')
 cat('\r')
  }
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




